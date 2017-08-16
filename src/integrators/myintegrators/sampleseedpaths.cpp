#include "sampleseedpaths.h"
#include "markovchainstate.h"
#include "integrand.h"
#include "samplingutils.h"
#include "appendsensor.h"
#include "appendbsdf.h"
#include "appendemitter.h"
#include "parallel.h"
#include <mitsuba/render/scene.h>
#include <aether/RandomSequence.h>
#include <vector>
#include <mutex>

MTS_NAMESPACE_BEGIN

static Float estimate(const Scene *scene,
					  const std::vector<RandomSequence<Vertex>> &paths,
					  const Float time,
					  std::vector<MarkovChainState> &seedPathsPool,
					  std::vector<Float> &seedPathsImportance) {
	Float totalImportance(0.f);
	std::vector<double> pdfs;		
	for (size_t i = 0; i < paths.size(); i++) {
		const RandomSequence<Vertex> &path = paths[i];
		if (!path.AllValid()) {
			continue;
		}
		Point2 screenPosition;
		if (!project(scene, path, screenPosition)) {
			continue;
		}
		if (i < paths.size() - 1) {
			if (occluded(scene, to_point(path[i].Value()),
								to_point(path[i + 1].Value()), time)) {
				continue;
			}
		}
		const TSpectrum<double, SPECTRUM_SAMPLES> f = evalBidir(scene, make_view(path));
		if (f.isZero() || !f.isValid()) {
			continue;
		}
		const double p = path.Pdf();
		if (p <= 0. || !std::isfinite(p)) {
			// TODO: issue warning here
			continue;
		}

		pdfs.clear();
		for (size_t j = 0; j < paths.size(); j++) {
			pdfs.push_back(paths[j].Pdf(path));
		}
		const double weight = misWeight(i, pdfs);
		if (weight <= 0.) {
			continue;
		}
		const double imp = importance(f);
		if (imp <= 0) {
			continue;
		}
		const double weightedImportance = (weight / p) * imp;
		SAssert(std::isfinite(weightedImportance));
		totalImportance += weightedImportance;
		seedPathsPool.push_back(MarkovChainState{f / imp, imp, screenPosition, path});
		seedPathsImportance.push_back(weightedImportance);
	}
	return totalImportance;
}

Float sampleSeedPaths(const Scene *scene, Sampler *sampler, 
		const int luminanceSamples, const int workUnits,
		const int maxDepth, const bool separateDirect,
		std::vector<MarkovChainState> &seeds) {
	std::vector<MarkovChainState> seedPathsPool;
	std::vector<Float> seedPathsImportance;
	Float importanceSum(0);
	UniDist uniDist(sampler);
	Raycast raycaster(scene);

	const int64_t numSamplesPerThread = luminanceSamples / NumSystemCores();
	const int64_t threadsNeedExtraSamples = numSamplesPerThread % NumSystemCores();
	std::mutex mutex;
	ParallelFor([&](const int workerId) {
		const int64_t numSamplesThisThread =
			numSamplesPerThread + ((workerId < threadsNeedExtraSamples) ? 1 : 0);
		std::vector<MarkovChainState> localSeedPathsPool;
		std::vector<Float> localSeedPathsImportance;
		for (int sampleIndex = 0; sampleIndex < numSamplesThisThread; sampleIndex++) {
			// Sample time
			// TODO: make this optional
			const Float time = sampler->next1D();
			// Sample camera subpath
			RandomSequence<Vertex> sensorSubpath;
			AppendPositionOnSensor(sensorSubpath, uniDist, scene->getSensor(), time);
			AppendDirectionFromSensor(sensorSubpath, uniDist, raycaster);
			sensorSubpath.Sample();

			// Loop until reaching specified maximum depth
			// e.g. if maxDepth == 2, we don't need sensorSubpath.Size() > 3
			for(;sensorSubpath.Size() <= maxDepth;) {
				AppendBSDF(sensorSubpath, uniDist, raycaster);
				sensorSubpath.Sample();
			}

			// Sample emitter subpath
			RandomSequence<Vertex> emitterSubpath;
			AppendPositionOnEmitter(emitterSubpath, uniDist, scene->getEmitters(), time);
			AppendDirectionFromEmitter(emitterSubpath, uniDist, raycaster);
			emitterSubpath.Sample();
			// Loop until reaching specified maximum depth
			// e.g. if maxDepth == 2, we don't need sensorSubpath.Size() > 3
			for(;emitterSubpath.Size() <= maxDepth;) {
				AppendBSDF(emitterSubpath, uniDist, raycaster);
				emitterSubpath.Sample();
			}

			for (int pathLength = separateDirect ? 3 : 2; pathLength <= maxDepth; pathLength++) {
				std::vector<RandomSequence<Vertex>> paths;
				for (int sensorSubpathSize = 1; sensorSubpathSize <= pathLength + 1; sensorSubpathSize++) {
					const int emitterSubpathSize = pathLength + 1 - sensorSubpathSize;
					auto sensorSubpathSlice = sensorSubpath.Slice(0, sensorSubpathSize);
					RandomSequence<Vertex> path;
					if (emitterSubpathSize != 1) {
						auto emitterSubpathSlice = reverse_(emitterSubpath.Slice(0, emitterSubpathSize));
						path = sensorSubpathSlice.Concat(emitterSubpathSlice);
					} else {
						// Special case: we want to do specialized direct importance sampling here
						AppendDirectSampleEmitter(sensorSubpathSlice, uniDist, scene->getEmitters());
						sensorSubpathSlice.Sample();
						path = sensorSubpathSlice;
					}
					SAssert(path.Size() == pathLength + 1);
					paths.push_back(path);
				}
				importanceSum += estimate(scene, paths, time, localSeedPathsPool, localSeedPathsImportance);
			}
		}

		std::lock_guard<std::mutex> lock(mutex);
		// merge seed pool
		seedPathsPool.insert(seedPathsPool.end(), localSeedPathsPool.begin(), localSeedPathsPool.end());
		seedPathsImportance.insert(
				seedPathsImportance.end(), localSeedPathsImportance.begin(), localSeedPathsImportance.end());
	}, NumSystemCores());

	if (seedPathsImportance.size() < workUnits || importanceSum <= Float(0)) {
		SLog(EError, "MLT initialization failed.");
	}

	SAssert(seedPathsPool.size() == seedPathsImportance.size());
	SLog(EInfo, "Taking %i seeds from %i...", workUnits, seedPathsPool.size());

	DiscreteDistribution seedPDF(seedPathsImportance.size());
	for (size_t i=0; i < seedPathsImportance.size(); ++i) {
		seedPDF.append(seedPathsImportance[i]);
	}
	seedPDF.normalize();

	seeds.clear();
	seeds.reserve(workUnits);
	for (size_t i=0; i < workUnits; ++i) {
		seeds.push_back(seedPathsPool.at(seedPDF.sample(sampler->next1D())));
	}

	return importanceSum / luminanceSamples;
}

MTS_NAMESPACE_END
