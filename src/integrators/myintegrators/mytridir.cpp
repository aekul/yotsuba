#include <mitsuba/render/scene.h>
#include <mitsuba/core/convertaether.h>
#include <mitsuba/render/unidist.h>
#include <mitsuba/render/renderqueue.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/core/fstream.h>
#include <mutex>
#include "integrand.h"
#include "samplingutils.h"
#include "appendsensor.h"
#include "appendbsdf.h"
#include "appendemitter.h"
#include "appendkeyhole.h"
#include "parallel.h"
#include "occlusioncache.h"
#include "progressreporter.h"

using namespace aether;

// Necessary for PDF caches
// TODO: make a macro for this or get rid of this totally
std::atomic<int> aether::Object::next_id{0};

MTS_NAMESPACE_BEGIN

class MyTriDirIntegrator : public Integrator {
public:
	MyTriDirIntegrator(const Properties &props) : Integrator(props) {
		/* Depth to begin using russian roulette */
		m_rrDepth = props.getInteger("rrDepth", 5);

		/* Longest visualized path depth (\c -1 = infinite).
		   A value of \c 1 will visualize only directly visible light sources.
		   \c 2 will lead to single-bounce (direct-only) illumination, and so on. */
		m_maxDepth = props.getInteger("maxDepth", 6);
		SAssert(m_maxDepth != -1);

		/* When this flag is set to true, contributions from directly
		 * visible emitters will not be included in the rendered image */
		m_hideEmitters = props.getBoolean("hideEmitters", false);
	}

	/// Unserialize from a binary data stream
	MyTriDirIntegrator(Stream *stream, InstanceManager *manager)
	 : Integrator(stream, manager) {
		m_hideEmitters = stream->readBool();
		configure();
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		Integrator::serialize(stream, manager);
		stream->writeBool(m_hideEmitters);
	}

	void configure() {
		Integrator::configure();
	}

	bool preprocess(const Scene *scene, RenderQueue *queue, const RenderJob *job,
			int sceneResID, int sensorResID, int samplerResID) {
		Integrator::preprocess(scene, queue, job, sceneResID, sensorResID, samplerResID);
		m_raycaster = std::make_unique<Raycast>(scene);
		return true;
	}

	void cancel() {
		m_running = false;
	}

	void configureSampler(const Scene *scene, Sampler *sampler) {
		/* Prepare the sampler for tile-based rendering */
		sampler->setFilmResolution(scene->getFilm()->getCropSize(), true);
	}

	bool render(Scene *scene, RenderQueue *queue, const RenderJob *job,
			int sceneResID, int sensorResID, int samplerResID) {
		ref<Scheduler> sched = Scheduler::getInstance();
		ref<Sensor> sensor = scene->getSensor();
		ref<Film> film = sensor->getFilm();
		ref<Sampler> sampler = scene->getSampler();
		size_t nCores = sched->getCoreCount();
		Log(EInfo, "Starting render job (%ix%i, " SIZE_T_FMT " %s, " SSE_STR ") ..",
			film->getCropSize().x, film->getCropSize().y,
			nCores, nCores == 1 ? "core" : "cores");

		// TODO: take cropOffset into account
		//Point2i cropOffset = film->getCropOffset();
		Vector2i cropSize = film->getCropSize();
		ref<ImageBlock> imageBlock = new ImageBlock(Bitmap::ESpectrum,
			cropSize, film->getReconstructionFilter());
		imageBlock->setOffset(Point2i(imageBlock->getBorderSize(), imageBlock->getBorderSize()));
		ref<Bitmap> bitmap = new Bitmap(Bitmap::ESpectrum, Bitmap::EFloat, film->getSize());
		imageBlock->clear();
		bitmap->clear();
		film->clear();
		m_running = true;

	    const int tileSize = scene->getBlockSize();
	    const int nXTiles = (cropSize.x + tileSize - 1) / tileSize;
	    const int nYTiles = (cropSize.y + tileSize - 1) / tileSize;
        MyProgressReporter reporter(nXTiles * nYTiles);
	    bool unfinished = false;
	    std::mutex mutex;

		ParallelFor([&](const Vector2i tile) {
			ref<Sampler> clonedSampler = sampler->clone();
        	const int x0 = tile.x * tileSize;
        	const int x1 = std::min(x0 + tileSize, cropSize.x);
        	const int y0 = tile.y * tileSize;
        	const int y1 = std::min(y0 + tileSize, cropSize.y);
			// For each pixel
			for (int y = y0; y < y1; y++) {
				for (int x = x0; x < x1; x++) {
					if (!m_running) {
						unfinished = true;
						return;
					}
					// For each sample
					for (int sampleId = 0; sampleId < sampler->getSampleCount(); sampleId++) {
						render(scene, clonedSampler, imageBlock.get(), x, y, mutex);
					}
				}
			}
			std::lock_guard<std::mutex> lock(mutex);
			bitmap->copyFrom(imageBlock->getBitmap());
			film->setBitmap(bitmap, Float(1) / sampler->getSampleCount());
			queue->signalRefresh(job);
            reporter.Update(1);
		}, Vector2i(nXTiles, nYTiles));
		TerminateWorkerThreads();

		if (unfinished) {
			return false;
		}

        reporter.Done();
		return true;
	}

	struct SplatElement {
		Point2 position;
		Spectrum contribution;
		int pathLength;
		int sensorSubpathIndex;
	};

	void render(const Scene *scene, Sampler *sampler, ImageBlock *imageBlock, const int x, const int y,
				std::mutex &mutex) {
		UniDist uniDist(sampler);

        // TODO: move triangle co-ordinates to scene file
        // cbox_keyhole
        //aether::Vector3 tri1_v0{-0.173878, 1.990000, -0.189689};
        //aether::Vector3 tri1_v1{ 0.153878, 1.990000, 0.139690};
        //aether::Vector3 tri1_v2{-0.173878, 1.990000,  0.139690};

        //aether::Vector3 tri2_v0{-0.173878, 1.990000, -0.189689};
        //aether::Vector3 tri2_v1{ 0.153878, 1.990000, -0.189689};
        //aether::Vector3 tri2_v2{ 0.153878, 1.990000, 0.139690};


        // veach_door
        aether::Vector3 tri1_v0{121.9446, 0, -73.0101};
        aether::Vector3 tri1_v1{119.1694, 121.5088, -56.649};
        aether::Vector3 tri1_v2{119.1694, 0.5088, -56.649};

        aether::Vector3 tri2_v0{121.9446, 0, -73.0101};
        aether::Vector3 tri2_v1{121.9446, 123.3378, -73.0101};
        aether::Vector3 tri2_v2{119.1694, 121.5088, -56.649};

        // Top of opening
        aether::Vector3 tri3_v0{121.9446, 123.3378, -73.0101};
        aether::Vector3 tri3_v1{62.1593, 121.5088, -74.4293};
        aether::Vector3 tri3_v2{119.1694, 121.5088, -56.649};

        // Back of door
        aether::Vector3 tri4_v0{121.9446, 123.3378, -77.0547};
        aether::Vector3 tri4_v1{58.7111, 123.3378, -77.0547};
        aether::Vector3 tri4_v2{58.7111, 0, -77.0547};

        aether::Vector3 tri5_v0{121.9446, 123.3378, -77.0547};
        aether::Vector3 tri5_v1{58.7111, 0, -77.0547 };
        aether::Vector3 tri5_v2{121.9446, 0, -77.0547};

        std::vector<std::array<aether::Vector3, 3>> keyhole_tris;
        keyhole_tris.push_back({{tri1_v0, tri1_v1, tri1_v2}});
        keyhole_tris.push_back({{tri2_v0, tri2_v1, tri2_v2}});
        keyhole_tris.push_back({{tri3_v0, tri3_v1, tri3_v2}});
        keyhole_tris.push_back({{tri4_v0, tri4_v1, tri4_v2}});
        keyhole_tris.push_back({{tri5_v0, tri5_v1, tri5_v2}});

        RandomSequence<Vertex> keyholePath;
        AppendKeyhole(keyholePath, uniDist, *m_raycaster, keyhole_tris);
        keyholePath.Sample();

        const ref_vector<Emitter> &emitters = scene->getEmitters();

        // Sample time
        // TODO: make this optional
        const Float time = sampler->next1D();
        OcclusionCache occlusionCache(scene, time);

        // Sample camera subpath
        RandomSequence<Vertex> sensorSubpath;
        AppendPositionOnSensor(sensorSubpath, uniDist, scene->getSensor(), time);
        AppendDirectionFromSensor(sensorSubpath, uniDist, *m_raycaster, x, y);
        sensorSubpath.Sample();

        // Loop until reaching specified maximum depth
        // e.g. if maxDepth == 2, we don't need sensorSubpath.Size() > 3
        for(;sensorSubpath.Size() <= m_maxDepth;) {
            AppendBSDF(sensorSubpath, uniDist, *m_raycaster);
            sensorSubpath.Sample();
        }

        // Sample emitter subpath
        RandomSequence<Vertex> emitterSubpath;
        AppendPositionOnEmitter(emitterSubpath, uniDist, emitters, time);
        AppendDirectionFromEmitter(emitterSubpath, uniDist, *m_raycaster);
        emitterSubpath.Sample();
        // Loop until reaching specified maximum depth
        // e.g. if maxDepth == 2, we don't need sensorSubpath.Size() > 3
		for(;emitterSubpath.Size() <= m_maxDepth;) {
			AppendBSDF(emitterSubpath, uniDist, *m_raycaster);
			emitterSubpath.Sample();
		}

		// Combine
        std::vector<SplatElement> splats;
        for (int pathLength = 2; pathLength <= m_maxDepth; pathLength++) {
            std::vector<RandomSequence<Vertex>> paths;
            for (int sensorSubpathSize = 1; sensorSubpathSize <= pathLength + 1; sensorSubpathSize++) {
                const int emitterSubpathSize = pathLength + 1 - sensorSubpathSize;
                auto sensorSubpathSlice = sensorSubpath.Slice(0, sensorSubpathSize);
                RandomSequence<Vertex> path;

                // Tri-directional path
                // - Shorten the sensor and emitter subpaths by 1
                // - Replace the removed vertices with the two keyhole vertices
                if (sensorSubpathSize > 1 && emitterSubpathSize >= 1) {
                    auto sensorTriDirSubpathSlice = sensorSubpath.Slice(0, sensorSubpathSize - 1);
                    auto emitterTriDirSubpathSlice = reverse_(emitterSubpath.Slice(0, emitterSubpathSize - 1));
                    paths.push_back(sensorTriDirSubpathSlice.Concat(keyholePath).Concat(emitterTriDirSubpathSlice));
                }

                // - Shorten the sensor subpaths by 2
                if (sensorSubpathSize > 2) {
                    auto sensorTriDirSubpathSlice = sensorSubpath.Slice(0, sensorSubpathSize - 2);
                    auto emitterTriDirSubpathSlice = reverse_(emitterSubpath.Slice(0, emitterSubpathSize));
                    paths.push_back(sensorTriDirSubpathSlice.Concat(keyholePath).Concat(emitterTriDirSubpathSlice));
                }
                // - Shorten the emitter subpaths by 2
                if (emitterSubpathSize >= 2) {
                    auto sensorTriDirSubpathSlice = sensorSubpath.Slice(0, sensorSubpathSize);
                    auto emitterTriDirSubpathSlice = reverse_(emitterSubpath.Slice(0, emitterSubpathSize - 2));
                    paths.push_back(sensorTriDirSubpathSlice.Concat(keyholePath).Concat(emitterTriDirSubpathSlice));
                }


                // Bi-directional paths
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
            estimate(scene, imageBlock, paths, time, splats, occlusionCache);
        }

		std::lock_guard<std::mutex> lock(mutex);

		if (!m_hideEmitters) {
			auto path = sensorSubpath.Slice(0, 2);
			if (path.AllValid()) {
				Spectrum contribution = estimateBidir(scene, make_view(path));
				Point2 position;
				if (project(scene, sensorSubpath, position)) {
					imageBlock->put(position, &contribution[0]);
				}
			}
		}

		for (const auto &splatElement : splats) {
			Point2 position = splatElement.position;
			Spectrum contribution = splatElement.contribution;
			imageBlock->put(position, &contribution[0]);
		}
	}

	void estimate(const Scene *scene,
				  ImageBlock *imageBlock,
				  const std::vector<RandomSequence<Vertex>> &paths,
				  const Float time,
				  std::vector<SplatElement> &splats,
				  OcclusionCache& occlusionCache) const {
		std::vector<double> pdfs;
		for (size_t i = 0; i < paths.size(); i++) {
			const RandomSequence<Vertex> &path = paths[i];
			if (!path.AllValid()) {
				continue;
			}
			Point2 position;
			if (!project(scene, path, position)) {
				continue;
			}
      		// TODO: improve this; check only the connection edges
      		bool occ = false;
			for (size_t j = 0, N = path.Size() - 1; j < N; j++) {
        		if (occlusionCache.Query(std::make_pair(path[j].Value(), path[j + 1].Value()))) {
          			occ = true;
          			break;
        		}
			}
      		if (occ) {
        		continue;
      		}
			Spectrum contribution = estimateBidir(scene, make_view(path));
			if (contribution.isZero()) {
				continue;
			}

			pdfs.clear();
			for (size_t j = 0; j < paths.size(); j++) {
				double pdf = paths[j].Pdf(path);
				pdfs.push_back(pdf);
			}
			double weight = misWeight(i, pdfs);
			Spectrum weightedContribution = weight * contribution;
			int pathLength = path.Size() - 1;
			splats.push_back(SplatElement{position, weightedContribution, pathLength, (int)i});
		}
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "MyTriDirIntegrator[" << endl
			<< "  hideEmitters = " << m_hideEmitters << endl
			<< "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
private:
	bool m_hideEmitters;
	int m_rrDepth;
	int m_maxDepth;

	std::unique_ptr<Raycast> m_raycaster;
	bool m_running;
};

MTS_IMPLEMENT_CLASS_S(MyTriDirIntegrator, false, SamplingIntegrator)
MTS_EXPORT_PLUGIN(MyTriDirIntegrator, "My tridir integrator");
MTS_NAMESPACE_END
