#include <mitsuba/render/scene.h>
#include <mitsuba/core/convertaether.h>
#include <mitsuba/render/unidist.h>
#include "integrand.h"
#include "samplingutils.h"
#include "appendsensor.h"
#include "appendbsdf.h"
#include "appendemitter.h"

using namespace aether;

// Necessary for PDF caches
// TODO: make a macro for this or get rid of this totally
std::atomic<int> aether::Object::next_id{0};

MTS_NAMESPACE_BEGIN

class MyPathIntegrator : public SamplingIntegrator {
public:
	MyPathIntegrator(const Properties &props) : SamplingIntegrator(props) {
		/* Depth to begin using russian roulette */
		m_rrDepth = props.getInteger("rrDepth", 5);

		/* Longest visualized path depth (\c -1 = infinite).
		   A value of \c 1 will visualize only directly visible light sources.
		   \c 2 will lead to single-bounce (direct-only) illumination, and so on. */
		m_maxDepth = props.getInteger("maxDepth", -1);

		/* When this flag is set to true, contributions from directly
		 * visible emitters will not be included in the rendered image */
		m_hideEmitters = props.getBoolean("hideEmitters", false);
	}

	/// Unserialize from a binary data stream
	MyPathIntegrator(Stream *stream, InstanceManager *manager)
	 : SamplingIntegrator(stream, manager) {
		m_hideEmitters = stream->readBool();
		configure();
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		SamplingIntegrator::serialize(stream, manager);
		stream->writeBool(m_hideEmitters);
	}

	void configure() {
		SamplingIntegrator::configure();
	}

	void configureSampler(const Scene *scene, Sampler *sampler) {
		SamplingIntegrator::configureSampler(scene, sampler);
	}

	bool preprocess(const Scene *scene, RenderQueue *queue, const RenderJob *job,
			int sceneResID, int sensorResID, int samplerResID) {
		SamplingIntegrator::preprocess(scene, queue, job, sceneResID, sensorResID, samplerResID);
		m_raycaster = std::make_unique<Raycast>(scene);
		return true;
	}

	Spectrum Li(const RayDifferential &ray, RadianceQueryRecord &rRec) const {
		UniDist uniDist(rRec.sampler);
		const Scene *scene = rRec.scene;
		Spectrum Li(0.f);

		// Sample camera subpath
		RandomSequence<Vertex> sensorSubpath;
		AppendPositionOnSensor(sensorSubpath, ray);
		AppendDirectionFromSensor(sensorSubpath, *m_raycaster, ray);
		sensorSubpath.Sample();
		if (!m_hideEmitters) {
			Li += estimate(scene, make_view(sensorSubpath));
		}

		for(;;) {
			if (!sensorSubpath.Back().Valid()) {
				// Invalidated vertex, probably russian roulette termination
				sensorSubpath.RemoveBack();
				break;
			}
			if (!sensorSubpath.Back().Get(intersection_).isValid()) {
				// Hit the environment, we're done
				break;
			}
			if (m_maxDepth != -1 && sensorSubpath.Size() >= m_maxDepth + 1) {
				// Reached specified maximum depth
				// e.g. if maxDepth == 2, we don't need sensorSubpath.Size() > 3
				break;
			}
			AppendBSDF(sensorSubpath, uniDist, *m_raycaster, m_rrDepth);
			sensorSubpath.Sample();
		}

		const ref_vector<Emitter> &emitters = rRec.scene->getEmitters();
		for (int pathLength = 2; pathLength < sensorSubpath.Size(); pathLength++) {
			// Sample emitter subpath
			RandomSequence<Vertex> emitterSamplingPath = sensorSubpath.Slice(0, pathLength);
			if (!emitterSamplingPath.Back().Get(intersection_).isValid()) {
				continue;
			}
			AppendDirectSampleEmitter(emitterSamplingPath, uniDist, emitters);
			emitterSamplingPath.Sample();
			auto bsdfSamplingPath = sensorSubpath.Slice(0, pathLength + 1);

			// Direct light source sampling
			if (!occluded(rRec.scene,
					to_point(emitterSamplingPath[emitterSamplingPath.Size() - 2].Value()),
					to_point(emitterSamplingPath[emitterSamplingPath.Size() - 1].Value()),
					ray.time)) {
				const Float pdfA = emitterSamplingPath.Pdf(emitterSamplingPath);
				const Float pdfB = bsdfSamplingPath.Pdf(emitterSamplingPath);
				Li += misWeight(pdfA, pdfB) * estimate(scene, make_view(emitterSamplingPath));
			}
			
			// BSDF importance sampling
			{
				const Float pdfA = bsdfSamplingPath.Pdf(bsdfSamplingPath);
				const Float pdfB = emitterSamplingPath.Pdf(bsdfSamplingPath);
				Li += misWeight(pdfA, pdfB) * estimate(scene, make_view(bsdfSamplingPath));
			}
		}

		return Li;
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "MyPathIntegrator[" << endl
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
};

MTS_IMPLEMENT_CLASS_S(MyPathIntegrator, false, SamplingIntegrator)
MTS_EXPORT_PLUGIN(MyPathIntegrator, "My path integrator");
MTS_NAMESPACE_END
