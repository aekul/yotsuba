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

class MyMIDirectIntegrator : public SamplingIntegrator {
public:
	MyMIDirectIntegrator(const Properties &props) : SamplingIntegrator(props) {
		/* When this flag is set to true, contributions from directly
		 * visible emitters will not be included in the rendered image */
		m_hideEmitters = props.getBoolean("hideEmitters", false);
	}

	/// Unserialize from a binary data stream
	MyMIDirectIntegrator(Stream *stream, InstanceManager *manager)
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
		// Call .Sample() to actually sample the path
		sensorSubpath.Sample();
		if (!m_hideEmitters) {
			Li += estimate(scene, make_view(sensorSubpath));
		}
		if (!sensorSubpath.Back().Get(intersection_).isValid()) {
			// Hit the environment, we're done
			return Li;
		}

		// Direct light source sampling
		auto emitterSamplingPath = sensorSubpath.Slice(0, 2);
		const ref_vector<Emitter> &emitters = rRec.scene->getEmitters();
		AppendDirectSampleEmitter(emitterSamplingPath, uniDist, emitters);
		emitterSamplingPath.Sample();

		// BSDF importance sampling
		AppendBSDF(sensorSubpath, uniDist, *m_raycaster);
		sensorSubpath.Sample();

		if (!occluded(rRec.scene,
					  to_point(emitterSamplingPath[1].Value()),
					  to_point(emitterSamplingPath[2].Value()), ray.time)) {
			const Float pdfA = emitterSamplingPath.Pdf(emitterSamplingPath);
			const Float pdfB = sensorSubpath.Pdf(emitterSamplingPath);
			Li += misWeight(pdfA, pdfB) * estimate(scene, make_view(emitterSamplingPath));
		}
	
		{
			const Float pdfA = sensorSubpath.Pdf(sensorSubpath);
			const Float pdfB = emitterSamplingPath.Pdf(sensorSubpath);
			Li += misWeight(pdfA, pdfB) * estimate(scene, make_view(sensorSubpath));
		}

		return Li;
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "MyMIDirectIntegrator[" << endl
			<< "  hideEmitters = " << m_hideEmitters << endl
			<< "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
private:
	bool m_hideEmitters;
	std::unique_ptr<Raycast> m_raycaster;
};

MTS_IMPLEMENT_CLASS_S(MyMIDirectIntegrator, false, SamplingIntegrator)
MTS_EXPORT_PLUGIN(MyMIDirectIntegrator, "My direct illumination integrator");
MTS_NAMESPACE_END
