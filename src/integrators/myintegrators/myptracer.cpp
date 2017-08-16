#include <mitsuba/render/scene.h>
#include <mitsuba/core/convertaether.h>
#include <mitsuba/render/unidist.h>
#include <mitsuba/render/renderqueue.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/core/fstream.h>
#include "integrand.h"
#include "samplingutils.h"

using namespace aether;

// Necessary for PDF caches
// TODO: make a macro for this or get rid of this totally
std::atomic<int> aether::Object::next_id{0};

MTS_NAMESPACE_BEGIN

class MyAdjointParticleTracer : public Integrator {
public:
    MyAdjointParticleTracer(const Properties &props) : Integrator(props) {
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
    MyAdjointParticleTracer(Stream *stream, InstanceManager *manager)
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
        m_sampleBSDF = Node<sample_bsdf_rr_t>(m_rrDepth);
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

        // For each pixel
        for (int y = 0; y < cropSize.y; y++) {
            for (int x = 0; x < cropSize.x; x++) {
                if (!m_running) {
                    return false;
                }
                // For each sample
                for (int sampleId = 0; sampleId < sampler->getSampleCount(); sampleId++) {
                    render(scene, sampler, imageBlock.get(), x, y);
                }
            }
            bitmap->copyFrom(imageBlock->getBitmap(), imageBlock->getOffset());
            film->setBitmap(bitmap, Float(1) / sampler->getSampleCount());
            queue->signalRefresh(job);
        }

        return true;
    }

    void render(const Scene *scene, Sampler *sampler, ImageBlock *imageBlock, const int x, const int y) {
        UniDist uniDist(sampler);

        // Sample time
        // TODO: make this optional
        const Float time = sampler->next1D();

        // Sample camera subpath
        RandomSequence<Vertex> sensorSubpath;
        sensorSubpath.Append(m_samplePositionOnCamera, uniDist, scene->getSensor(), time);
        sensorSubpath.Append(m_sampleDirectionFromCamera, uniDist, *m_raycaster, x, y);
        sensorSubpath.Sample();

        if (!m_hideEmitters) {
            Spectrum contribution = estimateBidir(scene, make_view(sensorSubpath));
            Point2 position;
            if (project(scene, sensorSubpath, position)) {
                imageBlock->put(position, &contribution[0]);
            }
        }

        // Loop until reaching specified maximum depth
        // e.g. if maxDepth == 2, we don't need sensorSubpath.Size() > 3
        for(;sensorSubpath.Size() <= m_maxDepth;) {
            sensorSubpath.Append(m_sampleBSDF, uniDist, *m_raycaster);
            sensorSubpath.Sample();
        }

        // Sample emitter subpath
        const ref_vector<Emitter> &emitters = scene->getEmitters();
        RandomSequence<Vertex> emitterSubpath;
        emitterSubpath.Append(m_samplePositionOnEmitter, uniDist, emitters, time);
        emitterSubpath.Append(m_sampleDirectionFromEmitter, uniDist, *m_raycaster);
        emitterSubpath.Sample();
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
                // e.g. if maxDepth == 2, we don't need emitterSubpath.Size() > 3
                break;
            }
            emitterSubpath.Append(m_sampleBSDF, uniDist, *m_raycaster);
            emitterSubpath.Sample();
        }

        auto sensorVertex = sensorSubpath.Slice(0, 1);

        for (int pathLength = 2; pathLength <= m_maxDepth; pathLength++) {
            auto emitterSubpathSlice = reverse_(emitterSubpath.Slice(0, pathLength));
            auto path = sensorVertex.Concat(emitterSubpathSlice);
            estimate(scene, imageBlock, path, time);
        }
    }

    void estimate(const Scene *scene,
                  ImageBlock *imageBlock,
                  const RandomSequence<Vertex> &path,
                  const Float time) const {
        if (!path.AllValid()) {
            return;
        }
        Point2 position;
        if (!project(scene, path, position)) {
            return;
        }
        if (occluded(scene, to_point(path[0].Value()),
                            to_point(path[1].Value()), time)) {
            return;
        }
        
        Spectrum contribution = estimateBidir(scene, make_view(path));
        if (contribution.isZero()) {
            return;
        }

        imageBlock->put(position, &contribution[0]);
    }

    std::string toString() const {
        std::ostringstream oss;
        oss << "MyAdjointParticleTracer[" << endl
            << "  hideEmitters = " << m_hideEmitters << endl
            << "]";
        return oss.str();
    }

    MTS_DECLARE_CLASS()
private:
    bool m_hideEmitters;
    int m_rrDepth;
    int m_maxDepth;

    // TODO: get rid of mutable
    mutable Node<sample_position_on_camera_t> m_samplePositionOnCamera;
    mutable Node<sample_direction_from_camera_t> m_sampleDirectionFromCamera;
    mutable Node<sample_bsdf_rr_t> m_sampleBSDF;
    mutable Node<sample_position_on_emitter_t> m_samplePositionOnEmitter;
    mutable Node<sample_direction_from_emitter_t> m_sampleDirectionFromEmitter;
    std::unique_ptr<Raycast> m_raycaster;
    bool m_running;
};

MTS_IMPLEMENT_CLASS_S(MyAdjointParticleTracer, false, SamplingIntegrator)
MTS_EXPORT_PLUGIN(MyAdjointParticleTracer, "My adjoint particle tracer");
MTS_NAMESPACE_END
