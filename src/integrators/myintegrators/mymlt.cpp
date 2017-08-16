#include "integrand.h"
#include "samplingutils.h"
#include "mutation.h"
#include "bidirmutation.h"
#include "lensperturbation.h"
#include "causticsperturbation.h"
#include "multichainperturbation.h"
#include "occlusioncache.h"
#include "parallel.h"
#include "sampleseedpaths.h"
#include "appendsensor.h"
#include "appendbsdf.h"
#include "appendemitter.h"
#include "progressreporter.h"

#include <mitsuba/render/scene.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/core/convertaether.h>
#include <mitsuba/render/unidist.h>
#include <mitsuba/render/renderqueue.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/core/fstream.h>

using namespace aether;

// Necessary for PDF caches
// TODO: make a macro for this or get rid of this totally
std::atomic<int> aether::Object::next_id{0};

MTS_NAMESPACE_BEGIN

class MyMLTIntegrator : public Integrator {
public:
    MyMLTIntegrator(const Properties &props) : Integrator(props) {
        /* Longest visualized path depth (\c -1 = infinite).
           A value of \c 1 will visualize only directly visible light sources.
           \c 2 will lead to single-bounce (direct-only) illumination, and so on. */
        m_maxDepth = props.getInteger("maxDepth", 6);
        SAssert(m_maxDepth != -1);

        /* When this flag is set to true, contributions from directly
         * visible emitters will not be included in the rendered image */
        m_hideEmitters = props.getBoolean("hideEmitters", false);

        /* Number of samples used to estimate the total luminance
           received by the scene's sensor */
        m_luminanceSamples = props.getInteger("luminanceSamples", 100000);

        /* Specifies the number of parallel work units required for
           multithreaded and network rendering. When set to <tt>-1</tt>, the
           amount will default to four times the number of cores. Note that
           every additional work unit entails a significant amount of
           communication overhead (a full-sized floating put image must be
           transmitted), hence it is important to set this value as low as
           possible, while ensuring that there are enough units to keep all
           workers busy. */
        m_workUnits = props.getInteger("workUnits", -1);

        /* This parameter can be used to specify the samples per pixel used to
           render the direct component. Should be a power of two (otherwise, it will
           be rounded to the next one). When set to zero or less, the
           direct illumination component will be hidden, which is useful
           for analyzing the component rendered by MLT. When set to -1,
           MLT will handle direct illumination as well */
        m_directSamples = props.getInteger("directSamples", -1);
        m_separateDirect = m_directSamples >= 0;

        m_specularRoughnessThreshold = props.getFloat("specularRoughnessThreshold", 0.02);

        m_misAccept = props.getBoolean("misAccept", false);

        /* Selectively enable/disable the bidirectional mutation */
        m_bidirectionalMutation = props.getBoolean("bidirectionalMutation", true);

        /* Selectively enable/disable the lens perturbation */
        m_lensPerturbation = props.getBoolean("lensPerturbation", true);

        /* Selectively enable/disable the caustic perturbation */
        m_causticPerturbation = props.getBoolean("causticPerturbation", true);

        /* Selectively enable/disable the multi-chain perturbation */
        m_multiChainPerturbation = props.getBoolean("multiChainPerturbation", true);
    }

    /// Unserialize from a binary data stream
    MyMLTIntegrator(Stream *stream, InstanceManager *manager)
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
        film->clear();
        m_running = true;
        ref<Bitmap> directBitmap = new Bitmap(Bitmap::ESpectrum, Bitmap::EFloat, film->getSize());
        directBitmap->clear();
        // Render direct-only illumination
        if (!renderDirect(scene, directBitmap, queue, job)) {
            return false;
        }

        // Compute needed number of threads
        const size_t sampleCount = sampler->getSampleCount();
        int workUnits = m_workUnits;
        if (workUnits <= 0) {
            const size_t desiredMutationsPerWorkUnit = 200000;
            const size_t cropArea  = (size_t) cropSize.x * cropSize.y;
            workUnits = ((desiredMutationsPerWorkUnit - 1) +
                (cropArea * sampleCount)) / desiredMutationsPerWorkUnit;
            Assert(workUnits >= 0 && workUnits <= (int) std::numeric_limits<int>::max());
            workUnits = std::max(workUnits, 1);
        }

        // Seeding Markov Chains
        std::vector<MarkovChainState> seeds;
        Float normalizationConstant = sampleSeedPaths(scene, sampler,
            m_luminanceSamples, workUnits, m_maxDepth, m_separateDirect, seeds);
        Log(EInfo, "Finished initialization.  Normalization constant:%f", normalizationConstant);
        SAssert(seeds.size() == workUnits);
        int nMutation = (cropSize.x * cropSize.y * sampleCount) / workUnits;
        int finishedWorkUnits = 0;
        bool unfinished = false;
        std::mutex mutex;

        // Setup MLT parameters
        // Jump sizes recommended by Eric Veach
        const Real minJump = 0.1, coveredArea = 0.05;
        // Simple heuristic for choosing a jump size: assumes that each
        // pixel on the camera subtends the same area on the sphere
        const PerspectiveCamera *camera = static_cast<const PerspectiveCamera *>(scene->getSensor());
        const Float degPerPixel = std::min(
                        camera->getXFov() / cropSize.x,
                        camera->getYFov() / cropSize.y),
                    radPerPixel = degPerPixel * M_PI / 180.0f;
        const Float causticsR1 = minJump,
                    causticsR2 = std::sqrt(coveredArea * cropSize.x*cropSize.y / M_PI); /* [Veach, p. 354] */
        const Real causticsTheta1 = radPerPixel * causticsR1;
        const Real causticsTheta2 = radPerPixel * causticsR2;
        const Real multiChainTheta1 = degToRad(0.0001f);
        const Real multiChainTheta2 = degToRad(0.1f);

        MyProgressReporter reporter(workUnits);

        ref<Bitmap> mltBitmap = new Bitmap(Bitmap::ESpectrum, Bitmap::EFloat, film->getSize());
        mltBitmap->clear();
        ref<Bitmap> combinedBitmap = new Bitmap(Bitmap::ESpectrum, Bitmap::EFloat, film->getSize());
        combinedBitmap->clear();
        ParallelFor([&](const int workUnitId) {
            ref<ImageBlock> imageBlock = new ImageBlock(Bitmap::ESpectrum,
                cropSize, film->getReconstructionFilter());
            imageBlock->setOffset(Point2i(imageBlock->getBorderSize(), imageBlock->getBorderSize()));
            imageBlock->clear();

            ref<Sampler> workUnitSampler = sampler->clone();
            UniDist uniDist(workUnitSampler.get());
            MarkovChainState currentState = seeds[workUnitId], proposalState;
            SAssert(currentState.path.Size() >= 2);
            // Time is fixed in this chain
            OcclusionCache occlusionCache(scene, currentState.path[1].Get(intersection_).time);
            // Setup mutation
            std::vector<std::shared_ptr<Mutation>> mutations;
            int minDepth = m_separateDirect ? 3 : 2;
            if (m_bidirectionalMutation) {
                mutations.push_back(std::make_shared<BidirectionalMutation>(
                    scene, uniDist, *m_raycaster, minDepth, m_maxDepth, m_specularRoughnessThreshold));
            }
            if (m_lensPerturbation) {
                mutations.push_back(std::make_shared<LensPerturbation>(
                    scene, minJump, coveredArea, m_specularRoughnessThreshold, uniDist, *m_raycaster));
            }
            if (m_causticPerturbation) {
                mutations.push_back(std::make_shared<CausticsPerturbation>(
                    causticsTheta1, causticsTheta2, m_specularRoughnessThreshold, uniDist, *m_raycaster));
            }
            if (m_multiChainPerturbation) {
                mutations.push_back(std::make_shared<MultiChainPerturbation>(
                    scene, minJump, coveredArea, multiChainTheta1, multiChainTheta2,
                    m_specularRoughnessThreshold, uniDist, *m_raycaster));
            }
            std::vector<std::shared_ptr<Mutation>> mutables;
            for (int mutationIndex = 0; mutationIndex < nMutation; mutationIndex++) {
                if (!m_running) {
                    return;
                }

                // Setup the list of mutations that can be used for this path
                mutables.clear();
                for (const auto &mutation : mutations) {
                    if (mutation->Mutable(currentState.path)) {
                        mutables.push_back(mutation);
                    }
                }

                // If there's no usable mutation, acceptable probability is zero (and we stuck here forever..)
                Assert(mutables.size() > 0);
                // Randomly choose one mutation
                const int selected = workUnitSampler->next1D() * mutables.size();
                // Actually mutate
                mutables[selected]->Mutate(scene, currentState, proposalState, occlusionCache);
                // Compute acceptance prob.
                const double a = acceptProbability(mutations, currentState, proposalState,
                    m_misAccept ? nullptr : mutables[selected].get());
                SAssert(a >= 0 && a <= 1);

                Assert(!m_separateDirect || proposalState.path.Size() > 3 || proposalState.path.Size() == 0);

                // Expectation trick
                splat(imageBlock, currentState, 1 - a);
                splat(imageBlock, proposalState, a);
                if (workUnitSampler->next1D() < a) {
                    mutables[selected]->Accepted();
                    std::swap(currentState, proposalState);
                }
            }

            std::lock_guard<std::mutex> lock(mutex);
            finishedWorkUnits++;
            Float previewFactor = Float(workUnits) / Float(finishedWorkUnits);
            mltBitmap->accumulate(imageBlock->getBitmap());
            combinedBitmap->clear();
            combinedBitmap->accumulate(mltBitmap);
            combinedBitmap->scale(previewFactor * normalizationConstant / sampler->getSampleCount());
            combinedBitmap->accumulate(directBitmap);
            film->setBitmap(combinedBitmap);
            queue->signalRefresh(job);

            reporter.Update(1);
        }, workUnits);
        TerminateWorkerThreads();
        reporter.Done();

        return !unfinished;
    }

    bool renderDirect(Scene *scene, ref<Bitmap> directBitmap, RenderQueue *queue, const RenderJob *job) {
        if (m_directSamples <= 0) {
            return true;
        }
        ref<Sensor> sensor = scene->getSensor();
        ref<Film> film = sensor->getFilm();
        const Vector2i cropSize = film->getCropSize();
        const int tileSize = scene->getBlockSize();
        const int nXTiles = (cropSize.x + tileSize - 1) / tileSize;
        const int nYTiles = (cropSize.y + tileSize - 1) / tileSize;
        bool unfinished = false;

        ref<Sampler> sampler = static_cast<Sampler *> (PluginManager::getInstance()->
            createObject(MTS_CLASS(Sampler), Properties("independent")));

        ref<ImageBlock> imageBlock = new ImageBlock(Bitmap::ESpectrum,
            cropSize, film->getReconstructionFilter());
        imageBlock->setOffset(Point2i(imageBlock->getBorderSize(), imageBlock->getBorderSize()));
        imageBlock->clear();

        std::mutex mutex;
        ParallelFor([&](const Vector2i tile) {
            ref<Sampler> clonedSampler = sampler->clone();
            const int x0 = tile.x * tileSize;
            const int x1 = std::min(x0 + tileSize, cropSize.x);
            const int y0 = tile.y * tileSize;
            const int y1 = std::min(y0 + tileSize, cropSize.y);
            for (int y = y0; y < y1; y++) {
                for (int x = x0; x < x1; x++) {
                    if (!m_running) {
                        unfinished = true;
                        return;
                    }
                    Spectrum Li(0.f);
                    for (int sampleId = 0; sampleId < m_directSamples; sampleId++) {
                        Li += directIllumination(scene, clonedSampler, x, y);
                    }
                    Li /= Float(m_directSamples);
                    std::lock_guard<std::mutex> lock(mutex);
                    // XXX: might have to fix this
                    imageBlock->put(Point2(x + 0.5, y + 0.5), &Li[0]);
                }
            }
            std::lock_guard<std::mutex> lock(mutex);
            directBitmap->copyFrom(imageBlock->getBitmap());
            film->setBitmap(directBitmap, Float(1) / sampler->getSampleCount());
            queue->signalRefresh(job);
        }, Vector2i(nXTiles, nYTiles));
        TerminateWorkerThreads();

        if (unfinished) {
            return false;
        }

        return true;
    }

    Spectrum directIllumination(const Scene *scene, Sampler *sampler, const int x, const int y) {
        UniDist uniDist(sampler);
        Spectrum Li(0.f);

        // Sample time
        // TODO: make this optional
        const Float time = sampler->next1D();
        // Sample camera subpath
        RandomSequence<Vertex> sensorSubpath;
        AppendPositionOnSensor(sensorSubpath, uniDist, scene->getSensor(), time);
        AppendDirectionFromSensor(sensorSubpath, uniDist, *m_raycaster, x, y);
        sensorSubpath.Sample();
        if (!m_hideEmitters) {
            Li += estimateBidir(scene, make_view(sensorSubpath));
        }
        if (!sensorSubpath.Back().Get(intersection_).isValid()) {
            // Hit the environment, we're done
            return Li;
        }

        auto emitterSamplingPath = sensorSubpath.Slice(0, 2);
        const ref_vector<Emitter> &emitters = scene->getEmitters();
        AppendDirectSampleEmitter(emitterSamplingPath, uniDist, emitters);
        emitterSamplingPath.Sample();

        AppendBSDF(sensorSubpath, uniDist, *m_raycaster);
        sensorSubpath.Sample();

        if (!occluded(scene,
                      to_point(emitterSamplingPath[1].Value()),
                      to_point(emitterSamplingPath[2].Value()),
                      time)) {
            const Float pdfA = emitterSamplingPath.Pdf(emitterSamplingPath);
            const Float pdfB = sensorSubpath.Pdf(emitterSamplingPath);
            Li += misWeight(pdfA, pdfB) * estimateBidir(scene, make_view(emitterSamplingPath));
        }
        
        {
            const Float pdfA = sensorSubpath.Pdf(sensorSubpath);
            const Float pdfB = emitterSamplingPath.Pdf(sensorSubpath);
            Li += misWeight(pdfA, pdfB) * estimateBidir(scene, make_view(sensorSubpath));
        }

        return Li;
    }

    void splat(ImageBlock *imageBlock, const MarkovChainState &state, const Float weight) const {
        if (weight > 0) {
            Spectrum weightedContribution = weight * Spectrum(state.contribution);
            imageBlock->put(state.screenPosition, &weightedContribution[0]);
        }
    }

    std::string toString() const {
        std::ostringstream oss;
        oss << "MyMLTIntegrator[" << endl
            << "  hideEmitters = " << m_hideEmitters << endl
            << "]";
        return oss.str();
    }

    MTS_DECLARE_CLASS()
private:
    bool m_hideEmitters;
    int m_maxDepth;
    int m_luminanceSamples;
    int m_workUnits;
    int m_directSamples;
    bool m_separateDirect;
    Float m_specularRoughnessThreshold;
    bool m_misAccept;

    bool m_bidirectionalMutation;
    bool m_lensPerturbation;
    bool m_causticPerturbation;
    bool m_multiChainPerturbation;

    std::unique_ptr<Raycast> m_raycaster;
    bool m_running;
};

MTS_IMPLEMENT_CLASS_S(MyMLTIntegrator, false, SamplingIntegrator)
MTS_EXPORT_PLUGIN(MyMLTIntegrator, "My MLT integrator");
MTS_NAMESPACE_END
