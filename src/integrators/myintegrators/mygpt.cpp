#include "integrand.h"
#include "samplingutils.h"
#include "appendsensor.h"
#include "appendbsdf.h"
#include "appendemitter.h"
#include "parallel.h"
#include "occlusioncache.h"
#include "poissonsolver.h"
#include "classification.h"
#include "intersect.h"
#include "perturbdirection.h"
#include "progressreporter.h"

#include <mitsuba/render/scene.h>
#include <mitsuba/core/convertaether.h>
#include <mitsuba/render/unidist.h>
#include <mitsuba/render/renderqueue.h>
#include <mitsuba/render/aethersensors_impl.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/core/fstream.h>
#include <mutex>

using namespace aether;

// Necessary for PDF caches
// TODO: make a macro for this or get rid of this totally
std::atomic<int> aether::Object::next_id{0};

namespace aether {

    // version of ConditionalPdf that doesn't allocate extra memory
    template <typename Mutation, typename T>
    Real ConditionalPdf(Node<Mutation>& mutation,
                        const RandomSequence<T>& output_seq,
                        const RandomSequence<T>& input_seq,
                        RandomSequence<T> &input_buffer) {
        if (!output_seq.AllValid() || output_seq.IsEmpty()) {
            return Real(0);
        }

        std::vector<int> shared_rv;
        shared_rv.reserve(input_seq.Size());
        // For each r.v. in input_seq, check if there is any r.v. in output_seq that shares the same strategy
        // If so, the output r.v. is copied from input r.v., and this results in a Dirac conditional pdf
        // For Metropolis-Hastings application we want to cancel the Diracs out 
        // between the numerator and denominator
        // Therefore we want to keep all the copied r.v..
        // If even one of them is removed in Node's operation, we set the pdf to zero
        for (size_t i = 0; i < input_seq.Size(); i++) {
            for (size_t j = 0; j < output_seq.Size(); j++) {
                if (output_seq.strategies[i].get() == input_seq.strategies[i].get()) {
                    shared_rv.push_back(j);
                    break;
                }
            }
            if (shared_rv.size() <= i) {
                // If not found we set the index to -1
                shared_rv.push_back(-1);
            }
        }

        SourceDistributionHelper source_distribution_helper;
        cast_to_constant(input_seq, input_buffer);
        return Pdf_(mutation, output_seq, input_buffer, shared_rv, source_distribution_helper);
    }
}

MTS_NAMESPACE_BEGIN

struct OffsetSampler {
    Real Uniform1D() {
        return Real(0);
    }

    std::array<Real, 2> Uniform2D() {
        return {{previousScreenPosition[0] + (offset[0] * Real(invResolution.x)),
                previousScreenPosition[1] + (offset[1] * Real(invResolution.y))}};
    }

    const Point2 &previousScreenPosition;
    const Point2i &offset;
    const Vector2 &invResolution;
};

struct gradient_perturbation_t {
    template <typename T>
    auto perturb(Context<T>& context, const RandomSequence<Vertex>& path) const {
        const Sensor *sensor = scene->getSensor();
        const Vector2 invResolution = sensor->getInvResolution();

        auto uniformScreenSpace = make_random_vector(make_random_var(u1), make_random_var(u2));

        Point2 previousScreenPosition;
        bool insideScreen = project(sensor, 
            to_point(path[0].Value()), to_point(path[1].Value()),
            path[1].Get(intersection_).time, previousScreenPosition);
        previousScreenPosition[0] *= Real(invResolution.x);
        previousScreenPosition[1] *= Real(invResolution.y);

        OffsetSampler offsetSampler{previousScreenPosition, offset, invResolution};
        auto uv = context.Uniform2D(offsetSampler);
        auto newScreenPosition = uniformScreenSpace.Sample(uv);

        insideScreen = insideScreen && (
            at<0>(newScreenPosition).Value() >= 0 && at<0>(newScreenPosition).Value() < 1 &&
            at<1>(newScreenPosition).Value() >= 0 && at<1>(newScreenPosition).Value() < 1);

        auto sensorSampler = sensor->makeSampler();
        auto dir = sensorSampler.Sample(sensor_sampling_tag::Direction{},
            path[0], context, uniDist, newScreenPosition);
        auto p = constant(path[0].Value());

        // Use Mitsuba's raycasting engine to obtain intersection information
        const Ray ray(to_point(p.Value()), to_vector(dir.Value()), path[1].Get(intersection_).time);
        const Intersection its = context.constant_call(raycaster, ray);
        // Do symbolic intersection
        auto intersectSample = intersect(p, dir, its);
        // TODO: setup its.emitter
        const Emitter *emitter = its.shape != nullptr ? its.shape->getEmitter() : nullptr;

        return optional_sample(
            insideScreen
            , intersectSample
            , emitter_ = emitter
        );
    }

    template <typename T>
    auto operator()(Context<T>& context, const RandomSequence<Vertex>& path) const {
        slice_(path, 0, 1, eyeSubpath);
        eyeSubpath.Append(perturb(context, path));
        int index = 1;
        for(;index < path.Size();) {
            for (; isSpecular(path[index], specularRoughnessThreshold); index++) {
                if (index >= path.Size() - 1) {
                    break;
                }
                // Check whether the original interaction was reflection or refraction
                // Keep reflect/refract
                const auto prevWi = to_vector((path[index - 1].Value() - path[index].Value()).normalized());
                const auto prevWo = to_vector((path[index + 1].Value() - path[index].Value()).normalized());
                const Intersection &prevIts = path[index].Get(intersection_);
                const auto prevLocalWi = prevIts.toLocal(prevWi);
                const auto prevLocalWo = prevIts.toLocal(prevWo);
                const bool isReflect = prevLocalWi.z * prevLocalWo.z > Float(0);
                AppendBSDF(eyeSubpath, uniDist, raycaster, -1,
                           isReflect ? BSDF::EReflection : BSDF::ETransmission);
            }

            if (index == path.Size() - 1) {
                break;
            }

            if (!isSpecular(path[index + 1], specularRoughnessThreshold)) {
                break;
            }

            // Perturb a very small direction to mimic dirac
            auto perturbDir = perturb_direction{degToRad(0.0001f), degToRad(0.001f), uniDist, raycaster};
            const aether::Vector3 previousDir =
                (path[index + 1].Value() - path[index].Value()).normalized();
            eyeSubpath.Append(perturbDir(context,
                path[index].Value(), previousDir, path[index].Get(intersection_).time));
            index++;
        }
        SAssert(path.Size() >= eyeSubpath.Size());
        slice_(path, eyeSubpath.Size(), path.Size() - eyeSubpath.Size(), pathSuffix);
        concat_(eyeSubpath, pathSuffix, proposalPath);
        return proposalPath;
    }

    const Scene *scene;
    UniDist &uniDist;
    const Raycast &raycaster;
    const Real specularRoughnessThreshold;
    const Point2i offset;

    // buffers
    mutable RandomSequence<Vertex> eyeSubpath;
    mutable RandomSequence<Vertex> pathSuffix;
    mutable RandomSequence<Vertex> proposalPath;
};


class MyGradientDomainPathTracer : public Integrator {
public:
    MyGradientDomainPathTracer(const Properties &props) : Integrator(props) {
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

        m_specularRoughnessThreshold = props.getFloat("specularRoughnessThreshold", 0.02);

        m_doL1 = props.getBoolean("doL1Reconstruction", true);
    }

    /// Unserialize from a binary data stream
    MyGradientDomainPathTracer(Stream *stream, InstanceManager *manager)
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
        // Setup image buffers
        Vector2i cropSize = film->getCropSize();
        ref<ImageBlock> primalImageBlock = new ImageBlock(Bitmap::ESpectrum,
            cropSize, film->getReconstructionFilter());
        ref<ImageBlock> dxImageBlock = new ImageBlock(Bitmap::ESpectrum,
            cropSize, film->getReconstructionFilter());
        ref<ImageBlock> dyImageBlock = new ImageBlock(Bitmap::ESpectrum,
            cropSize, film->getReconstructionFilter());
        ref<ImageBlock> visibleEmitterImage = new ImageBlock(Bitmap::ESpectrum,
            cropSize, film->getReconstructionFilter());
        primalImageBlock->setOffset(
                Point2i(primalImageBlock->getBorderSize(), primalImageBlock->getBorderSize()));
        primalImageBlock->clear();
        dxImageBlock->setOffset(Point2i(dxImageBlock->getBorderSize(), dxImageBlock->getBorderSize()));
        dxImageBlock->clear();
        dxImageBlock->setWarn(false);
        dyImageBlock->setOffset(Point2i(dyImageBlock->getBorderSize(), dyImageBlock->getBorderSize()));
        dyImageBlock->clear();
        dyImageBlock->setWarn(false);
        visibleEmitterImage->setOffset(
                Point2i(visibleEmitterImage->getBorderSize(), visibleEmitterImage->getBorderSize()));
        visibleEmitterImage->clear();
        ref<Bitmap> bitmap = new Bitmap(Bitmap::ESpectrum, Bitmap::EFloat, film->getSize());
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
                        render(scene,
                               clonedSampler,
                               visibleEmitterImage.get(),
                               primalImageBlock.get(),
                               dxImageBlock.get(),
                               dyImageBlock.get(),
                               x,
                               y,
                               mutex);
                    }
                }
            }
            std::lock_guard<std::mutex> lock(mutex);
            bitmap->copyFrom(primalImageBlock->getBitmap());
            film->setBitmap(bitmap);
            queue->signalRefresh(job);
            reporter.Update(1);
        }, Vector2i(nXTiles, nYTiles));
        TerminateWorkerThreads();

        if (unfinished) {
            return false;
        }

        reporter.Done();

        //dxImageBlock->getBitmap()->write(fs::path("dx.exr"));
        //dyImageBlock->getBitmap()->write(fs::path("dy.exr"));

        Log(EInfo, "Rendering finished.  Executing poisson solver.");
        ps::Image3 throughputImage(cropSize.x, cropSize.y);
        ps::Image3 dxImage(cropSize.x, cropSize.y);
        ps::Image3 dyImage(cropSize.x, cropSize.y);
        ps::Image3 visibleEmitterImage3(cropSize.x, cropSize.y);
        // Transform the data for the solver.
        for (int y = 0; y < cropSize.y; y++) {
            for (int x = 0; x < cropSize.x; x++) {
                for (int c = 0; c < 3; c++) {
                    const int width = primalImageBlock->getBitmap()->getWidth();
                    //const int border = primalImageBlock->getBorderSize();
                    const int index = 3 * (y * width + x) + c;
                    throughputImage.At(x, y)[c] = primalImageBlock->getBitmap()->getFloatData()[index];
                    dxImage.At(x, y)[c] = dxImageBlock->getBitmap()->getFloatData()[index];
                    dyImage.At(x, y)[c] = dyImageBlock->getBitmap()->getFloatData()[index];
                    visibleEmitterImage3.At(x, y)[c] =
                        visibleEmitterImage->getBitmap()->getFloatData()[index];
                }
            }
        }
        ps::Image3 reconstructedImage(cropSize.x, cropSize.y);
        ps::Solve(throughputImage, dxImage, dyImage, visibleEmitterImage3, 0.04, m_doL1, reconstructedImage);
        // Transform the data back
        for (int y = 0; y < cropSize.y; y++) {
            for (int x = 0; x < cropSize.x; x++) {
                for (int c = 0; c < 3; c++) {
                    const int width = bitmap->getWidth();
                    const int index = 3 * (y * width + x) + c;
                    bitmap->getFloatData()[index] = reconstructedImage.At(x, y)[c];
                }
            }
        }
        film->setBitmap(bitmap);
        queue->signalRefresh(job);

        return true;
    }

    bool isValid(const double pdf) {
        return pdf > 0 && std::isfinite(pdf);
    }

    void render(const Scene *scene, Sampler *sampler,
            ImageBlock *visibleEmitterImage, ImageBlock *primalImageBlock,
            ImageBlock *gradientXImageBlock, ImageBlock *gradientYImageBlock,
            const int x, const int y, std::mutex &mutex) {
        UniDist uniDist(sampler);

        // Sample time
        // TODO: make this optional
        const Float time = sampler->next1D();

        const Point2i offsets[] = {
            Point2i(-1,  0),
            Point2i( 1,  0),
            Point2i( 0, -1),
            Point2i( 0,  1)
        };
        constexpr int numOffsets = sizeof(offsets) / sizeof(Point2i);
        // Keep in mind that offset paths are never actually traced, they merely store the strategies for
        // MIS computation
        RandomSequence<Vertex> offsetSubpaths[numOffsets];
        // Set up the gradient shifts
        Node<gradient_perturbation_t> perturbGradient[numOffsets] = {
            Node<gradient_perturbation_t>(scene, uniDist, *m_raycaster, m_specularRoughnessThreshold, offsets[0]),
            Node<gradient_perturbation_t>(scene, uniDist, *m_raycaster, m_specularRoughnessThreshold, offsets[1]),
            Node<gradient_perturbation_t>(scene, uniDist, *m_raycaster, m_specularRoughnessThreshold, offsets[2]),
            Node<gradient_perturbation_t>(scene, uniDist, *m_raycaster, m_specularRoughnessThreshold, offsets[3])
        };

        // Sample camera subpath
        RandomSequence<Vertex> sensorSubpath;
        AppendPositionOnSensor(sensorSubpath, uniDist, scene->getSensor(), time);
        AppendDirectionFromSensor(sensorSubpath, uniDist, *m_raycaster, x, y);
        sensorSubpath.Sample();
        for (int i = 0; i < numOffsets; i++) {
            AppendPositionOnSensor(offsetSubpaths[i], uniDist, scene->getSensor(), time);
            AppendDirectionFromSensor(offsetSubpaths[i], uniDist, *m_raycaster, x, y);
        }

        // Loop until reaching specified maximum depth
        // e.g. if maxDepth == 2, we don't need sensorSubpath.Size() > 3
        for(;sensorSubpath.Size() <= m_maxDepth;) {
            AppendBSDF(sensorSubpath, uniDist, *m_raycaster);
            for (int i = 0; i < numOffsets; i++) {
                AppendBSDF(offsetSubpaths[i], uniDist, *m_raycaster);
            }
            sensorSubpath.Sample();
        }

        // Combine
        const ref_vector<Emitter> &emitters = scene->getEmitters();
        OcclusionCache occlusionCache(scene, time);
        Spectrum primary(0.f);
        Spectrum gradients[numOffsets];
        for (int i = 0; i < numOffsets; i++) {
            gradients[i] = Spectrum(0.f);
        }
        RandomSequence<Vertex> buffer;
        for (int pathLength = 2; pathLength <= m_maxDepth; pathLength++) {
            // Sample emitter subpath
            auto emitterSamplingPath = sensorSubpath.Slice(0, pathLength);
            if (!emitterSamplingPath.Back().Get(intersection_).isValid()) {
                continue;
            }
            AppendDirectSampleEmitter(emitterSamplingPath, uniDist, emitters);
            emitterSamplingPath.Sample();
            auto bsdfSamplingPath = sensorSubpath.Slice(0, pathLength + 1);

            TSpectrum<double, SPECTRUM_SAMPLES> emitterEval(0.f);
            const double pdfEmitter = emitterSamplingPath.Pdf(emitterSamplingPath);
            if (isValid(pdfEmitter) &&
                    !occlusionCache.Query(
                        std::make_pair(emitterSamplingPath[emitterSamplingPath.Size() - 2].Value(),
                                       emitterSamplingPath[emitterSamplingPath.Size() - 1].Value())
                        )) {
                emitterEval = evalBidir(scene, emitterSamplingPath);
                if (!emitterEval.isZero()) {
                    const double pdfBSDF = bsdfSamplingPath.Pdf(emitterSamplingPath);
                    const double weight = misWeight(pdfEmitter, pdfBSDF);
                    Spectrum emitterPrimary = weight * Spectrum(emitterEval / pdfEmitter);
                    primary += emitterPrimary;
                }
            }
            const double pdfBSDF = bsdfSamplingPath.Pdf(bsdfSamplingPath);
            TSpectrum<double, SPECTRUM_SAMPLES> bsdfEval(0.f);
            if (isValid(pdfBSDF)) {
                bsdfEval = evalBidir(scene, bsdfSamplingPath);
                if (!bsdfEval.isZero()) {
                    const double pdfEmitter = emitterSamplingPath.Pdf(bsdfSamplingPath);
                    const double weight = misWeight(pdfBSDF, pdfEmitter);
                    Spectrum bsdfPrimary = weight * Spectrum(bsdfEval / pdfBSDF);
                    primary += bsdfPrimary;
                }
            }

            if (!isValid(pdfEmitter) && !isValid(pdfBSDF)) {
                continue;
            }

            for (int i = 0; i < numOffsets; i++) {
                auto emitterSamplingPathShift = emitterSamplingPath.Mutate(perturbGradient[i]);
                auto bsdfSamplingPathShift = bsdfSamplingPath.Mutate(perturbGradient[i]);
                // For MIS computation
                auto offsetEmitterSamplingPath = offsetSubpaths[i].Slice(0, pathLength);
                AppendDirectSampleEmitter(offsetEmitterSamplingPath, uniDist, emitters);
                auto offsetBSDFSamplingPath = offsetSubpaths[i].Slice(0, pathLength + 1);

                // Compute emitter sampling contribution
                TSpectrum<double, SPECTRUM_SAMPLES> shiftEmitterEval(0.);
                Float emitterJacobian =
                    jacobian(perturbGradient[i], emitterSamplingPath, emitterSamplingPathShift, buffer);
                if (emitterSamplingPathShift.AllValid() &&
                        !occluded(emitterSamplingPathShift, occlusionCache)) {
                    shiftEmitterEval = emitterJacobian * evalBidir(scene, emitterSamplingPathShift);
                }

                {
                    // 4 possible strategies here: two from this path, two from offset path
                    if (pdfEmitter > 0 && std::isfinite(pdfEmitter)) {
                        const double pdfBSDF = bsdfSamplingPath.Pdf(emitterSamplingPath);
                        const double pdfOffsetEmitter =
                            offsetEmitterSamplingPath.Pdf(emitterSamplingPath) * emitterJacobian;
                        const double pdfOffsetBSDF =
                            offsetBSDFSamplingPath.Pdf(emitterSamplingPath) * emitterJacobian;
                        const double weight =
                            misWeight(0, {pdfEmitter, pdfBSDF, pdfOffsetEmitter, pdfOffsetBSDF});
                        Spectrum emitterGradient =
                            weight * Spectrum((shiftEmitterEval - emitterEval) / pdfEmitter);
                        gradients[i] += emitterGradient;
                    }
                }

                // Compute bsdf sampling contribution
                TSpectrum<double, SPECTRUM_SAMPLES> shiftBSDFEval(0.);
                Float bsdfJacobian =
                    jacobian(perturbGradient[i], bsdfSamplingPath, bsdfSamplingPathShift, buffer);
                if (bsdfSamplingPathShift.AllValid() && !occluded(bsdfSamplingPathShift, occlusionCache)) {
                    shiftBSDFEval = bsdfJacobian * evalBidir(scene, bsdfSamplingPathShift);
                } else {
                    bsdfJacobian = 0.;
                }

                {
                    // 4 possible strategies here: two from this path, two from offset path
                    if (pdfBSDF > 0 && std::isfinite(pdfBSDF)) {
                        const double pdfEmitter = emitterSamplingPath.Pdf(bsdfSamplingPath);
                        const double pdfOffsetEmitter =
                            offsetEmitterSamplingPath.Pdf(bsdfSamplingPath) * bsdfJacobian;
                        const double pdfOffsetBSDF =
                            offsetBSDFSamplingPath.Pdf(bsdfSamplingPath) * bsdfJacobian;
                        const double weight =
                            misWeight(1, {pdfEmitter, pdfBSDF, pdfOffsetEmitter, pdfOffsetBSDF});
                        Spectrum bsdfGradient = weight * Spectrum((shiftBSDFEval - bsdfEval) / pdfBSDF);
                        gradients[i] += bsdfGradient;
                    }
                }
            }
        }

        std::lock_guard<std::mutex> lock(mutex);
        // Splat contributions
        Point2 screenPosition;
        bool inside = project(scene, sensorSubpath, screenPosition);
        if (!inside) {
            return;
        }
        const Float invSampleCount = Float(1) / sampler->getSampleCount();
        primary *= invSampleCount;
        for (int i = 0; i < numOffsets; i++) {
            gradients[i] *= invSampleCount;
        }
        gradients[0] *= Float(-1);
        gradients[2] *= Float(-1);
        Assert(!primary.isNaN());
        Assert(!gradients[0].isNaN());
        Assert(!gradients[1].isNaN());
        Assert(!gradients[2].isNaN());
        Assert(!gradients[3].isNaN());
        primalImageBlock->put(screenPosition, &primary[0]);
        gradientXImageBlock->put(screenPosition - Vector2(1, 0), &(gradients[0][0]));
        gradientXImageBlock->put(screenPosition, &(gradients[1][0]));
        gradientYImageBlock->put(screenPosition - Vector2(0, 1), &(gradients[2][0]));
        gradientYImageBlock->put(screenPosition, &(gradients[3][0]));

        if (!m_hideEmitters) {
            auto path = sensorSubpath.Slice(0, 2);
            if (path.AllValid()) {
                Spectrum contribution = estimateBidir(scene, make_view(path));
                visibleEmitterImage->put(screenPosition, &contribution[0]);
            }
        }
    }

    template<typename Mutation>
    double jacobian(Mutation &mutation,
            const RandomSequence<Vertex> &primaryPath,
            const RandomSequence<Vertex> &shiftPath,
            RandomSequence<Vertex> &buffer) const {
        const double pPrimaryGivenShift = ConditionalPdf(mutation, primaryPath, shiftPath, buffer);
        const double pShiftGivenPrimary = ConditionalPdf(mutation, shiftPath, primaryPath, buffer);
        if (!std::isfinite(pShiftGivenPrimary) || !std::isfinite(pPrimaryGivenShift) ||
                pPrimaryGivenShift <= 0. || pShiftGivenPrimary <= 0.) {
            return 0.;
        }
        return pPrimaryGivenShift / pShiftGivenPrimary;
    }

    bool occluded(const RandomSequence<Vertex> &path, OcclusionCache &occlusionCache) const {
        for (size_t i = 0; i < path.Size() - 1; i++) {
            if (occlusionCache.Query(std::make_pair(path[i].Value(), path[i + 1].Value()))) {
                return true;
            }
        }
        return false;
    }

    std::string toString() const {
        std::ostringstream oss;
        oss << "MyGradientDomainPathTracer[" << endl
            << "  hideEmitters = " << m_hideEmitters << endl
            << "]";
        return oss.str();
    }

    MTS_DECLARE_CLASS()
private:
    bool m_hideEmitters;
    int m_rrDepth;
    int m_maxDepth;
    Real m_specularRoughnessThreshold;
    bool m_doL1;

    std::unique_ptr<Raycast> m_raycaster;
    bool m_running;
};

MTS_IMPLEMENT_CLASS_S(MyGradientDomainPathTracer, false, SamplingIntegrator)
MTS_EXPORT_PLUGIN(MyGradientDomainPathTracer, "My gradient domain path tracer");
MTS_NAMESPACE_END
