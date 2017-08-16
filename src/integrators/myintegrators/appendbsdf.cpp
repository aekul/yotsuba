#include "appendbsdf.h"
#include "intersect.h"
#include <mitsuba/render/scene.h>
#include <mitsuba/render/aetherbsdfs_impl.h>
#include <mitsuba/render/randomframe.h>

MTS_NAMESPACE_BEGIN

// Importance sample bsdf, possiblely terminate with Russian roulette
struct sample_bsdf_rr_t {
    template <typename T>
    auto operator()(Context<T>& context,
                    const RandomSequence<Vertex>& path,
                    UniDist& uniDist,
                    const Raycast& raycaster) const {
        // Need at least the current vertex and the incoming direction
        assert(path.Size() >= 2);
        const Intersection &currentIts = path.Back().Get(intersection_);
        auto shadingFrame = make_random_frame(currentIts);
        // Get incoming direction (pointing toward previous vertex)
        auto dirToPrevWorldVal = 
            (path[path.Size() - 2].Value() - path[path.Size() - 1].Value()).normalized();
        // Cast the incoming direction to a symbolic constant
        auto dirToPrevWorld = constant(dirToPrevWorldVal);
        // Transform the incoming direction to local frame
        auto dirToPrevLocal = to_local(shadingFrame, dirToPrevWorld);

        const BSDF *bsdf = currentIts.isValid() ? currentIts.getBSDF() : nullptr;
        // bsdfSampler randomly chooses one of the bsdf components based on incoming direction
        auto bsdfSampler = bsdf == nullptr ? bsdf_component_sampler{} :
            bsdf->makeSampler(currentIts.toLocal(to_vector(dirToPrevWorldVal)),
                              currentIts, componentMask);

        // Obtain next direction on local frame
        auto dirToNextLocal = bsdfSampler.Sample(dirToPrevLocal, context, uniDist);
        // Transform the local direction to world coordinates
        auto dirToNextWorld = to_world(shadingFrame, dirToNextLocal);

        auto currentPos = constant(path.Back().Value());
        // Do not raycast if the current vertex is not on a surface
        bool sampleNext = path.Back().Valid() && currentIts.isValid();
        if (sampleNext && rrDepth != -1 && path.Size() > rrDepth && path.Back().Valid()) {
            // Russian roulette
            const Real rrProbability = std::min((Real) bsdf->getAlbedo(currentIts).max(), (Real) 0.95);
            auto sampleNextRV = discrete_dynamic(std::vector<int>{true, false},
                                                 std::vector<Real>{rrProbability, 1.f - rrProbability});
            sampleNext = context.Sample(sampleNextRV, context.Uniform1D(uniDist));
        }

        // Use Mitsuba's raycasting engine to obtain intersection information
        const Ray ray(to_point(currentPos.Value()), to_vector(dirToNextWorld.Value()), currentIts.time);
        const Intersection nextIts = context.constant_call(raycaster, ray, sampleNext);
        // Do symbolic intersection
        auto intersectSample = intersect(currentPos, dirToNextWorld, nextIts);
        // TODO: setup nextIts.emitter
        const Emitter *emitter = nextIts.shape != nullptr ? nextIts.shape->getEmitter() : nullptr;
        return optional_sample(
            sampleNext
            , intersectSample
            , emitter_ = emitter
        );
    }

    const int rrDepth;
    const uint32_t componentMask;
};

void AppendBSDF(RandomSequence<Vertex> &path, UniDist &uniDist, const Raycast &raycaster,
        const int rrDepth, const uint32_t componentMask) {
    Node<sample_bsdf_rr_t> sampleBSDF(rrDepth, componentMask);
    path.Append(sampleBSDF, uniDist, raycaster);
}

MTS_NAMESPACE_END
