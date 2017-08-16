#include "lensperturbation.h"
#include "intersect.h"
#include "appendbsdf.h"
#include "perturblens.h"

#include <mitsuba/render/aethersensors_impl.h>
#include <mitsuba/core/statistics.h>

MTS_NAMESPACE_BEGIN

static StatsCounter statsAccepted("My lens perturbation",
        "Accepted mutations", EPercentage);

struct lens_perturbation_t {
    LensPerturbation *mutation;

    template <typename T>
    auto operator()(Context<T>& context, const RandomSequence<Vertex>& path) const {
        const Scene *scene = mutation->scene;
        const Real minJump = mutation->minJump;
        const Real coveredArea = mutation->coveredArea;
        UniDist &uniDist = mutation->uniDist;
        const Raycast &raycaster = mutation->raycaster;

        auto eyeSubpath = path.Slice(0, 1);
        // Perturb the screen position of the path
        auto perturb = perturb_lens{scene, minJump, coveredArea, uniDist, raycaster};
        eyeSubpath.Append(perturb(context, path));
        // Propagate through specular interactions
        for (int i = 1; isSpecular(path[i], mutation->specularRoughnessThreshold); i++) {
            // Find out whether the next interaction is reflection or refraction
            // Use it to choose the corresponding BSDF layer
            const auto prevWi = to_vector((path[i - 1].Value() - path[i].Value()).normalized());
            const auto prevWo = to_vector((path[i + 1].Value() - path[i].Value()).normalized());
            const Intersection &prevIts = path[i].Get(intersection_);
            const auto prevLocalWi = prevIts.toLocal(prevWi);
            const auto prevLocalWo = prevIts.toLocal(prevWo);
            const bool isReflect = prevLocalWi.z * prevLocalWo.z > Float(0);
            AppendBSDF(eyeSubpath, uniDist, raycaster, -1,
                       isReflect ? BSDF::EReflection : BSDF::ETransmission);
        }
        SAssert(path.Size() >= eyeSubpath.Size());
        auto suffix = path.Slice(eyeSubpath.Size(), path.Size() - eyeSubpath.Size());
        auto ret = eyeSubpath.Concat(suffix);
        SAssert(ret.Size() == path.Size());
        return ret;
    }
};

void LensPerturbation::Mutate(const Scene *scene, const MarkovChainState &currentState,
        MarkovChainState &proposalState, OcclusionCache &occlusionCache) {
    statsAccepted.incrementBase(1);
    Node<lens_perturbation_t> strategy{this};
    Mutation::Mutate(strategy, scene, currentState, proposalState, occlusionCache);
}

double LensPerturbation::Pdf(const RandomSequence<Vertex> &target, const RandomSequence<Vertex> &given) {
    Node<lens_perturbation_t> strategy{this};
    return ConditionalPdf(strategy, target, given);
}

bool LensPerturbation::Mutable(const RandomSequence<Vertex> &path) const {
    if (path.Size() < 2 || !path.AllValid()) {
        return false;
    }

    int firstDiffuseIndex = 1;
    while (isSpecular(path[firstDiffuseIndex], specularRoughnessThreshold)) {
        firstDiffuseIndex++;
        if (firstDiffuseIndex >= path.Size()) {
            return false;
        }
    }

    if (firstDiffuseIndex >= path.Size() - 1) {
        return true;
    }

    return isDiffuse(path[firstDiffuseIndex]    , specularRoughnessThreshold) &&
         !isSpecular(path[firstDiffuseIndex + 1], specularRoughnessThreshold);
}

void LensPerturbation::Accepted() const {
    ++statsAccepted;
}

MTS_NAMESPACE_END
