#include "mutation.h"

MTS_NAMESPACE_BEGIN

struct CausticsPerturbation : public Mutation {
    CausticsPerturbation(const Real theta_1, const Real theta_2,
            const Float specularRoughnessThreshold, UniDist &uniDist, const Raycast &raycaster)
        : theta_1(theta_1), theta_2(theta_2),
          specularRoughnessThreshold(specularRoughnessThreshold), uniDist(uniDist), raycaster(raycaster) {}

    void Mutate(const Scene *scene, const MarkovChainState &currentState,
            MarkovChainState &proposalState, OcclusionCache &occlusionCache) override;

    double Pdf(const RandomSequence<Vertex> &target, const RandomSequence<Vertex> &given) override;

    bool Mutable(const RandomSequence<Vertex> &path) const override;

    void Accepted() const override;

    const Real theta_1;
    const Real theta_2;
    const Float specularRoughnessThreshold;
    UniDist &uniDist;
    const Raycast &raycaster;
};

MTS_NAMESPACE_END
