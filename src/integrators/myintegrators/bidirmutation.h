#pragma once

#include "mutation.h"
#include <map>

MTS_NAMESPACE_BEGIN

struct BidirectionalMutation : public Mutation {
	BidirectionalMutation(const Scene *scene, UniDist &uniDist,
			const Raycast &raycaster, const int minDepth, const int maxDepth,
			const Real specularRoughnessThreshold)
		: scene(scene), uniDist(uniDist), raycaster(raycaster), 
		minLength(minDepth + 2), maxLength(maxDepth + 2),
		specularRoughnessThreshold(specularRoughnessThreshold) {}

	void Mutate(const Scene *scene, const MarkovChainState &currentState,
			MarkovChainState &proposalState, OcclusionCache &occlusionCache) override;

	double Pdf(const RandomSequence<Vertex> &target, const RandomSequence<Vertex> &given) override;

	bool Mutable(const RandomSequence<Vertex> &path) const override {
		return true;
	}

	void Accepted() const override;

	const Scene *scene;
	UniDist &uniDist;
	const Raycast &raycaster;
	const int minLength;
	const int maxLength;
	const Real specularRoughnessThreshold;

	// memory buffers for speed improvement
	std::vector<int> deleteLength;
	std::vector<aether::Real> deleteLengthWeight;
	std::vector<int> desiredLength;
	std::vector<aether::Real> desiredLengthWeight;
	std::vector<int> deleteBegin;
	std::vector<int> subpathInsertSize;

	// dynamic r.v. caches
	// not sure if std::map is better than std::unordered_map
	std::map<int, discrete_random_var_dynamic<int>> desiredLengthRV;
	std::map<std::pair<int, int>, discrete_random_var_dynamic<int>> deleteLengthRV;
	std::map<int, discrete_random_var_dynamic<int>> eyeSubpathInsertSizeRV;

	// RandomSequence buffers
	RandomSequence<Vertex> eyeSubpath;
	RandomSequence<Vertex> lightSubpath, reversedLightSubpath;
	RandomSequence<Vertex> proposalPath;
};

MTS_NAMESPACE_END
