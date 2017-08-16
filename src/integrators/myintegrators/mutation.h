#pragma once

#include "markovchainstate.h"
#include "classification.h"
#include <mitsuba/core/spectrum.h>
#include <mitsuba/render/scene.h>
#include <aether/RandomSequence.h>

MTS_NAMESPACE_BEGIN

struct Mutation;

double acceptProbability(const std::vector<std::shared_ptr<Mutation>> &mutations,
						const MarkovChainState &currentState,
						const MarkovChainState &proposalState,
						const Mutation *selectedMutation);

struct Mutation {
	// Mutate and fill proposalState
	virtual void Mutate(const Scene *scene, const MarkovChainState &currentState,
		MarkovChainState &proposalState, OcclusionCache &occlusionCache) = 0;

	// Return conditional pdf
	virtual double Pdf(const RandomSequence<Vertex> &target, const RandomSequence<Vertex> &given) = 0;

	virtual bool Mutable(const RandomSequence<Vertex> &path) const = 0;

	// Use for gathering statistics
	virtual void Accepted() const {};

	template<typename T>
	void Mutate(Node<T> &mutationStrategy, const Scene *scene, const MarkovChainState &currentState,
					MarkovChainState &proposalState, OcclusionCache &occlusionCache) {
		RandomSequence<Vertex> proposalPath = currentState.path.Mutate(mutationStrategy);
		proposalState = f(scene, proposalPath, occlusionCache);
	}
};

MTS_NAMESPACE_END
