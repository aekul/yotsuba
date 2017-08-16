#include "multichainperturbation.h"
#include "intersect.h"
#include "appendbsdf.h"
#include "perturblens.h"
#include "perturbdirection.h"

#include <mitsuba/render/aethersensors_impl.h>
#include <mitsuba/core/statistics.h>

MTS_NAMESPACE_BEGIN

static StatsCounter statsAccepted("My multi-chain perturbation",
		"Accepted mutations", EPercentage);

struct multichain_perturbation_t {
	MultiChainPerturbation *mutation;

	template <typename T>
	auto operator()(Context<T>& context, const RandomSequence<Vertex>& path) const {
	    const Scene *scene = mutation->scene;
		const Real minJump = mutation->minJump;
		const Real coveredArea = mutation->coveredArea;
		const Real theta_1 = mutation->theta_1;
		const Real theta_2 = mutation->theta_2;
		UniDist &uniDist = mutation->uniDist;
		const Raycast &raycaster = mutation->raycaster;

		auto eyeSubpath = path.Slice(0, 1);
		// Perturb the screen position of the path
		auto perturb = perturb_lens{scene, minJump, coveredArea, uniDist, raycaster};
		eyeSubpath.Append(perturb(context, path));
		int index = 1;
		// Propagate through multiple specular chains
		for(;index < path.Size();) {
			for (; isSpecular(path[index], mutation->specularRoughnessThreshold); index++) {
				// Find out whether the next interaction is reflection or refraction
				// Use it to choose the corresponding BSDF layer
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

			if (!isSpecular(path[index + 1], mutation->specularRoughnessThreshold)) {
				break;
			}

			auto perturbDir = perturb_direction{theta_1, theta_2, uniDist, raycaster};
	    	const aether::Vector3 previousDir =
	    		(path[index + 1].Value() - path[index].Value()).normalized();
			eyeSubpath.Append(perturbDir(context,
				path[index].Value(), previousDir, path[index].Get(intersection_).time));
			index++;
		}
		
		SAssert(path.Size() >= eyeSubpath.Size());
		auto suffix = path.Slice(eyeSubpath.Size(), path.Size() - eyeSubpath.Size());
		return eyeSubpath.Concat(suffix);
	}
};

void MultiChainPerturbation::Mutate(const Scene *scene, const MarkovChainState &currentState,
		MarkovChainState &proposalState, OcclusionCache &occlusionCache) {
	statsAccepted.incrementBase(1);
	Node<multichain_perturbation_t> strategy{this};
	Mutation::Mutate(strategy, scene, currentState, proposalState, occlusionCache);
}

double MultiChainPerturbation::Pdf(const RandomSequence<Vertex> &target, const RandomSequence<Vertex> &given) {
	Node<multichain_perturbation_t> strategy{this};
	return ConditionalPdf(strategy, target, given);
}

bool MultiChainPerturbation::Mutable(const RandomSequence<Vertex> &path) const {
	if (!path.AllValid()) {
		return false;
	}

	for (int i = 1; i < path.Size() - 1; i++) {
		if (!isSpecular(path[i], specularRoughnessThreshold) &&
				!isSpecular(path[i + 1], specularRoughnessThreshold)) {
			return false;
		}
		if (!isSpecular(path[i], specularRoughnessThreshold)) {
			return true;
		}
	}

	// a pure specular path, let lens perturbation handles it
	return false;
}

void MultiChainPerturbation::Accepted() const {
	++statsAccepted;
}

MTS_NAMESPACE_END
