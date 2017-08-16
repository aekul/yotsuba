#include "causticsperturbation.h"
#include "intersect.h"
#include "appendbsdf.h"
#include "perturbdirection.h"

#include <mitsuba/render/aetheremitters_impl.h>
#include <mitsuba/core/statistics.h>

MTS_NAMESPACE_BEGIN

static StatsCounter statsAccepted("My caustics perturbation",
		"Accepted mutations", EPercentage);

struct caustics_perturbation_t {
	CausticsPerturbation *mutation;

	template <typename T>
	auto operator()(Context<T>& context, const RandomSequence<Vertex>& path) const {
		UniDist &uniDist = mutation->uniDist;
		const Raycast &raycaster = mutation->raycaster;
		if (path.Size() < 3) {
			// If there is no surface vertex, reject the mutation
			return RandomSequence<Vertex>{};
		}

		// Find the first diffuse vertex counting from camera
		int splitPosition = 0;
		for (int i = 2; i < path.Size(); i++) {
			if (i == path.Size() - 1 || 
					isDiffuse(path[i], mutation->specularRoughnessThreshold)) {
				splitPosition = i;
				break;
			}
		}
		SAssert(splitPosition >= 2);

		// Compute the pertubation angle spread (see Veach's thesis)
		const Real lengthE = (path[0].Value() - path[1].Value()).norm();
		Real lengthL(0);
		for (int i = 1; i < splitPosition; i++) {
			lengthL += (path[i].Value() - path[i + 1].Value()).norm();
		}
		SAssert(lengthL > 0);
		const Real factor = lengthE / lengthL;
		const Real theta_1 = mutation->theta_1 * factor;
		const Real theta_2 = mutation->theta_2 * factor;

		// Now slice, perturb, and append to construct the new path
		auto lightSubpath = reverse_(path.Slice(splitPosition, path.Size() - splitPosition));
		auto perturbDir = perturb_direction{theta_1, theta_2, uniDist, raycaster};
		const aether::Vector3 previousDir =
			(path[splitPosition - 1].Value() - path[splitPosition].Value()).normalized();
		lightSubpath.Append(perturbDir(context,
			path[splitPosition].Value(), previousDir, path[splitPosition].Get(intersection_).time));
		int index = splitPosition - 1;
		// Propagate through specular interactions
		while (lightSubpath.Size() < path.Size() - 1) {
			// If we hit a non-specular surface, reject the path since this is irreversible
			if (!isSpecular(path[index], mutation->specularRoughnessThreshold)) {
				return RandomSequence<Vertex>{};
			}
			// Find out whether the next interaction is reflection or refraction
			// Use it to choose the corresponding BSDF layer
			const auto prevWi = to_vector((path[index - 1].Value() - path[index].Value()).normalized());
			const auto prevWo = to_vector((path[index + 1].Value() - path[index].Value()).normalized());
			const Intersection &prevIts = path[index].Get(intersection_);
			const auto prevLocalWi = prevIts.toLocal(prevWi);
			const auto prevLocalWo = prevIts.toLocal(prevWo);
			const bool isReflect = prevLocalWi.z * prevLocalWo.z > Float(0);
			AppendBSDF(lightSubpath, uniDist, raycaster, -1,
					   isReflect ? BSDF::EReflection : BSDF::ETransmission);
			index--;
		}
		auto ret = path.Slice(0, 1).Concat(reverse_(lightSubpath));
		SAssert(ret.Size() == path.Size());
		return ret;
	}
};

void CausticsPerturbation::Mutate(const Scene *scene, const MarkovChainState &currentState,
		MarkovChainState &proposalState, OcclusionCache &occlusionCache) {
	SAssert(!isSpecular(currentState.path[1], specularRoughnessThreshold));
	statsAccepted.incrementBase(1);
	Node<caustics_perturbation_t> strategy{this};
	Mutation::Mutate(strategy, scene, currentState, proposalState, occlusionCache);
}

double CausticsPerturbation::Pdf(const RandomSequence<Vertex> &target, const RandomSequence<Vertex> &given) {
	Node<caustics_perturbation_t> strategy{this};
	return ConditionalPdf(strategy, target, given);
}

bool CausticsPerturbation::Mutable(const RandomSequence<Vertex> &path) const {
	if (!path.AllValid()) {
		return false;
	}
	if (path.Size() >= 3 && isDiffuse(path[1], specularRoughnessThreshold)) {
		// Only do LS*D+E path
		return true;
	}
	return false;
}

void CausticsPerturbation::Accepted() const {
	++statsAccepted;
}

MTS_NAMESPACE_END
