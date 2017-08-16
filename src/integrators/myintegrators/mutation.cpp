#include "mutation.h"

MTS_NAMESPACE_BEGIN

double acceptProbability(const std::vector<std::shared_ptr<Mutation>> &mutations,
                         const MarkovChainState &currentState,
                         const MarkovChainState &proposalState,
                         const Mutation *selectedMutation) {
    if (!proposalState.Valid() || proposalState.importance <= 0 || !proposalState.contribution.isValid()) {
        return 0;
    }

    if (!currentState.Valid() || currentState.importance <= 0 || !currentState.contribution.isValid()) {
        return 1;
    }

    const auto &proposalPath = proposalState.path;
    const auto &currentPath = currentState.path;

    double pCurrentGivenProposal = 0.;
    int proposalMutableCount = 0;
    for (const auto &mutation : mutations) {
        if (mutation->Mutable(proposalPath)) {
            proposalMutableCount++;
            if ((selectedMutation == nullptr || selectedMutation == mutation.get())) {
                pCurrentGivenProposal += mutation->Pdf(currentPath, proposalPath);
            }
        }
    }
    if (proposalMutableCount == 0) {
        return 0.;
    }
    pCurrentGivenProposal /= (double)proposalMutableCount;

    double pProposalGivenCurrent = 0.;
    int currentMutableCount = 0;
    for (const auto &mutation : mutations) {
        if (mutation->Mutable(currentPath)) {
            currentMutableCount++;
            if ((selectedMutation == nullptr || selectedMutation == mutation.get())) {
                pProposalGivenCurrent += mutation->Pdf(proposalPath, currentPath);
            }
        }
    }
    if (currentMutableCount == 0) {
        return 0.;
    }
    pProposalGivenCurrent /= (double)currentMutableCount;

    if (!std::isfinite(pCurrentGivenProposal) || !std::isfinite(pProposalGivenCurrent) ||
        pProposalGivenCurrent <= 0. || pCurrentGivenProposal <= 0.) {
        return 0;
    }

    return std::min(
            (proposalState.importance / pProposalGivenCurrent) / 
            (currentState.importance / pCurrentGivenProposal),
                1.);
}

MTS_NAMESPACE_END
