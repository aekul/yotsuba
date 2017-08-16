#include "bidirmutation.h"
#include "appendsensor.h"
#include "appendbsdf.h"
#include "appendemitter.h"

#include <mitsuba/core/statistics.h>

MTS_NAMESPACE_BEGIN

static StatsCounter statsAccepted("My bidirectional mutation",
        "Accepted mutations", EPercentage);

// Slightly modified bidirectional mutation to match Mitsuba's implementation
struct bidir_mutation_t {
    BidirectionalMutation *mutation;

    template <typename T>
    auto operator()(Context<T>& context, const RandomSequence<Vertex>& path) {
        // Set up common variables
        const Scene *scene = mutation->scene;
        UniDist &uniDist = mutation->uniDist;
        const Raycast &raycaster = mutation->raycaster;
        const int minLength = mutation->minLength;
        const int maxLength = mutation->maxLength;

        // Sample the discrete random variables that determine the portion of the path to delete
        // and the portion of the path to insert
        const int   pathLength      = path.Size() + 1;
        const auto& desiredLengthRV = getDesiredLengthRV(pathLength);
        const int   newPathLength   = context.Sample(desiredLengthRV, context.Uniform1D(uniDist));
        const int   newPathSize     = newPathLength - 1;
        if (newPathLength < minLength || newPathLength > maxLength) {
            return RandomSequence<Vertex>{};
        }
        // Speed up by early exit
        // This function always returns true when the context is SamplingContext
        // When the context is PdfContext, it takes the size of the target RandomSequence and
        // returns whether it's the same
        if (!context.MatchTargetSize(newPathSize)) {
            return RandomSequence<Vertex>{};
        }

        // Setup the discrete random variables for the path lengthes for deletion and insertion
        const auto& deleteLengthRV = getDeleteLengthRV(pathLength, newPathLength);
        const int   deleteLength   = context.Sample(deleteLengthRV, context.Uniform1D(uniDist));
        const int   insertLength   = newPathLength - pathLength + deleteLength;
        const auto& deleteBeginRV  = getDeleteBeginRV(path, deleteLength, insertLength);
        if (deleteBeginRV.Size() == 0) {
            return RandomSequence<Vertex>{};
        }
        const int   deleteBegin    = context.Sample(deleteBeginRV, context.Uniform1D(uniDist));
        const int   deleteEnd      = deleteBegin + deleteLength;
        SAssert(insertLength >= 0);

        // insertLength represents the number of edges we want to add
        // If insertLength == 1, we only need to do connection
        const int   insertSize     = insertLength - 1;
        int eyeSubpathInsertSize = 0;
        if (insertLength > 0) {
            const auto& eyeSubpathInsertSizeRV = getEyeSubpathInsertSizeRV(insertSize);
            eyeSubpathInsertSize = context.Sample(eyeSubpathInsertSizeRV, context.Uniform1D(uniDist));
        }
        const int lightSubpathInsertSize = insertSize - eyeSubpathInsertSize;

        // Set up buffers (we don't want to reallocate these paths every time we mutate)
        RandomSequence<Vertex> &eyeSubpath   = mutation->eyeSubpath;
        RandomSequence<Vertex> &lightSubpath = mutation->lightSubpath;
        RandomSequence<Vertex> &reversedLightSubpath = mutation->reversedLightSubpath;
        RandomSequence<Vertex> &proposalPath = mutation->proposalPath;

        // Delete the unwanted portion of old path
        SAssert(deleteBegin + 1 >= 0);
        slice_(path, 0, deleteBegin + 1, eyeSubpath);
        if (deleteEnd < path.Size()) {
            slice_(path, deleteEnd, path.Size() - deleteEnd, lightSubpath);
            reverse_(lightSubpath, reversedLightSubpath);
        } else {
            SAssert(deleteEnd == path.Size());
            clear(lightSubpath);
            clear(reversedLightSubpath);
        }

        // Insert the new subpath
        bool hasCameraPosition = !(deleteBegin == -1 && deleteLength >= 2);
        SAssert(path.Size() >= 2);
        const Float time = path[1].Get(intersection_).time;
        for (int i = 0; i < eyeSubpathInsertSize; ++i) {
            if (eyeSubpath.Size() == 0) {
                AppendPositionOnSensor(eyeSubpath, uniDist, scene->getSensor(), time);
                SAssert(!hasCameraPosition);
                hasCameraPosition = true;
            } else if (eyeSubpath.Size() == 1) {
                AppendDirectionFromSensor(eyeSubpath, uniDist, raycaster);
            } else {
                AppendBSDF(eyeSubpath, uniDist, raycaster);
            }
        }

        if (lightSubpathInsertSize == 1 && reversedLightSubpath.Size() == 0) {
            // Special case: do direct importance sampling
            AppendDirectSampleEmitter(eyeSubpath, uniDist, scene->getEmitters());
        } else {
            for (int i = 0; i < lightSubpathInsertSize; ++i) {
                if (reversedLightSubpath.Size() == 0) {
                    AppendPositionOnEmitter(reversedLightSubpath, uniDist, scene->getEmitters(), time);
                } else if (reversedLightSubpath.Size() == 1) {
                    AppendDirectionFromEmitter(reversedLightSubpath, uniDist, raycaster);
                } else if (reversedLightSubpath.Size() == newPathSize - 1) {
                    SAssert(false); // shouldn't happen for pinhole camera
                    AppendPositionOnSensor(reversedLightSubpath, uniDist, scene->getSensor(), time);
                    SAssert(!hasCameraPosition);
                    hasCameraPosition = true;
                } else {
                    AppendBSDF(reversedLightSubpath, uniDist, raycaster);
                }
            }
        }
        SAssert(hasCameraPosition);

        reverse_(reversedLightSubpath, lightSubpath);
        concat_(eyeSubpath, lightSubpath, proposalPath);
        SAssert(proposalPath.Size() >= 2 && proposalPath.Size() == newPathSize);
        return proposalPath;
    }


    void constructDesiredLength(const int pathLength,
            std::vector<int> &desiredLengthBuffer, std::vector<Real> &desiredLengthWeight) {
        desiredLengthBuffer.clear();
        desiredLengthWeight.clear();
        for (int length = mutation->minLength; length <= mutation->maxLength; length++) {
            const int offset = std::abs(length - pathLength);
            desiredLengthBuffer.push_back(length);
            desiredLengthWeight.push_back(pow(2.0, -offset));
        }
    }

    const discrete_random_var_dynamic<int>& getDesiredLengthRV(const int pathLength) {
        std::vector<int> &desiredLengthBuffer = mutation->desiredLength;
        std::vector<Real> &desiredLengthWeight = mutation->desiredLengthWeight;
        auto &rvCache = mutation->desiredLengthRV;
        auto it = rvCache.find(pathLength);
        if (it == rvCache.end()) {
            constructDesiredLength(pathLength, desiredLengthBuffer, desiredLengthWeight);
            return rvCache.insert({pathLength,
                discrete_dynamic(desiredLengthBuffer, desiredLengthWeight)}).first->second;
        }

        return it->second;
    }

    void constructDeleteLength(const int pathLength, const int desiredLength,
            std::vector<int> &deleteLengthBuffer, std::vector<Real> &deleteLengthWeight) {
        deleteLengthBuffer.clear();
        deleteLengthWeight.clear();
        int minDeletion = std::max((pathLength == desiredLength) ? 2 : 1, pathLength-desiredLength+1);
        for (int length = minDeletion; length <= pathLength; length++) {
            const int offset = std::abs(length - 2);
            deleteLengthBuffer.push_back(length);
            deleteLengthWeight.push_back(pow(2.0, -offset));
        }
    }

    const discrete_random_var_dynamic<int>& getDeleteLengthRV(
            const int pathLength, const int desiredLength) {
        std::vector<int>  &deleteLengthBuffer = mutation->deleteLength;
        std::vector<Real> &deleteLengthWeight = mutation->deleteLengthWeight;
        auto &rvCache = mutation->deleteLengthRV;
        auto rvKey = std::make_pair(pathLength, desiredLength);
        auto it = rvCache.find(rvKey);
        if (it == rvCache.end()) {
            constructDeleteLength(pathLength, desiredLength, deleteLengthBuffer, deleteLengthWeight);
            return rvCache.insert({rvKey,
                discrete_dynamic(deleteLengthBuffer, deleteLengthWeight)}).first->second;
        }

        return it->second;
    }

    discrete_random_var_dynamic<int> getDeleteBeginRV(const RandomSequence<Vertex> &path,
            const int deleteLength, const int insertLength) {
        std::vector<int> &deleteBeginBuffer = mutation->deleteBegin;
        deleteBeginBuffer.clear();
        const Real specularRoughnessThreshold = mutation->specularRoughnessThreshold;
        const int pathLength = path.Size() + 1;
        int minIndex = 0;
        int maxIndex = pathLength - deleteLength - 1;
        if (deleteLength == 1 || insertLength == 1) {
            maxIndex--;
        }
        for (int index = minIndex; index <= maxIndex; index++) {
            const bool leftSpec = index > 0 && index < path.Size() - 1 &&
                                  isSpecular(path[index], specularRoughnessThreshold);
            const int rightIndex = index + deleteLength;
            const bool rightSpec = rightIndex > 0 && rightIndex < path.Size() - 1 &&
                                   isSpecular(path[rightIndex], specularRoughnessThreshold);
            if (!leftSpec && !rightSpec) {
                deleteBeginBuffer.push_back(index);
            }
        }

        return discrete_dynamic(deleteBeginBuffer);
    }

    discrete_random_var_dynamic<int> getEyeSubpathInsertSizeRV(const int insertSize) {
        auto &rvCache = mutation->eyeSubpathInsertSizeRV;
        auto it = rvCache.find(insertSize);
        if (it == rvCache.end()) {
            return rvCache.insert({insertSize, discrete_dynamic(insertSize + 1)}).first->second;
        }
        return it->second;
    }

};

void BidirectionalMutation::Mutate(const Scene *scene, const MarkovChainState &currentState,
        MarkovChainState &proposalState, OcclusionCache &occlusionCache) {
    statsAccepted.incrementBase(1);
    Node<bidir_mutation_t> strategy{this};
    Mutation::Mutate(strategy, scene, currentState, proposalState, occlusionCache);
}

double BidirectionalMutation::Pdf(const RandomSequence<Vertex> &target, const RandomSequence<Vertex> &given) {
    Node<bidir_mutation_t> strategy{this};
    return ConditionalPdf(strategy, target, given);
}

void BidirectionalMutation::Accepted() const {
    ++statsAccepted;
}

MTS_NAMESPACE_END
