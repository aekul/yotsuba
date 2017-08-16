#pragma once

#include <mitsuba/render/scene.h>
#include <vector>

MTS_NAMESPACE_BEGIN

struct MarkovChainState;

Float sampleSeedPaths(const Scene *scene, Sampler *sampler,
    const int luminanceSamples, const int workUnits,
    const int maxDepth, const bool separateDirect,
    std::vector<MarkovChainState> &seeds);

MTS_NAMESPACE_END
