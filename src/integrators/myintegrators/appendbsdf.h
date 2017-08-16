#pragma once

#include <mitsuba/mitsuba.h>
#include <mitsuba/render/vertex.h>
#include <mitsuba/render/bsdf.h>
#include <aether/RandomSequence.h>

MTS_NAMESPACE_BEGIN

struct UniDist;
struct Raycast;

void AppendBSDF(RandomSequence<Vertex> &path, UniDist &uniDist, const Raycast &raycaster,
                const int rrDepth = -1, const uint32_t componentMask = BSDF::EAll);

MTS_NAMESPACE_END
