#pragma once

#include <mitsuba/mitsuba.h>
#include <mitsuba/render/vertex.h>
#include <mitsuba/render/bsdf.h>
#include <aether/RandomSequence.h>

MTS_NAMESPACE_BEGIN

struct UniDist;
struct Raycast;

void AppendKeyhole(RandomSequence<Vertex> &path, UniDist &uniDist, const Raycast &raycaster,
                   const std::vector<std::array<aether::Vector3, 3>>& tris);

MTS_NAMESPACE_END
