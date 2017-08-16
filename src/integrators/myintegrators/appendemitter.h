#pragma once

#include <mitsuba/mitsuba.h>
#include <mitsuba/render/vertex.h>
#include <aether/RandomSequence.h>

MTS_NAMESPACE_BEGIN

struct UniDist;
class Emitter;
struct Raycast;

void AppendPositionOnEmitter(RandomSequence<Vertex> &path, UniDist &uniDist,
                             const ref_vector<Emitter> &emitters, const Float time);
void AppendDirectionFromEmitter(RandomSequence<Vertex> &path, UniDist &uniDist,
                                const Raycast &raycaster);
void AppendDirectSampleEmitter(RandomSequence<Vertex> &path, UniDist &uniDist,
                               const ref_vector<Emitter> &emitters);

MTS_NAMESPACE_END
