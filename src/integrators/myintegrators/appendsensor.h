#pragma once

#include <mitsuba/mitsuba.h>
#include <mitsuba/render/vertex.h>
#include <aether/RandomSequence.h>

MTS_NAMESPACE_BEGIN

struct UniDist;
class Sensor;
struct Raycast;

void AppendPositionOnSensor(RandomSequence<Vertex> &path, UniDist &uniDist, const Sensor *sensor,
                            const Float time);
void AppendPositionOnSensor(RandomSequence<Vertex> &path, const Ray &ray);
void AppendDirectionFromSensor(RandomSequence<Vertex> &path, UniDist &uniDist, const Raycast &raycaster);
void AppendDirectionFromSensor(RandomSequence<Vertex> &path, UniDist &uniDist, const Raycast &raycaster,
                               const int x, const int y);
void AppendDirectionFromSensor(RandomSequence<Vertex> &path, const Raycast &raycaster, const Ray &ray);

MTS_NAMESPACE_END
