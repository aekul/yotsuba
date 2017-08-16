#pragma once

#include <mitsuba/mitsuba.h>
#include <mitsuba/render/vertex.h>
#include <aether/RandomSequence.h>

MTS_NAMESPACE_BEGIN

bool isSpecular(const Vertex &vertex, const Float specularRoughnessThreshold);
bool isDiffuse(const Vertex &vertex, const Float specularRoughnessThreshold);

MTS_NAMESPACE_END
