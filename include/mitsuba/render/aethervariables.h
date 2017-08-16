#pragma once
#if !defined(__MITSUBA_RENDER_AETHERVARIABLES_H_)
#define __MITSUBA_RENDER_AETHERVARIABLES_H_

#include <aether/NamedParameter.h>

MTS_NAMESPACE_BEGIN

// TODO: get rid of this
using namespace aether;

NAMED_PARAM(valid);
NAMED_PARAM(is_triangle);
NAMED_PARAM(intersection);
NAMED_PARAM(emitter);
NAMED_PARAM(sensor);
NAMED_PARAM(sph_center);
NAMED_PARAM(sph_radius);

MTS_NAMESPACE_END

#endif // __MITSUBA_RENDER_AETHERVARIABLES_H_
