#pragma once
#if !defined(__MITSUBA_RENDER_VERTEX_H_)
#define __MITSUBA_RENDER_VERTEX_H_

#include <aether/RandomSequence.h>
#include <aether/Shapes.h>
#include <mitsuba/core/convertaether.h>
#include <mitsuba/render/aethervariables.h>
#include <mitsuba/render/shape.h>

MTS_NAMESPACE_BEGIN

// TODO: get rid of this
using namespace aether;

using element_data_t = decltype(make_named_struct(
  field<bool>(valid_)
  , field<aether::Vector3>(v0_)
  , field<aether::Vector3>(v1_)
  , field<aether::Vector3>(v2_)
  , field<aether::Vector3>(sph_center_)
  , field<Real>(sph_radius_)
  , field<Real>(is_triangle_)
  , field<mitsuba::Intersection>(intersection_)
  , field<const mitsuba::Emitter*>(emitter_)
  , field<const mitsuba::Sensor*>(sensor_)
));

using Vertex = Element<aether::Vector3, element_data_t>;

MTS_NAMESPACE_END

#endif /* __MITSUBA_RENDER_VERTEX_H_ */
