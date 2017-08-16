#pragma once

#include <aether/fwd/Math.h>

MTS_NAMESPACE_BEGIN

#if !defined(__MITSUBA_CORE_CONVERT_AETHER_H_)
#define __MITSUBA_CORE_CONVERT_AETHER_H_

inline aether::Vector2 to_vector2(const mitsuba::Point2 &x) {
	return aether::Vector2(x[0], x[1]);
}

inline aether::Vector2 to_vector2(const mitsuba::Vector2 &x) {
	return aether::Vector2(x[0], x[1]);
}

inline aether::Vector3 to_vector3(const mitsuba::Point &x) {
  return aether::Vector3(x[0], x[1], x[2]);
}

inline aether::Vector3 to_vector3(const mitsuba::Vector &x) {
  return aether::Vector3(x[0], x[1], x[2]);
}

inline aether::Vector3 to_vector3(const mitsuba::Normal &x) {
  return aether::Vector3(x[0], x[1], x[2]);
}

inline Point to_point(const aether::Vector3 &x) {
	return Point(x[0], x[1], x[2]);
}

inline Vector to_vector(const aether::Vector3 &x) {
	return Vector(x[0], x[1], x[2]);
}

inline Normal to_normal(const aether::Vector3 &x) {
	return Normal(x[0], x[1], x[2]);
}


MTS_NAMESPACE_END

#endif // __MITSUBA_CORE_CONVERT_AETHER_H_
