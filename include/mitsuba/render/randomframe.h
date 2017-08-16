#pragma once
#if !defined(__MITSUBA_RENDER_RANDOM_FRAME_H_)
#define __MITSUBA_RENDER_RANDOM_FRAME_H_

#include <mitsuba/core/transform.h>
#include <mitsuba/render/shape.h>

MTS_NAMESPACE_BEGIN

template<typename S, typename T, typename N>
struct random_frame {
    S s;
    T t;
    N n;
};

template<typename S, typename T, typename N>
inline auto make_random_frame(const S &s, const T &t, const N &n) {
    return random_frame<S, T, N>{s, t, n};
}

inline auto make_random_frame(const Frame &frame) {
  return make_random_frame(constant(to_vector3(frame.s)),
                           constant(to_vector3(frame.t)),
                           constant(to_vector3(frame.n)));
}

inline auto make_random_frame(const Intersection &its) {
    return make_random_frame(its.shFrame);
}

inline auto make_random_frame(const Transform &trafo) {
    const Matrix4x4 &m = trafo.getMatrix();
    return make_random_frame(
        constant(aether::Vector3(m.m[0][0], m.m[1][0], m.m[2][0])),
        constant(aether::Vector3(m.m[0][1], m.m[1][1], m.m[2][1])),
        constant(aether::Vector3(m.m[0][2], m.m[1][2], m.m[2][2])));
}

inline auto make_random_frame(const mitsuba::Vector &N) {
  return make_random_frame(mitsuba::Frame(N));
}

template<typename N>
inline auto make_random_frame(const N &n) {
  auto a = make_random_vector(make_random_var(zero), -at<2>(n), at<1>(n));
  auto b = make_random_vector(at<2>(n), make_random_var(zero), -at<0>(n));
  auto s = normalize(pattern(
    when(dot(a, a) - dot(b, b) > make_random_var(zero), a)
    , otherwise(b)
  ));
  auto t = cross(n, s);
  return make_random_frame(s, t, n);
}

template<typename RandomFrame, typename LocalDir>
inline auto to_world(const RandomFrame &frame, const LocalDir &dir) {
    auto to_world_matrix = make_random_matrix(frame.s, frame.t, frame.n);
    return to_world_matrix * dir;
}

template<typename RandomFrame, typename WorldDir>
inline auto to_local(const RandomFrame &frame, const WorldDir &dir) {
    return make_random_vector(dot(frame.s, dir),
                              dot(frame.t, dir),
                              dot(frame.n, dir));
}

MTS_NAMESPACE_END

#endif // __MITSUBA_RENDER_RANDOM_FRAME_H_
