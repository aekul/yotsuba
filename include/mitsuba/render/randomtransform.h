#pragma once
#if !defined(__MITSUBA_RENDER_RANDOM_TRANSFORM_H_)
#define __MITSUBA_RENDER_RANDOM_TRANSFORM_H_

#include <mitsuba/core/transform.h>
#include <mitsuba/render/shape.h>
#include <mitsuba/render/randomframe.h>

MTS_NAMESPACE_BEGIN

template<typename Frame, typename Translation, typename W, typename Woffset>
struct random_transform {
    Frame frame;
    Translation translation;
    W w;
    Woffset woffset;
};

template<typename Frame, typename Translation, typename W, typename Woffset>
inline auto make_random_transform(const Frame &f, const Translation &t, const W &w, const Woffset &wo) {
    return random_transform<Frame, Translation, W, Woffset>{f, t, w, wo};
}

inline auto make_random_transform(const Transform &trafo) {
    const Matrix4x4 &m = trafo.getMatrix();
    return make_random_transform(
        make_random_frame(trafo),
        constant(aether::Vector3(m.m[0][3], m.m[1][3], m.m[2][3])),
        constant(aether::Vector3(m.m[3][0], m.m[3][1], m.m[3][2])),
        constant(Real(m.m[3][3])));
}

template<typename RandomTransform, typename LocalPosition>
inline auto transform_position(const RandomTransform &transform, const LocalPosition &position) {
    auto frame = transform.frame;
    auto to_world_matrix = make_random_matrix(frame.s, frame.t, frame.n);
    auto homogenous_position = to_world_matrix * position + transform.translation;
    auto w = dot(transform.w, position) + transform.woffset;
    auto invW = rcp(w);
    return homogenous_position * invW;
}

MTS_NAMESPACE_END

#endif // __MITSUBA_RENDER_RANDOM_FRAME_H_
