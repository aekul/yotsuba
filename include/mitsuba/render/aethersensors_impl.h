#pragma once
#if !defined(__MITSUBA_RENDER_AETHER_SENSORS_IMPL_H_)
#define __MITSUBA_RENDER_AETHER_SENSORS_IMPL_H_

#include <mitsuba/render/sensor.h>
#include <mitsuba/render/aethersensors.h>
#include <mitsuba/render/aethercommon.h>
#include <mitsuba/render/randomframe.h>
#include <mitsuba/render/randomtransform.h>

MTS_NAMESPACE_BEGIN

template <typename ContextType>
auto perspective_sensor::SamplePosition(aether::Context<ContextType>& context, mitsuba::UniDist &uniDist, const Float &time) const {
    const Transform trafo = sensor->getWorldTransform()->eval(time);
    auto pt = constant(to_vector3(trafo(Point(0.0f))));
    // Fill intersection information for integrand evaluation and coming up sampling
    Intersection its;
    its.p = to_point(pt.Value());
    //its.geoFrame = Frame(trafo(Vector(0.0f, 0.0f, 1.0f)));
    //its.shFrame = its.geoFrame;
    its.time = time;
    return sample(
      pt
      , intersection_ = its
    );
}

template <typename ContextType>
auto perspective_sensor::SampleDirection(const Vertex &previousVertex, aether::Context<ContextType>& context, mitsuba::UniDist &uniDist) const {
    SAssert(sensor != nullptr);
    const Float invAspect = Float(1) / sensor->getAspect();
    auto screen_x = make_random_var(u1) * constant(Real(-2) * tanFov) + constant(tanFov);
    auto screen_y = make_random_var(u2) * constant(Real(-2) * tanFov * invAspect) + constant(tanFov * invAspect);
    auto screen_z = make_random_var(one);
    auto uv = context.Uniform2D(uniDist);
    auto localDir = make_random_vector(screen_x, screen_y, screen_z).Sample(uv[0], uv[1]);
    const Transform trafo = sensor->getWorldTransform()->eval(previousVertex.Get(intersection_).time);
    auto frame = make_random_frame(trafo);
    auto worldDir = normalize(at<0>(localDir) * frame.s + at<1>(localDir) * frame.t + at<2>(localDir) * frame.n);
    return sample(worldDir);
}

template <typename ContextType>
auto perspective_sensor::SampleDirection(const Vertex &previousVertex,
        aether::Context<ContextType>& context, mitsuba::UniDist &uniDist, const int x, const int y) const {
    SAssert(sensor != nullptr);
    const Float invAspect = Float(1) / sensor->getAspect();
    const Vector2 invResolution = sensor->getInvResolution();

    auto screen_x = make_random_var(u1) * constant(Real(-2) * tanFov) + constant(tanFov);
    auto screen_y = make_random_var(u2) * constant(Real(-2) * tanFov * invAspect) + constant(tanFov * invAspect);
    auto screen_z = make_random_var(one);
    auto uv = context.Uniform2D(uniDist);
    // Hack: if both uvs are zero, we know it's from PdfContext, let them still be zero
    // context.PiecewiseConstant2D relies on this assumption so we need to preserve it.
    // A better approach is to return something that is not just a real number but also with a flag.
    // Maybe use the sign of uv?
    if (uv[0] != Real(0) || uv[1] != Real(0)) {
        uv[0] = (uv[0] + x) * invResolution.x;
        uv[1] = (uv[1] + y) * invResolution.y;
    }
    auto localDir = make_random_vector(screen_x, screen_y, screen_z).Sample(uv[0], uv[1]);
    const Transform trafo = sensor->getWorldTransform()->eval(previousVertex.Get(intersection_).time);
    auto frame = make_random_frame(trafo);
    auto worldDir = normalize(at<0>(localDir) * frame.s + at<1>(localDir) * frame.t + at<2>(localDir) * frame.n);
    return sample(worldDir);
}

template <typename ContextType, typename UV>
auto perspective_sensor::SampleDirection(const Vertex &previousVertex,
        aether::Context<ContextType>& context, mitsuba::UniDist &uniDist, const UV &uv) const {
    SAssert(sensor != nullptr);
    const Float invAspect = Float(1) / sensor->getAspect();

    auto screen_x = at<0>(uv) * constant(Real(-2) * tanFov) + constant(tanFov);
    auto screen_y = at<1>(uv) * constant(Real(-2) * tanFov * invAspect) + constant(tanFov * invAspect);
    auto screen_z = make_random_var(one);
    auto localDir = make_random_vector(screen_x, screen_y, screen_z);
    const Transform trafo = sensor->getWorldTransform()->eval(previousVertex.Get(intersection_).time);
    auto frame = make_random_frame(trafo);
    auto worldDir = normalize(at<0>(localDir) * frame.s + at<1>(localDir) * frame.t + at<2>(localDir) * frame.n);
    return sample(worldDir);
}

MTS_NAMESPACE_END

#endif // __MITSUBA_RENDER_AETHER_SENSORS_IMPL_H_
