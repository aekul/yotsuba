#pragma once

#include <mitsuba/render/randomframe.h>

MTS_NAMESPACE_BEGIN

struct perturb_direction {
    const Real theta_1;
    const Real theta_2;
    UniDist &uniDist;
    const Raycast &raycaster;

    template <typename T>
    auto operator()(Context<T>& context,
                    const aether::Vector3 &previousPosition,
                    const aether::Vector3 &previousDirection,
                    const Float time) {
        // Generate an offset direction
        auto theta = constant(theta_2) * exp(-log(constant(theta_2 / theta_1)) * u1);
        auto phi = two * get_pi() * u2;
        auto sin_theta = sin(theta);
        auto cos_theta = cos(theta);
        auto sin_phi = sin(phi);
        auto cos_phi = cos(phi);
        auto offset = make_random_vector(
            sin_theta * cos_phi,
            sin_theta * sin_phi,
            cos_theta);
        auto frame = make_random_frame(to_vector(previousDirection));
        auto uv = context.Uniform2D(uniDist);
        auto p = constant(previousPosition);
        auto dir = to_world(frame, offset).Sample(uv[0], uv[1]);

        // Use Mitsuba's raycasting engine to obtain intersection information
        const Ray ray(to_point(p.Value()), to_vector(dir.Value()), time);
        const Intersection its = context.constant_call(raycaster, ray);
        // Do symbolic intersection
        auto intersectSample = intersect(p, dir, its);

        // TODO: setup its.emitter
        const Emitter *emitter = its.shape != nullptr ? its.shape->getEmitter() : nullptr;

        return sample(
            intersectSample
            , emitter_ = emitter
        );
    }
};

MTS_NAMESPACE_END
