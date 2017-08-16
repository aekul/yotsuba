#pragma once

#include <mitsuba/render/aethervariables.h>
#include <mitsuba/render/aethercommon.h>
#include <aether/RandomVar.h>

// TODO: get rid of this
using namespace aether;

MTS_NAMESPACE_BEGIN

template<int ID, typename IsTriangle, typename Org, typename Dir>
auto intersect(const aether::Vector3 &v0,
               const aether::Vector3 &v1,
               const aether::Vector3 &v2,
               const aether::Vector3 &sph_center,
               const aether::Real sph_radius,
               const IsTriangle &is_triangle,
               const Org &org,
               const Dir &dir) {
    // If we do not hit a triangle, intersect with a sphere
    return pattern(
        when(is_triangle == make_random_var(zero),
            intersect_sphere(
                named_constant<ID>(sph_center_, sph_center),
                named_constant<ID>(sph_radius_, sph_radius),
                org, dir
            )
        ),
        otherwise(
            intersect_triangle<ID>(v0, v1, v2, org, dir)
        )
    );
}

template <typename Org, typename Dir>
auto intersect(const Org &org, const Dir &dir, const Intersection &nextIts) {
    auto v0 = to_vector3(nextIts.triangle_v0);
    auto v1 = to_vector3(nextIts.triangle_v1);
    auto v2 = to_vector3(nextIts.triangle_v2);
    auto sph_center = to_vector3(nextIts.sph_center);
    auto sph_radius = aether::Real(nextIts.sph_radius);
    auto is_tri = named_constant<0>(is_triangle_, Real(nextIts.isTriangle));
    auto pt = intersect<0>(v0, v1, v2, sph_center, sph_radius, is_tri, org, dir);
    return sample(
        pt
        , v0_ = v0
        , v1_ = v1
        , v2_ = v2
        , sph_center_ = sph_center
        , sph_radius_ = sph_radius
        , is_triangle_ = nextIts.isTriangle
        , intersection_ = nextIts
    );
}

MTS_NAMESPACE_END
