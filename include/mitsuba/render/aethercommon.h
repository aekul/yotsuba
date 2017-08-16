#pragma once
#if !defined(__MITSUBA_RENDER_AETHER_COMMON_H_)
#define __MITSUBA_RENDER_AETHER_COMMON_H_

#include <aether/RandomVar.h>

MTS_NAMESPACE_BEGIN

template <typename SO, typename SR, typename Org, typename Dir>
constexpr auto intersect_sphere(const SO &so, // sphere origin
                                const SR &sr, // sphere radius
                                const Org &org,
                                const Dir &dir) {
    // intersect ray(o, d) with sphere(so, sr)
    // Solve t^2*d.d + 2*t*(so-o).d + (so-o).(so-o) = SR^2 
    auto od = so - org;
    auto proj = dot(od, dir);
    auto det = sqrt(sq(proj) - dot(od, od) + sq(sr));
    auto t = pattern(
        when(length(od) > sr, proj - det)
        , otherwise(proj + det)
    );
    auto pt = org + dir * t;
    return pt;
}

constexpr auto CosHemisphereRV() {
	constexpr auto r = sqrt(one - u1);
	constexpr auto phi = two * get_pi() * u2;
	constexpr auto x = r * cos(phi);
	constexpr auto y = r * sin(phi);
	constexpr auto z = sqrt(u1);
	return make_random_vector(make_random_var(x), make_random_var(y), make_random_var(z));
}

MTS_NAMESPACE_END

#endif // __MITSUBA_RENDER_AETHER_COMMON_H_
