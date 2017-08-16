#pragma once
#if !defined(__MITSUBA_RENDER_AETHER_BSDFS_IMPL_H_)
#define __MITSUBA_RENDER_AETHER_BSDFS_IMPL_H_

#include <mitsuba/render/aetherbsdfs.h>
#include <mitsuba/render/aethercommon.h>

MTS_NAMESPACE_BEGIN

template <typename DirectionToPrevious, typename ContextType>
auto cosine_hemisphere_sampling::Sample(const DirectionToPrevious&, 
                                        aether::Context<ContextType>& context,
                                        mitsuba::UniDist &uniDist) const {
    auto uv = context.Uniform2D(uniDist);
    return sample(CosHemisphereRV().Sample(uv[0], uv[1]));
}

inline auto BeckmannRV(const double alpha) {
    auto phi = two * get_pi() * u2;
    auto tanThetaSqr = constant(- alpha * alpha) * make_random_var(log(one - u1));
    auto cosTheta = make_random_var(one) / sqrt(make_random_var(one) + tanThetaSqr);
    auto sinTheta = sqrt(make_random_var(one) - sq(cosTheta));
    return make_random_vector(
            sinTheta * cos(phi),
            sinTheta * sin(phi),
            cosTheta
        );
}

template<typename Wi, typename M>
auto reflect(const Wi &wi, const M &m) {
    return two * dot(wi, m) * m - wi;
}

template<typename Wi, typename M>
auto refract(const Real eta, const Wi &wi, const M &m) {
    // http://www.graphics.cornell.edu/~bjw/microfacetbsdf.pdf
    // Eq 40, note that there is a typo and the eta inside the square root should be squared
    auto c = dot(wi, m);

    // eta = eta_interior / eta_exterior
    // eta_ = eta_incoming / eta_transmittance
    // if wi[2] > 0, we are at exterior, so eta_ = 1/eta, otherwise eta_ = eta
    auto eta_ = constant(at<2>(wi).Value() > Real(0) ? (Real(1) / eta) : eta);

    auto dsq = make_random_var(one) + sq(eta_) * (sq(c) - make_random_var(one));
    auto d = sqrt(dsq);
    auto e = pattern(when(at<2>(wi) > make_random_var(zero), d), otherwise(-d));
    auto ret = (eta_ * c - e) * m - eta_ * wi;
    return ret;
}

template <typename DirectionToPrevious, typename ContextType>
auto beckmann_sampling::Sample(const DirectionToPrevious &directionToPrevious, 
                               aether::Context<ContextType>& context,
                               mitsuba::UniDist &uniDist) const {
    auto uv = context.Uniform2D(uniDist);
    auto m = BeckmannRV(alpha).Sample(uv[0], uv[1]);

    auto dir = pattern(
        when(constant(Real(mode)) == constant(Real(REFLECT)), reflect(directionToPrevious, m)),
        otherwise(refract(eta, directionToPrevious, m))
    );
    return sample(dir);
}

MTS_NAMESPACE_END

#endif // __MITSUBA_RENDER_AETHER_BSDFS_IMPL_H_
