#pragma once
#if !defined(__MITSUBA_RENDER_AETHER_EMITTERS_IMPL_H_)
#define __MITSUBA_RENDER_AETHER_EMITTERS_IMPL_H_

#include <mitsuba/core/convertaether.h>
#include <mitsuba/render/aetheremitters.h>
#include <mitsuba/render/aethercommon.h>
#include <mitsuba/render/randomframe.h>
#include <mitsuba/render/emitter.h>
#include <aether/UniformTriangle.h>

MTS_NAMESPACE_BEGIN

template <typename ContextType>
auto trimesh_emitter::SampleDirect(const Vertex &previousVertex, aether::Context<ContextType>& context, mitsuba::UniDist &uniDist) const {
	return SamplePosition(context, uniDist, previousVertex.Get(intersection_).time);
}

template <typename ContextType>
auto trimesh_emitter::SamplePosition(aether::Context<ContextType>& context, mitsuba::UniDist &uniDist, const Float &time) const {
	// TODO: deal with moving triangle meshes
	SAssert(trianglePtrs != nullptr);
	SAssert(triangleWeights != nullptr);
	auto rv = discrete_dynamic(*trianglePtrs, *triangleWeights);
	auto trianglePtr = context.Sample(rv, context.Uniform1D(uniDist));
	const aether::Vector3 v0 = to_vector3(positions[trianglePtr->idx[0]]);
	const aether::Vector3 v1 = to_vector3(positions[trianglePtr->idx[1]]);
	const aether::Vector3 v2 = to_vector3(positions[trianglePtr->idx[2]]);
	auto uv = context.Uniform2D(uniDist);
	auto pt = uniform_triangle(v0, v1, v2).Sample(uv[0], uv[1]);
	// Fill intersection information for integrand evaluation and coming up sampling
	Intersection its;
	// Important to show that this intersection is valid
	its.t = 0.f;
	its.p = to_point(pt.Value());
	its.geoFrame = Frame(to_normal((v1 - v0).cross(v2 - v0).normalized()));
	if (normals == nullptr) {
		its.shFrame = its.geoFrame;
	} else {
		// Do phong normal interpolation
		const Normal n0 = normals[trianglePtr->idx[0]];
		const Normal n1 = normals[trianglePtr->idx[1]];
		const Normal n2 = normals[trianglePtr->idx[2]];
		// We can also project its.p back to its barycentric coordinate, which is less hacky than the following
    	const Float a = sqrt(uv[0]);
    	const Float b0 = 1.f - a;
	    const Float b1 = uv[1] * a;
		its.shFrame = Frame(b0 * n0 + b1 * n1 + (1.f - b0 - b1) * n2);
	}
	its.time = time;
	return sample(
	  pt
	  , v0_ = v0
	  , v1_ = v1
	  , v2_ = v2
	  , is_triangle_ = Real(true)
	  , intersection_ = its
	);
}

template <typename ContextType>
auto trimesh_emitter::SampleDirection(const Vertex &previousVertex, aether::Context<ContextType>& context, mitsuba::UniDist &uniDist) const {
	auto uv = context.Uniform2D(uniDist);
	auto localDir = CosHemisphereRV().Sample(uv[0], uv[1]);
	auto frame = make_random_frame(previousVertex.Get(intersection_));
	auto worldDir = to_world(frame, localDir);
	return sample(worldDir);
}

constexpr auto SphericalRV() {
	constexpr auto phi = two * get_pi() * u1;
	constexpr auto sinPhi = sin(phi);
	constexpr auto cosPhi = cos(phi);
	constexpr auto theta = get_pi() * u2;
	constexpr auto sinTheta = sin(theta);
	constexpr auto cosTheta = cos(theta);
	return make_random_vector(
	        make_random_var(sinPhi * sinTheta),
	        make_random_var(cosTheta),
	        make_random_var(- cosPhi * sinTheta)
	      );
}

template <typename ContextType>
auto envmap_emitter::SampleDirect(const Vertex &previousVertex, 
	          aether::Context<ContextType>& context, mitsuba::UniDist &uniDist) const {
	SAssert(distribution2d != nullptr);
	SAssert(emitter != nullptr);
	auto uv = context.PiecewiseConstant2D(uniDist, distribution2d);
	const Float time = previousVertex.Get(intersection_).time;
	const Transform &trafo =
		emitter->getWorldTransform()->eval(time);
	auto frame = make_random_frame(trafo);
	auto dir = to_world(frame, (SphericalRV().Sample(uv[0], uv[1])));
	auto p = constant(previousVertex.Value());
	auto bo = named_constant<0>(sph_center_, to_vector3(boundSphereOrigin));
	auto br = named_constant<0>(sph_radius_, boundSphereRadius);
	auto pt = intersect_sphere(bo, br, p, dir);

	Intersection its;
	its.isTriangle = false;
	its.sph_center = boundSphereOrigin;
	its.sph_radius = boundSphereRadius;
	its.time = time;
	return sample(
	  pt
	  , sph_center_ = to_vector3(boundSphereOrigin)
	  , sph_radius_ = boundSphereRadius
	  , is_triangle_ = Real(false)
	  , intersection_ = its
	);
}

template <typename ContextType>
auto envmap_emitter::SamplePosition(aether::Context<ContextType>& context, mitsuba::UniDist &uniDist, const Float &time) const {
	SAssert(emitter != nullptr);
	auto uv = context.Uniform2D(uniDist);
	auto d = SphericalRV().Sample(uv);
	auto bo = named_constant<0>(sph_center_, to_vector3(boundSphereOrigin));
	auto br = named_constant<0>(sph_radius_, boundSphereRadius);
	auto pt = bo + d * br;

	Intersection its;
	its.isTriangle = false;
	its.sph_center = boundSphereOrigin;
	its.sph_radius = boundSphereRadius;
	its.time = time;
	return sample(
	  pt
	  , sph_center_ = to_vector3(boundSphereOrigin)
	  , sph_radius_ = boundSphereRadius
	  , is_triangle_ = Real(false)
	  , intersection_ = its
	);
}

template <typename ContextType>
auto envmap_emitter::SampleDirection(const Vertex &previousVertex, aether::Context<ContextType>& context, mitsuba::UniDist &uniDist) const {
#if 1
	// Should be this but it doesn't compile
	SAssert(distribution2d != nullptr);
	SAssert(emitter != nullptr);
	//auto uv = context.PiecewiseConstant2D(uniDist, distribution2d);
  auto uv = context.Uniform2D(uniDist);
	const Transform &trafo =
		emitter->getWorldTransform()->eval(previousVertex.Get(intersection_).time);
	auto frame = make_random_frame(trafo);
	auto dir = to_world(frame, -(SphericalRV().Sample(uv[0], uv[1])));
	return sample(dir);
#endif
	//auto uv = context.Uniform2D(uniDist);
	//auto localDir = CosHemisphereRV().Sample(uv[0], uv[1]);
	//auto frame = make_random_frame(previousVertex.Get(intersection_));
	//auto worldDir = to_world(frame, localDir);
	//return sample(worldDir);
}

template<typename CosCutoffMinusOne>
auto ConeRV(CosCutoffMinusOne &cosCutOffMinusOne) {
    auto u1r = make_random_var(u1);
    auto cosTheta = make_random_var(one) + u1r * cosCutOffMinusOne;
    auto sinTheta = sqrt(one - sq(cosTheta));
    auto phi = two * get_pi() * u2;
    auto sinPhi = sin(phi);
    auto cosPhi = cos(phi);
    return make_random_vector(
            cosPhi * sinTheta,
            sinPhi * sinTheta,
            cosTheta
        );
}

template <typename ContextType>
auto sphere_emitter::SampleDirect(const Vertex &previousVertex, 
	          aether::Context<ContextType>& context, mitsuba::UniDist &uniDist) const {
	auto p = constant(previousVertex.Value());
	auto so = named_constant<0>(sph_center_, to_vector3(origin));
	auto sr = named_constant<0>(sph_radius_, radius);
	auto toCenter = to_vector3(origin) - previousVertex.Value();
	auto distance = toCenter.norm();
	auto sinAlpha = radius / distance;
	auto cosAlpha = sqrt(1.f - sinAlpha * sinAlpha);
	auto toCenterN = toCenter / distance;
	auto frame = make_random_frame(to_vector(toCenterN));
	auto cosCutoffMinusOne = constant(cosAlpha - 1.0);//cosAlpha - make_random_var(one);
	auto uv = context.Uniform2D(uniDist);
	auto dir = to_world(frame, (ConeRV(cosCutoffMinusOne).Sample(uv[0], uv[1])));
	auto pt = intersect_sphere(so, sr, p, dir);
	// Fill intersection information for integrand evaluation
	Intersection its;
	its.isTriangle = false;
	its.sph_center = origin;
	its.sph_radius = radius;
	its.geoFrame = Frame(to_normal((pt.Value() - to_vector3(origin)).normalized()));
	its.shFrame = its.geoFrame;
	its.t = (Float)length(pt.Value() - previousVertex.Value());
	its.time = previousVertex.Get(intersection_).time;
	return sample(
	  pt
	  , sph_center_ = to_vector3(origin)
	  , sph_radius_ = radius
	  , is_triangle_ = Real(false)
	  , intersection_ = its
	);
}

template <typename ContextType>
auto sphere_emitter::SamplePosition(aether::Context<ContextType>& context, mitsuba::UniDist &uniDist, const Float &time) const {
	auto uv = context.Uniform2D(uniDist);
	auto d = SphericalRV().Sample(uv);
	auto bo = named_constant<0>(sph_center_, to_vector3(origin));
	auto br = named_constant<0>(sph_radius_, radius);
	auto pt = bo + br * d;

	Intersection its;
	its.isTriangle = false;
	its.sph_center = origin;
	its.sph_radius = radius;
	its.geoFrame = Frame(to_normal(d.Value()));
	its.shFrame = its.geoFrame;
	// Important to show that this intersection is not invalid
	its.t = Float(0);
	its.time = time;
	return sample(
	  pt
	  , sph_center_ = to_vector3(origin)
	  , sph_radius_ = radius
	  , is_triangle_ = Real(false)
	  , intersection_ = its
	);
}

template <typename ContextType>
auto sphere_emitter::SampleDirection(const Vertex &previousVertex, aether::Context<ContextType>& context, mitsuba::UniDist &uniDist) const {
	auto uv = context.Uniform2D(uniDist);
	auto localDir = CosHemisphereRV().Sample(uv[0], uv[1]);
	auto frame = make_random_frame(previousVertex.Get(intersection_));
	auto worldDir = to_world(frame, localDir);
	return sample(worldDir);
}

MTS_NAMESPACE_END

#endif // __MITSUBA_RENDER_AETHER_EMITTERS_IMPL_H_
