#include "appendkeyhole.h"
#include "intersect.h"
#include <mitsuba/render/scene.h>
#include <mitsuba/render/randomframe.h>

#include <aether/UniformTriangle.h>

MTS_NAMESPACE_BEGIN

template<int N, typename T, typename Org, typename Dir>
auto intersect(Context<T> &context, const Raycast &raycaster,
			   const Org &org, const Dir &dir, const float time) {
	// intersect one side of the direction with the scene
	Ray ray{to_point(org.Value()), to_vector(dir.Value()), time};
	const Intersection its = context.constant_call(raycaster, ray);
	auto its_v0 = to_vector3(its.triangle_v0);
	auto its_v1 = to_vector3(its.triangle_v1);
	auto its_v2 = to_vector3(its.triangle_v2);
	auto sph_center = to_vector3(its.sph_center);
	auto sph_radius = aether::Real(its.sph_radius);
	auto pt = intersect_triangle<N>(its_v0, its_v1, its_v2, org, dir);

	const Emitter *emitter = its.shape != nullptr ? its.shape->getEmitter() : nullptr;

	return optional_sample(
		its.isValid()
		, pt
		, v0_ = its_v0
		, v1_ = its_v1
		, v2_ = its_v2
		, sph_center_ = sph_center
		, sph_radius_ = sph_radius
		, is_triangle_ = its.isTriangle
		, intersection_ = its
		, emitter_ = emitter
	);
}

// Sample a position on a triangle mesh (representing a keyhole), extend the path to generate two vertices
struct sample_keyhole_t {
	template <typename T>
	auto operator()(Context<T>& context, const RandomSequence<Vertex>& path, UniDist& uniDist,
					const Raycast& raycaster, const std::vector<std::array<aether::Vector3, 3>>& tris) const {
	
		// discrete_dynamic(tris) creates a symbolic discrete random variable
		// note that we uniformly sample the triangles here, which ignores the area of the triangles
		// it is possible to pass an optional weight array to discrete_dynamic to sample according to areas
		auto tris_rv = discrete_dynamic(tris);
		// sample a triangle
		auto tri = context.Sample(tris_rv, context.Uniform1D(uniDist));

		// sample a point on the triangle
		auto v0 = constant(tri[0]);
		auto v1 = constant(tri[1]);
		auto v2 = constant(tri[2]);
		auto uv = context.Uniform2D(uniDist);
		auto keyhole_pt = uniform_triangle(v0, v1, v2).Sample(uv[0], uv[1]);

		// sample a direction wi on a unit sphere
		// note that we need to use a new set of uniforms (u3, u4)
		auto r = sqrt(aether::one - sq(u3));
		auto phi = two * get_pi() * u4;
		auto x = make_random_var(r * cos(phi));
		auto y = make_random_var(r * sin(phi));
		auto z = make_random_var(u3);
		auto ab = context.Uniform2D(uniDist);
		auto dir_local = make_random_vector(x, y, z).Sample(ab[0], ab[1]);

		// transform the direction to world space
		auto e1 = v0 - v2;
		auto e2 = v1 - v2;
		auto N = normalize(cross(e1, e2));
		auto basis = coordinate_basis(N);
		auto dir_world = basis * dir_local;

		Float time = context.Uniform1D(uniDist);

		// intersect the two sides of the direction with the scene
		auto sample0 = intersect<0>(context, raycaster, keyhole_pt,  dir_world, time);
		auto sample1 = intersect<1>(context, raycaster, keyhole_pt, -dir_world, time);

		// sample_tuple combines the two path vertices
		return sample_tuple(sample0, sample1);
	}
};

void AppendKeyhole(RandomSequence<Vertex> &path, UniDist& uniDist,
				   const Raycast& raycaster, const std::vector<std::array<aether::Vector3, 3>>& tris) {
	Node<sample_keyhole_t> sampleKeyhole;
	path.Append(sampleKeyhole, uniDist, raycaster, tris);
}

MTS_NAMESPACE_END
