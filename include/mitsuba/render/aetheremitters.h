#pragma once
#if !defined(__MITSUBA_RENDER_AETHER_EMITTERS_H_)
#define __MITSUBA_RENDER_AETHER_EMITTERS_H_

#include <mitsuba/core/triangle.h>
#include <mitsuba/render/unidist.h>
#include <mitsuba/render/aethervariables.h>
#include <mitsuba/render/vertex.h>
#include <mitsuba/core/convertaether.h>
#include <aether/Context.h>
#include <aether/CompositeRandomVar.h>

// TODO: get rid of this
using namespace aether;

MTS_NAMESPACE_BEGIN

class Emitter;

namespace emitter_sampling_tag {
  struct Direct {};
  struct Position {};
  struct Direction {};
}

struct trimesh_emitter {
  trimesh_emitter() : trianglePtrs(nullptr), triangleWeights(nullptr), positions(nullptr) {}

  trimesh_emitter(const std::vector<const Triangle*> *trianglePtrs,
  			          const std::vector<aether::Real> *triangleWeights,
  				        const Point *positions,
                  const Normal *normals) :
    trianglePtrs(trianglePtrs), triangleWeights(triangleWeights), positions(positions), normals(normals) {}

  template <typename ContextType>
  auto SampleDirect(const Vertex &, aether::Context<ContextType>& context, mitsuba::UniDist &uniDist) const;

  template <typename ContextType>
  auto SamplePosition(aether::Context<ContextType>& context, mitsuba::UniDist &uniDist, const Float &time) const;

  template <typename ContextType>
  auto SampleDirection(const Vertex &, aether::Context<ContextType>& context, mitsuba::UniDist &uniDist) const;

  const std::vector<const Triangle*> *trianglePtrs;
  const std::vector<Real> *triangleWeights;
  const Point *positions;
  const Normal *normals;
};

struct envmap_emitter {
  envmap_emitter() : emitter(nullptr), distribution2d(nullptr), boundSphereRadius(0.f) {}

  envmap_emitter(const Emitter *emitter,
                 const Distribution2D *distribution2d,
                 const mitsuba::Point &boundSphereOrigin,
                 const Real boundSphereRadius) :
    emitter(emitter), distribution2d(distribution2d),
    boundSphereOrigin(boundSphereOrigin), boundSphereRadius(boundSphereRadius) {}

  template <typename ContextType>
  auto SampleDirect(const Vertex &, aether::Context<ContextType>& context, mitsuba::UniDist &uniDist) const;

  template <typename ContextType>
  auto SamplePosition(aether::Context<ContextType>& context, mitsuba::UniDist &uniDist, const Float &time) const;

  template <typename ContextType>
  auto SampleDirection(const Vertex &, aether::Context<ContextType>& context, mitsuba::UniDist &uniDist) const;

  const Emitter *emitter;
  const Distribution2D *distribution2d;
  mitsuba::Point boundSphereOrigin;
  aether::Real boundSphereRadius;
};

struct sphere_emitter {
  sphere_emitter() {}
  sphere_emitter(const mitsuba::Point &origin,
                 const aether::Real radius) : origin(origin), radius(radius) {}

  template <typename ContextType>
  auto SampleDirect(const Vertex &, aether::Context<ContextType>& context, mitsuba::UniDist &uniDist) const;

  template <typename ContextType>
  auto SamplePosition(aether::Context<ContextType>& context, mitsuba::UniDist &uniDist, const Float &time) const;

  template <typename ContextType>
  auto SampleDirection(const Vertex &, aether::Context<ContextType>& context, mitsuba::UniDist &uniDist) const;

  mitsuba::Point origin;
  aether::Real radius;
};

using emitter_composite_t = aether::CompositeRandomVar<trimesh_emitter, envmap_emitter, sphere_emitter>;

MTS_NAMESPACE_END

namespace aether {

template <>
struct SampleCall<mitsuba::trimesh_emitter> {
  template <int I, typename Cond, typename ContextType>
  auto operator()(const mitsuba::trimesh_emitter& emitter, _int<I>, Cond, mitsuba::emitter_sampling_tag::Direct, 
                  const mitsuba::Vertex &previousVertex, aether::Context<ContextType>& context, mitsuba::UniDist &uniDist) const {
    return emitter.SampleDirect(previousVertex, context, uniDist);
  }

  template <int I, typename Cond, typename ContextType>
  auto operator()(const mitsuba::trimesh_emitter& emitter, _int<I>, Cond, mitsuba::emitter_sampling_tag::Position,
                  aether::Context<ContextType>& context, mitsuba::UniDist &uniDist, const mitsuba::Float &time) const {
    return emitter.SamplePosition(context, uniDist, time);
  }

  template <int I, typename Cond, typename ContextType>
  auto operator()(const mitsuba::trimesh_emitter& emitter, _int<I>, Cond, mitsuba::emitter_sampling_tag::Direction,
                  const mitsuba::Vertex &previousVertex, aether::Context<ContextType>& context, mitsuba::UniDist &uniDist) const {
    return emitter.SampleDirection(previousVertex, context, uniDist);
  }
};

template <>
struct SampleCall<mitsuba::envmap_emitter> {
  template <int I, typename Cond, typename ContextType>
  auto operator()(const mitsuba::envmap_emitter& emitter, _int<I>, Cond, mitsuba::emitter_sampling_tag::Direct, 
                  const mitsuba::Vertex &previousVertex, aether::Context<ContextType>& context, mitsuba::UniDist &uniDist) const {
    return emitter.SampleDirect(previousVertex, context, uniDist);
  }

  template <int I, typename Cond, typename ContextType>
  auto operator()(const mitsuba::envmap_emitter& emitter, _int<I>, Cond, mitsuba::emitter_sampling_tag::Position,
                  aether::Context<ContextType>& context, mitsuba::UniDist &uniDist, const mitsuba::Float &time) const {
    return emitter.SamplePosition(context, uniDist, time);
  }

  template <int I, typename Cond, typename ContextType>
  auto operator()(const mitsuba::envmap_emitter& emitter, _int<I>, Cond, mitsuba::emitter_sampling_tag::Direction,
                  const mitsuba::Vertex &previousVertex, aether::Context<ContextType>& context, mitsuba::UniDist &uniDist) const {
    return emitter.SampleDirection(previousVertex, context, uniDist);
  }
};

template <>
struct SampleCall<mitsuba::sphere_emitter> {
  template <int I, typename Cond, typename ContextType>
  auto operator()(const mitsuba::sphere_emitter& emitter, _int<I>, Cond, mitsuba::emitter_sampling_tag::Direct, 
                  const mitsuba::Vertex &previousVertex, aether::Context<ContextType>& context, mitsuba::UniDist &uniDist) const {
    return emitter.SampleDirect(previousVertex, context, uniDist);
  }

  template <int I, typename Cond, typename ContextType>
  auto operator()(const mitsuba::sphere_emitter& emitter, _int<I>, Cond, mitsuba::emitter_sampling_tag::Position,
                  aether::Context<ContextType>& context, mitsuba::UniDist &uniDist, const mitsuba::Float &time) const {
    return emitter.SamplePosition(context, uniDist, time);
  }

  template <int I, typename Cond, typename ContextType>
  auto operator()(const mitsuba::sphere_emitter& emitter, _int<I>, Cond, mitsuba::emitter_sampling_tag::Direction,
                  const mitsuba::Vertex &previousVertex, aether::Context<ContextType>& context, mitsuba::UniDist &uniDist) const {
    return emitter.SampleDirection(previousVertex, context, uniDist);
  }
};

} // aether

#endif // __MITSUBA_RENDER_AETHER_EMITTERS_H_
