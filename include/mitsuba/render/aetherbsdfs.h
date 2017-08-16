#pragma once
#if !defined(__MITSUBA_RENDER_AETHER_BSDFS_H_)
#define __MITSUBA_RENDER_AETHER_BSDFS_H_

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

struct cosine_hemisphere_sampling {
  template <typename DirectionToPrevious, typename ContextType>
  auto Sample(const DirectionToPrevious&, aether::Context<ContextType>& context, mitsuba::UniDist &uniDist) const;

  inline bool operator==(const cosine_hemisphere_sampling &) const {
    return true;
  }
};

struct beckmann_sampling {
  enum Mode {
    REFLECT,
    REFRACT
  };

  template <typename DirectionToPrevious, typename ContextType>
  auto Sample(const DirectionToPrevious&, aether::Context<ContextType>& context, mitsuba::UniDist &uniDist) const;

  inline bool operator==(const beckmann_sampling &other) const {
    return alpha == other.alpha && eta == other.eta && mode == other.mode;
  }

  Real alpha;
  Real eta;
  Mode mode;
};

using bsdf_composite_t = aether::CompositeRandomVar<cosine_hemisphere_sampling, beckmann_sampling>;

struct bsdf_component_sampler {
  void AddComponent(const bsdf_composite_t &component, const Real weight) {
    // TODO: make this more efficient
    components.push_back(component);
    weights.push_back(weight);
  }

  template <typename DirectionToPrevious, typename ContextType>
  auto Sample(const DirectionToPrevious& directionToPrevious, aether::Context<ContextType>& context, mitsuba::UniDist &uniDist) const {
    if (components.size() == 0) {
      auto component = bsdf_composite_t{cosine_hemisphere_sampling{}};
      return component.Sample(directionToPrevious, context, uniDist);
    }
    auto componentSampler = discrete_dynamic(components, weights);
    auto component = context.Sample(componentSampler, context.Uniform1D(uniDist));
    return component.Sample(directionToPrevious, context, uniDist);
  }

  std::vector<bsdf_composite_t> components;
  std::vector<Real> weights;
};

MTS_NAMESPACE_END

namespace aether {

template <>
struct SampleCall<mitsuba::cosine_hemisphere_sampling> {
  template <int I, typename Cond, typename DirectionToPrevious, typename ContextType>
  auto operator()(const mitsuba::cosine_hemisphere_sampling& bsdf, _int<I>, Cond, const DirectionToPrevious &directionToPrevious,
                  aether::Context<ContextType>& context, mitsuba::UniDist &uniDist) const {
    return bsdf.Sample(directionToPrevious, context, uniDist);
  }
};

template <>
struct SampleCall<mitsuba::beckmann_sampling> {
  template <int I, typename Cond, typename DirectionToPrevious, typename ContextType>
  auto operator()(const mitsuba::beckmann_sampling& bsdf, _int<I>, Cond, const DirectionToPrevious &directionToPrevious,
                  aether::Context<ContextType>& context, mitsuba::UniDist &uniDist) const {
    return bsdf.Sample(directionToPrevious, context, uniDist);
  }
};

} // aether

#endif // __MITSUBA_RENDER_AETHER_EMITTERS_H_
