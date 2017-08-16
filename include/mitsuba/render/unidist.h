#pragma once
#if !defined(__MITSUBA_RENDER_UNIDIST_H_)
#define __MITSUBA_RENDER_UNIDIST_H_

#include <mitsuba/render/sampler.h>
#include <aether/Sampler.h>

MTS_NAMESPACE_BEGIN

struct UniDist : public aether::BaseSampler {
  UniDist(mitsuba::Sampler* sampler) 
    : sampler(sampler)
  {}

  aether::Real operator()() override {
    return Uniform1D();
  }

  aether::Real Uniform1D(int) override {
    return sampler->next1D();
  }

  aether::Real Uniform1D() override {
    return sampler->next1D();
  }

  std::array<aether::Real, 2> Uniform2D(int) override {
    Point2 uv = sampler->next2D();
    return std::array<aether::Real, 2>{{uv[0], uv[1]}};
  }

  std::array<aether::Real, 2> Uniform2D() override {
    Point2 uv = sampler->next2D();
    return std::array<aether::Real, 2>{{uv[0], uv[1]}};
  }

  mitsuba::Sampler* sampler;
};

MTS_NAMESPACE_END

#endif // __MITSUBA_RENDER_UNIDIST_H_
