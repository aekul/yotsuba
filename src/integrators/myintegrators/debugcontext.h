#ifndef AETHER_CONTEXT_DEBUG_CONTEXT_H
#define AETHER_CONTEXT_DEBUG_CONTEXT_H

#include <aether/Context.h>
#include <aether/RandomVar.h>
#include <aether/Strategy.h>

namespace aether {

struct DebugContext : Context<DebugContext> {
  template <typename Fn, typename... Args>
  constexpr auto constant_call(Fn&& fn, Args&&... args);

  template <typename T, typename... Args>
  auto Sample(const discrete_random_var_dynamic<T>& rv, Args&&... args);

  template <typename Sampler>
  Real Uniform1D(Sampler& sampler);

  template <typename Sampler>
  std::array<Real, 2> Uniform2D(Sampler& sampler);

  template <typename Sampler>
  std::array<Real, 2> PiecewiseConstant2D(Sampler& sampler, const Distribution2D *dist2D);
};

inline DebugContext make_debug_context() {
  return {};
}

template <typename Fn, typename... Args>
constexpr auto DebugContext::constant_call(Fn&& fn, Args&&... args) {
  using result_t = decltype(std::forward<Fn>(fn)(std::forward<Args>(args)...)); 
  return result_t{};
}

template <typename Sampler>
Real DebugContext::Uniform1D(Sampler& sampler) {
  return 0;
}

template <typename Sampler>
std::array<Real, 2> DebugContext::Uniform2D(Sampler& sampler) {
  return std::array<Real, 2>{{0.0, 0.0}};
}

template <typename Sampler>
std::array<Real, 2> DebugContext::PiecewiseConstant2D(Sampler& sampler, const Distribution2D *dist2D) {
  return std::array<Real, 2>{{0.0, 0.0}};
}

template <typename T, typename... Args>
auto DebugContext::Sample(const discrete_random_var_dynamic<T>& rv, Args&&... args) {
  return rv.Sample(std::forward<Args>(args)...);
}

} // end namespace aether

#endif
