#pragma once
#if !defined(__MITSUBA_RENDER_AETHER_SENSORS_H_)
#define __MITSUBA_RENDER_AETHER_SENSORS_H_

#include <mitsuba/render/unidist.h>
#include <mitsuba/render/aethervariables.h>
#include <mitsuba/render/vertex.h>
#include <mitsuba/core/convertaether.h>
#include <aether/Context.h>

// TODO: get rid of this
using namespace aether;

MTS_NAMESPACE_BEGIN

class Sensor;

namespace sensor_sampling_tag {
  //struct Direct {}; // TODO: implement
  struct Position {};
  struct Direction {};
}

struct perspective_sensor {
  perspective_sensor() : sensor(nullptr), tanFov(0.f) {}

  perspective_sensor(const Sensor *sensor, const Float tanFov) :
    sensor(sensor), tanFov(tanFov) {}

  template <typename ContextType>
  auto SamplePosition(aether::Context<ContextType>& context, mitsuba::UniDist &uniDist, const Float &time) const;

  template <typename ContextType>
  auto SampleDirection(const Vertex &, aether::Context<ContextType>& context, mitsuba::UniDist &uniDist) const;

  template <typename ContextType>
  auto SampleDirection(const Vertex &, aether::Context<ContextType>& context, mitsuba::UniDist &uniDist, const int x, const int y) const;

  template <typename ContextType, typename UV>
  auto SampleDirection(const Vertex &, aether::Context<ContextType>& context, mitsuba::UniDist &uniDist, const UV &uv) const;

  // placeholders until we have another sensor
  template <typename ContextType>
  auto Sample(sensor_sampling_tag::Position, aether::Context<ContextType>& context, mitsuba::UniDist &uniDist, const Float &time) const {
    return SamplePosition(context, uniDist, time);
  }

  template <typename ContextType>
  auto Sample(sensor_sampling_tag::Direction, const Vertex &previousVertex, aether::Context<ContextType>& context, mitsuba::UniDist &uniDist) const {
    return SampleDirection(previousVertex, context, uniDist);
  }

  template <typename ContextType>
  auto Sample(sensor_sampling_tag::Direction, const Vertex &previousVertex,
      aether::Context<ContextType>& context, mitsuba::UniDist &uniDist, const int x, const int y) const {
    return SampleDirection(previousVertex, context, uniDist, x, y);
  }

  template <typename ContextType, typename UV>
  auto Sample(sensor_sampling_tag::Direction, const Vertex &previousVertex,
      aether::Context<ContextType>& context, mitsuba::UniDist &uniDist, const UV &uv) const {
    return SampleDirection(previousVertex, context, uniDist, uv);
  }

  const Sensor *sensor;
  Float tanFov;
};

using sensor_composite_t = perspective_sensor;

MTS_NAMESPACE_END

namespace aether {

// TODO: define SampleCall

} // aether

#endif // __MITSUBA_RENDER_AETHER_SENSORS_H_
