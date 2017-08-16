#include "appendsensor.h"
#include "intersect.h"
#include <mitsuba/render/scene.h>
#include <mitsuba/render/aethersensors_impl.h>

using namespace aether;

MTS_NAMESPACE_BEGIN

struct sample_position_on_sensor_t {
    template <typename T>
    auto operator()(Context<T>& context, const RandomSequence<Vertex>& path, const Ray& ray) const {
        return sample(constant(to_vector3(ray.o)));
    }

    template <typename T>
    auto operator()(Context<T>& context, const RandomSequence<Vertex>& path, UniDist& uniDist,
                    const Sensor *sensor, const Float &time) const {
        SAssert(sensor != nullptr);
        auto sampler = sensor->makeSampler();
        auto pt = sampler.Sample(sensor_sampling_tag::Position{}, context, uniDist, time);
        return sample(
            pt
            , sensor_ = sensor
        );
    }
};

struct sample_direction_from_sensor_t {
    template <typename T>
    auto operator()(Context<T>& context, const RandomSequence<Vertex>& path, const Raycast& raycaster,
                    const Ray& ray) const {
        // Use Mitsuba's raycasting engine to obtain intersection information
        const Intersection its = context.constant_call(raycaster, ray);
        // Do symbolic intersection
        auto intersectSample = intersect(constant(to_vector3(ray.o)), constant(to_vector3(ray.d)), its);
        // TODO: setup its.emitter
        const Emitter *emitter = its.shape != nullptr ? its.shape->getEmitter() : nullptr;

        return sample(
            intersectSample
            , emitter_ = emitter
        );
    }

    template <typename T>
    auto operator()(Context<T>& context, const RandomSequence<Vertex>& path, UniDist& uniDist,
                    const Raycast& raycaster, const int x, const int y) const {
        auto currentPos = constant(path.Back().Value());
        auto sensor = path.Back().Get(sensor_);
        SAssert(sensor != nullptr);
        auto sampler = sensor->makeSampler();
        auto dir = sampler.Sample(sensor_sampling_tag::Direction{}, path.Back(), context, uniDist, x, y);
        // Use Mitsuba's raycasting engine to obtain intersection information
        const Ray ray(to_point(currentPos.Value()), to_vector(dir.Value()),
                      path.Back().Get(intersection_).time);
        const Intersection its = context.constant_call(raycaster, ray);
        // Do symbolic intersection
        auto intersectSample = intersect(currentPos, dir, its);
        // TODO: setup its.emitter
        const Emitter *emitter = its.shape != nullptr ? its.shape->getEmitter() : nullptr;

        return sample(
            intersectSample
            , emitter_ = emitter
        );
    }

    template <typename T>
    auto operator()(Context<T>& context, const RandomSequence<Vertex>& path, UniDist& uniDist,
                    const Raycast& raycaster) const {
        auto currentPos = constant(path.Back().Value());
        auto sensor = path.Back().Get(sensor_);
        SAssert(sensor != nullptr);
        auto sampler = sensor->makeSampler();
        auto dir = sampler.Sample(sensor_sampling_tag::Direction{}, path.Back(), context, uniDist);
        // Use Mitsuba's raycasting engine to obtain intersection information
        const Ray ray(to_point(currentPos.Value()), to_vector(dir.Value()),
                      path.Back().Get(intersection_).time);
        const Intersection its = context.constant_call(raycaster, ray);
        // Do symbolic intersection
        auto intersectSample = intersect(currentPos, dir, its);
        // TODO: setup its.emitter
        const Emitter *emitter = its.shape != nullptr ? its.shape->getEmitter() : nullptr;

        return sample(
            intersectSample
            , emitter_ = emitter
        );
    }
};

void AppendPositionOnSensor(RandomSequence<Vertex> &path, UniDist &uniDist, const Sensor *sensor,
                            const Float time) {
    Node<sample_position_on_sensor_t> samplePositionOnSensor;
    path.Append(samplePositionOnSensor, uniDist, sensor, time);
}

void AppendPositionOnSensor(RandomSequence<Vertex> &path, const Ray &ray) {
    Node<sample_position_on_sensor_t> samplePositionOnSensor;
    path.Append(samplePositionOnSensor, ray);   
}

void AppendDirectionFromSensor(RandomSequence<Vertex> &path, UniDist &uniDist, const Raycast &raycaster) {
    Node<sample_direction_from_sensor_t> sampleDirectionFromSensor;
    path.Append(sampleDirectionFromSensor, uniDist, raycaster);
}

void AppendDirectionFromSensor(RandomSequence<Vertex> &path, UniDist &uniDist, const Raycast &raycaster, const int x, const int y) {
    Node<sample_direction_from_sensor_t> sampleDirectionFromSensor;
    path.Append(sampleDirectionFromSensor, uniDist, raycaster, x, y);
}

void AppendDirectionFromSensor(RandomSequence<Vertex> &path, const Raycast &raycaster, const Ray &ray) {
    Node<sample_direction_from_sensor_t> sampleDirectionFromSensor;
    path.Append(sampleDirectionFromSensor, raycaster, ray); 
}

MTS_NAMESPACE_END
