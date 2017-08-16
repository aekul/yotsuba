#include "appendemitter.h"
#include "intersect.h"
#include <mitsuba/render/scene.h>
#include <mitsuba/render/aetheremitters_impl.h>

MTS_NAMESPACE_BEGIN

// sample an emitter and sample a point on the emitter using direct importance sampling
struct direct_sample_emitter_t {
    template <typename T>
    auto operator()(Context<T>& context, const RandomSequence<Vertex>& path,
                    UniDist& uniDist, const ref_vector<Emitter> &emitters) const {
        // discrete_dynamic(emitters) creates a symbolic discrete random variable
        auto emitters_rv = discrete_dynamic(emitters);
        // context.Sample samples an emitter
        auto emitter = context.Sample(emitters_rv, context.Uniform1D(uniDist));
        // create a sampler for generating points on light source
        auto sampler = emitter->makeSampler();
        // sample the point using direct importance sampling
        auto pt = sampler.Sample(emitter_sampling_tag::Direct{}, path.Back(), context, uniDist);

        return sample(
            pt
            , emitter_ = emitter.get()
        );
    }
};

// sample an emitter and sample a point on the emitter using area sampling
struct sample_position_on_emitter_t {
    template <typename T>
    auto operator()(Context<T>& context, const RandomSequence<Vertex>& path,
                    UniDist& uniDist, const ref_vector<Emitter> &emitters, const Float &time) const {
        // discrete_dynamic(emitters) creates a symbolic discrete random variable
        auto emitters_rv = discrete_dynamic(emitters);
        // context.Sample samples an emitter
        auto emitter = context.Sample(emitters_rv, context.Uniform1D(uniDist));
        // create a sampler for generating points on light source
        auto sampler = emitter->makeSampler();
        // sample the point using area sampling
        auto pt = sampler.Sample(emitter_sampling_tag::Position{}, context, uniDist, time);

        return sample(
            pt
            , emitter_ = emitter.get()
        );
    }
};

// given a point on the emitter, sample a direction
struct sample_direction_from_emitter_t {
    template <typename T>
    auto operator()(Context<T>& context, const RandomSequence<Vertex>& path,
                    UniDist& uniDist, const Raycast& raycaster) const {
        // make sure we start from a path with only the emitter node
        SAssert(path.Size() == 1);
        auto emitter = path.Back().Get(emitter_);
        auto sampler = emitter != nullptr ? emitter->makeSampler() : emitter_composite_t{};
        // sample a direction
        auto dirToNextWorld =
            sampler.Sample(emitter_sampling_tag::Direction{}, path.Back(), context, uniDist);

        auto currentPos = constant(path.Back().Value());
        const Intersection &currentIts = path.Back().Get(intersection_);
        // Use Mitsuba's raycasting engine to obtain intersection information
        const Ray ray(to_point(currentPos.Value()), to_vector(dirToNextWorld.Value()), currentIts.time);
        const Intersection nextIts = context.constant_call(raycaster, ray);
        // Do symbolic intersection
        auto intersectSample = intersect(currentPos, dirToNextWorld, nextIts);

        return intersectSample;
    }
};

void AppendPositionOnEmitter(RandomSequence<Vertex> &path, UniDist &uniDist,
                             const ref_vector<Emitter> &emitters, const Float time) {
    Node<sample_position_on_emitter_t> samplePositionOnEmitter;
    path.Append(samplePositionOnEmitter, uniDist, emitters, time);
}

void AppendDirectionFromEmitter(RandomSequence<Vertex> &path, UniDist &uniDist,
                                const Raycast &raycaster) {
    Node<sample_direction_from_emitter_t> sampleDirectionFromEmitter;
    path.Append(sampleDirectionFromEmitter, uniDist, raycaster);
}

void AppendDirectSampleEmitter(RandomSequence<Vertex> &path, UniDist &uniDist,
                               const ref_vector<Emitter> &emitters) {
    Node<direct_sample_emitter_t> directSampleEmitter;
    path.Append(directSampleEmitter, uniDist, emitters);
}

MTS_NAMESPACE_END
