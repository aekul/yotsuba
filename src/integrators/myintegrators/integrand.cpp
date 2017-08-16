#include "integrand.h"
#include <mitsuba/render/scene.h>
#include <mitsuba/render/emitter.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/sensor.h>

using namespace aether;

MTS_NAMESPACE_BEGIN

inline mitsuba::Vector directionToPrevious(const Vertex &previousVertex, const Vertex &vertex) {
    return to_vector((previousVertex.Value() - vertex.Value()).normalized());
}

inline mitsuba::Vector directionToNext(const Vertex &vertex, const Vertex &nextVertex) {
    return to_vector((nextVertex.Value() - vertex.Value()).normalized());
}

Spectrum Le(const Scene *scene, const Vertex &previousVertex, const Vertex &vertex) {
    const Emitter *emitter = vertex.Get(emitter_);
    const Intersection &its = vertex.Get(intersection_);
    const mitsuba::Vector directionToLight = -directionToPrevious(previousVertex, vertex);
    if (emitter != nullptr) {
        // Sampled emitter
        return emitter->eval(its, -directionToLight);
    }
    if (its.isValid()) {
        // if its is an emitter, emitter shouldn't be nullptr
        SAssert(!its.isEmitter());
        if (its.isEmitter()) {
            // Hit an area light
            return its.Le(-directionToLight);
        }
        // Hit a non-emitter object
        return Spectrum(0.f);
    }
    // Hit background
    // It might seem not making sense using its.p, but evalEnvironment doesn't use ray.o anyway
    return scene->evalEnvironment(RayDifferential(its.p, directionToLight, its.time));
}

Spectrum eval(const Scene *scene, const RandomSequenceView<Vertex> &path, const int index) {
    SAssert(path.Size() > 1);
    SAssert(index >= 1);
    if (index == path.Size() - 1) {
        return Le(scene, path[index - 1], path[index]);
    }

    SAssert(index + 1 < path.Size());
    Intersection its = path[index].Get(intersection_);
    if (!its.isValid()) {
        // Not a valid surface
        return Spectrum(0.f);
    }
    SAssert(its.shape != nullptr);
    its.wi = its.toLocal(directionToPrevious(path[index - 1], path[index]));
    const mitsuba::Vector toNext = directionToNext(path[index], path[index + 1]);
    BSDFSamplingRecord bRec(its, its.toLocal(toNext));
    bRec.framePerturbed = true;
    // TODO: support ray differentials
    const BSDF *bsdf = its.getBSDF();
    SAssert(bsdf != nullptr);
    const Spectrum bsdfValue = bsdf->eval(bRec);
    if (bsdfValue.isZero()) {
        return Spectrum(0.f);
    }
    const Intersection &nextIts = path[index + 1].Get(intersection_);
    const Float nextCosine = nextIts.isValid() ? 
        fabs(dot(-toNext, nextIts.geoFrame.n)) : Float(1.f);
    const Float distSq = 
        distanceSquared(to_point(path[index].Value()),
                        to_point(path[index + 1].Value()));
    if (distSq == Float(0)) {
        return Spectrum(0.f);
    }
    // bsdf->eval already includes one of the cosines
    const Float geometryTerm = nextCosine / distSq;
    if (!std::isfinite(geometryTerm) || !bsdfValue.isValid()) {
        return Spectrum(0.f);
    }

    return geometryTerm * bsdfValue * eval(scene, path, index + 1);
}

TSpectrum<double, SPECTRUM_SAMPLES> eval(const Scene *scene, const RandomSequenceView<Vertex> &path) {
    SAssert(path.AllValid());
    if (path.Size() <= 1) {
        // Need to at least have a sensor and an emitter
        return TSpectrum<double, SPECTRUM_SAMPLES>(0.f);
    }
    return TSpectrum<double, SPECTRUM_SAMPLES>(eval(scene, path, 1));
}

Spectrum evalSensor(const Scene *scene, const RandomSequenceView<Vertex> &path) {
    // Should we obtain sensor from the path?
    const Sensor *sensor = scene->getSensor();
    SAssert(path[0].Get(sensor_) != nullptr);
    PositionSamplingRecord pRec(path[1].Get(intersection_).time);
    pRec.p = to_point(path[0].Value());
    DirectionSamplingRecord dRec;
    dRec.d = to_vector((path[1].Value() - path[0].Value()).normalized());
    dRec.measure = ESolidAngle;
    Spectrum We = sensor->evalDirection(dRec, pRec);
    const Intersection &nextIts = path[1].Get(intersection_);
    const Float nextCosine = nextIts.isValid() ? 
        fabs(dot(-dRec.d, nextIts.geoFrame.n)) : Float(1.f);
    // bsdf->eval already includes one of the cosines
    const Float geometryTerm = nextCosine /
        distanceSquared(to_point(path[0].Value()),
                        to_point(path[1].Value()));
    return geometryTerm * We;
}

TSpectrum<double, SPECTRUM_SAMPLES> evalBidir(const Scene *scene, const RandomSequenceView<Vertex> &path) {
    SAssert(path.AllValid());
    if (path.Size() <= 1) {
        // Need to at least have a sensor and an emitter
        return TSpectrum<double, SPECTRUM_SAMPLES>(0.f);
    }

    return TSpectrum<double, SPECTRUM_SAMPLES>(evalSensor(scene, path) * eval(scene, path, 1));
}


bool occluded(const Scene *scene, const Point &a, const Point &b, const Float time) {
    Vector dir = b - a;
    Float distance = dir.length();
    if (distance == Float(0)) {
        return false;
    }
    dir /= distance;
    Ray ray(a, dir, Epsilon, distance * (1 - ShadowEpsilon), time);
    return scene->rayIntersect(ray);
}

MTS_NAMESPACE_END
