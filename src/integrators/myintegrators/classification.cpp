#include "classification.h"

#include <mitsuba/render/scene.h>

MTS_NAMESPACE_BEGIN

bool isSpecular(const Vertex &vertex, const Float specularRoughnessThreshold) {
    if (!vertex.Valid()) {
        return false;
    }

    const auto &intersection = vertex.Get(intersection_);
    if (!intersection.isValid()) {
        return false;
    }

    if (intersection.shape == nullptr) {
        return false;
    }

    const BSDF *bsdf = intersection.getBSDF();
    if (bsdf == nullptr) {
        return false;
    }

    // TODO: make this probabilistic
    for (int i=0; i<bsdf->getComponentCount(); ++i) {
        if (bsdf->getType(i) & BSDF::ESmooth) {
            if (bsdf->getRoughness(intersection, i) > specularRoughnessThreshold) {
                return false;
            }
        }
    }

    // Phew!
    return true;
}

bool isDiffuse(const Vertex &vertex, const Float specularRoughnessThreshold) {
    if (!vertex.Valid()) {
        return false;
    }

    const auto &intersection = vertex.Get(intersection_);
    if (!intersection.isValid()) {
        return false;
    }

    if (intersection.shape == nullptr) {
        return false;
    }

    const BSDF *bsdf = intersection.getBSDF();
    if (bsdf == nullptr) {
        return false;
    }

    // TODO: make this probabilistic
    for (int i = 0; i < bsdf->getComponentCount(); ++i) {
        if (bsdf->getType(i) & BSDF::ESmooth) {
            if (bsdf->getRoughness(intersection, i) > specularRoughnessThreshold) {
                // Phew!
                return true;
            }
        }
    }

    return false;
}

MTS_NAMESPACE_END
