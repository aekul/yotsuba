#ifndef OCCLUSION_CACHE_H__
#define OCCLUSION_CACHE_H__

#include <unordered_map>
#include "integrand.h"

MTS_NAMESPACE_BEGIN

using OcclusionKey = std::pair<aether::Vector3, aether::Vector3>;

inline std::size_t hash(const aether::Vector3 &v) {
    std::size_t seed = std::hash<aether::Real>{}(v[0]);
    boost::hash_combine(seed, v[1]);
    boost::hash_combine(seed, v[2]);
    return seed; 
}

struct OcclusionKeyHash {
    inline std::size_t operator()(const OcclusionKey& key) const {
        std::size_t seed = hash(key.first);
        boost::hash_combine(seed, hash(key.second));
        return seed;
    }
};

struct OcclusionCache {
    OcclusionCache(const Scene *scene, const Float time) : scene(scene), time(time) {}

    inline bool Query(const OcclusionKey &key) {
        auto search = cacheMap.find(key);
        if (search == cacheMap.end()) {
            bool ret = occluded(scene, to_point(key.first), to_point(key.second), time);
            cacheMap.insert({key, ret});
            return ret;
        }
        return search->second;
    }

    const Scene *scene;
    const Float time;
    std::unordered_map<OcclusionKey, bool, OcclusionKeyHash> cacheMap;
};

bool evalOcclusion(const RandomSequenceView<Vertex> &pathView, OcclusionCache &occlusionCache);

MTS_NAMESPACE_END

#endif // OCCLUSION_CACHE_H__
