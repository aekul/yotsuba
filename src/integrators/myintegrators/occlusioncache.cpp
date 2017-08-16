#include "occlusioncache.h"

MTS_NAMESPACE_BEGIN

bool evalOcclusion(const RandomSequenceView<Vertex> &pathView, OcclusionCache &occlusionCache) {
    for (size_t i = 0; i < pathView.Size() - 1; i++) {
        if (occlusionCache.Query(std::make_pair(pathView[i].Value(), pathView[i + 1].Value()))) {
            return true;
        }
    }
    return false;
}

MTS_NAMESPACE_END
