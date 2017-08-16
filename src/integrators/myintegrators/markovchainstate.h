#pragma once

#include <mitsuba/core/spectrum.h>
#include <mitsuba/core/vector.h>
#include <aether/RandomSequence.h>
#include "samplingutils.h"
#include "integrand.h"
#include "occlusioncache.h"

MTS_NAMESPACE_BEGIN

struct MarkovChainState {
	TSpectrum<double, SPECTRUM_SAMPLES> contribution;
	double importance;
	Point2 screenPosition;
	RandomSequence<Vertex> path;

	bool Valid() const {
		return path.store.size() > 0 && path.AllValid();
	}
};

inline double importance(const Spectrum &contribution) {
	return contribution.getLuminance();
}

inline double importance(const TSpectrum<double, SPECTRUM_SAMPLES> &contribution) {
	return importance(Spectrum(contribution));
}

inline MarkovChainState f(const Scene *scene, const RandomSequence<Vertex> &path,
						  OcclusionCache &occlusionCache) {
	RandomSequenceView<Vertex> pathView{path};
	if (pathView.Size() == 0 || !pathView.AllValid()) {
		return MarkovChainState{};
	}
	Point2 screenPosition;
	if (!project(scene, path, screenPosition)) {
		return MarkovChainState{};
	}
	const TSpectrum<double, SPECTRUM_SAMPLES> value = evalBidir(scene, pathView);
	if (value.isZero()) {
		return MarkovChainState{};
	}
	if (!value.isValid()) {
		return MarkovChainState{};
	}
	if (evalOcclusion(pathView, occlusionCache)) {
		return MarkovChainState{};
	}
	const double imp = importance(value);
	if (imp <= 0.) {
		return MarkovChainState{};
	}
	return MarkovChainState{value / imp, imp, screenPosition, path};
}

MTS_NAMESPACE_END
