#pragma once
#if !defined(__MITSUBA_SAMPLING_UTILS_H_)
#define __MITSUBA_SAMPLING_UTILS_H_

#include <mitsuba/core/spectrum.h>
#include <mitsuba/render/scene.h>
#include <mitsuba/render/vertex.h>
#include <aether/RandomSequence.h>

#include <vector>

MTS_NAMESPACE_BEGIN

class Sensor;

Spectrum estimate(const Scene *scene, const RandomSequenceView<Vertex> &path);
Spectrum estimateBidir(const Scene *scene, const RandomSequenceView<Vertex> &path);
Float misWeight(Float pdfA, Float pdfB);

template <typename T>
double misWeightImpl(const size_t index, const T &pdfs) {
	// power heuristics
	const double sq_numerator = pdfs[index];
	if (sq_numerator <= 0. || !std::isfinite(sq_numerator)) {
		// This shouldn't happen though
		return 0.;
	}
	double denominator = 0.0;
	for (size_t i = 0; i < pdfs.size(); i++) {
		if (pdfs[i] <= 0. || !std::isfinite(pdfs[i])) {
			continue;
		}
		double ratio = pdfs[i] / sq_numerator;
		denominator += (ratio * ratio);
	}
	SAssert(denominator >= 0.);

	return 1.0 / denominator;
}

template<int N>
double misWeight(const size_t index, const std::array<double, N> &pdfs) {
	return misWeightImpl(index, pdfs);
}

inline double misWeight(const size_t index, const std::vector<double> &pdfs) {
	return misWeightImpl(index, pdfs);
}

Spectrum estimate(const Scene *scene, const RandomSequenceView<Vertex> &path);
Spectrum estimateBidir(const Scene *scene, const RandomSequenceView<Vertex> &path);
Float misWeight(Float pdfA, Float pdfB);
double misWeight(const size_t index, const std::vector<double> &pdfs);
bool project(const Sensor *sensor, const Point &pCamera, const Point &pSurface, const Float time, Point2 &position);
bool project(const Scene *scene, const RandomSequence<Vertex> &path, Point2 &position);

struct PowerHeuristic {
  Float operator()(std::size_t index, const std::vector<double>& pdfs) const;
};

MTS_NAMESPACE_END

#endif // __MITSUBA_SAMPLING_UTILS_H_
