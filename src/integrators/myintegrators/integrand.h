#if !defined(__INTEGRAND_H)
#define __INTEGRAND_H

#include <mitsuba/core/spectrum.h>
#include <mitsuba/render/vertex.h>
#include <aether/RandomSequence.h>

MTS_NAMESPACE_BEGIN

TSpectrum<double, SPECTRUM_SAMPLES> eval(const Scene *scene, const aether::RandomSequenceView<Vertex> &path);
TSpectrum<double, SPECTRUM_SAMPLES> evalBidir(const Scene *scene,
		const aether::RandomSequenceView<Vertex> &path);
bool occluded(const Scene *scene, const Point &a, const Point &b, const Float time);

MTS_NAMESPACE_END

#endif // __INTEGRAND_H
