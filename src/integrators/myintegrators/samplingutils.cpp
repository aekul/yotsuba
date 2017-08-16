#include "samplingutils.h"
#include "integrand.h"

MTS_NAMESPACE_BEGIN

Spectrum estimate(const Scene *scene, const RandomSequenceView<Vertex> &path) {
    SAssert(path.AllValid());
    // integrand
    TSpectrum<double, SPECTRUM_SAMPLES> Li = eval(scene, path);
    if (!Li.isValid() || Li.isZero()) {
        return Spectrum(0.f);
    }
    double pdf = double(path.Pdf(path));
    if (pdf == 0. || !std::isfinite(pdf)) {
        // TODO: issue warnings here
        return Spectrum(0.f);
    }

    return Spectrum(Li / pdf);
}

Spectrum estimateBidir(const Scene *scene, const RandomSequenceView<Vertex> &path) {
    SAssert(path.AllValid());
    // integrand
    TSpectrum<double, SPECTRUM_SAMPLES> Li = evalBidir(scene, path);
    if (!Li.isValid() || Li.isZero()) {
        return Spectrum(0.f);
    }
    double pdf = double(path.Pdf(path));
    if (pdf == 0. || !std::isfinite(pdf)) {
        // TODO: issue warnings here
        return Spectrum(0.f);
    }

    return Spectrum(Li / pdf);
}

Float misWeight(Float pdfA, Float pdfB) {
    if (pdfA <= Float(0) || !std::isfinite(pdfA)) {
        return Float(0);
    }
    if (pdfB <= Float(0) || !std::isfinite(pdfB)) {
        return Float(1);
    }
    // use double & ratio for numerical stability
    const double ratio = pdfB / pdfA;

    return 1.0 / (1.0 + ratio * ratio);
}

bool project(const Sensor *sensor, const Point &pCamera, const Point &pSurface, const Float time, Point2 &position) {
    PositionSamplingRecord pRec(time);
    pRec.p = pCamera;
    DirectionSamplingRecord dRec;
    dRec.d = normalize(pSurface - pCamera);
    return sensor->getSamplePosition(pRec, dRec, position);
}

bool project(const Scene *scene, const RandomSequence<Vertex> &path, Point2 &position) {
    SAssert(path.Size() >= 2);
    // Should we obtain sensor from path?
    const Sensor *sensor = scene->getSensor();
    return project(sensor, to_point(path[0].Value()), to_point(path[1].Value()), path[1].Get(intersection_).time, position);
}

Float PowerHeuristic::operator()(std::size_t index, const std::vector<double>& pdfs) const {
  return misWeight(index, pdfs);
}

MTS_NAMESPACE_END
