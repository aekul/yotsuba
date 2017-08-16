#pragma once

MTS_NAMESPACE_BEGIN

struct perturb_lens {
	const Scene *scene;
	const Real r1Value;
	const Real coveredAreaValue;
	UniDist &uniDist;
	const Raycast &raycaster;

	template <typename T>
	auto operator()(Context<T>& context, const RandomSequence<Vertex> &path) const {
		const Sensor *sensor = scene->getSensor();
		const Vector2 resolution = sensor->getResolution();
		const Vector2 invResolution = sensor->getInvResolution();
		// Generate an offset on screen
		auto r1 = constant(r1Value);
		auto r2 = constant(
				Real(std::sqrt(coveredAreaValue * Real(resolution[0]) * Real(resolution[1]) / M_PI)));
		auto r = r2 * exp(-log(r2 / r1) * u1);
		auto phi = two * get_pi() * u2;
		auto offset = make_random_vector(r * cos(phi) * constant(Real(invResolution.x)),
										 r * sin(phi) * constant(Real(invResolution.y)));

		// Get previous screen position
		Point2 previousScreenPosition;
		bool insideScreen = project(sensor, 
			to_point(path[0].Value()), to_point(path[1].Value()),
			path[1].Get(intersection_).time, previousScreenPosition);
		previousScreenPosition[0] *= Real(invResolution.x);
		previousScreenPosition[1] *= Real(invResolution.y);

		// Apply the screen offset
		auto uv = context.Uniform2D(uniDist);
		auto newScreenPosition = (constant(to_vector2(previousScreenPosition)) + offset).Sample(uv[0], uv[1]);
		insideScreen = insideScreen && (
			at<0>(newScreenPosition).Value() >= 0 && at<0>(newScreenPosition).Value() < 1 &&
			at<1>(newScreenPosition).Value() >= 0 && at<1>(newScreenPosition).Value() < 1);

		// Generate position and direction
		auto sensorSampler = sensor->makeSampler();
		auto p = constant(path[0].Value());
		auto dir = sensorSampler.Sample(sensor_sampling_tag::Direction{},
			path[0], context, uniDist, newScreenPosition);

		// Use Mitsuba's raycasting engine to obtain intersection information
		const Ray ray(to_point(p.Value()), to_vector(dir.Value()), path[1].Get(intersection_).time);
		const Intersection its = context.constant_call(raycaster, ray);
		// Do symbolic intersection
		auto intersectSample = intersect(p, dir, its);

		// TODO: setup its.emitter
		const Emitter *emitter = its.shape != nullptr ? its.shape->getEmitter() : nullptr;

		return optional_sample(
			insideScreen
			, intersectSample
			, emitter_ = emitter
		);
	}
};

MTS_NAMESPACE_END
