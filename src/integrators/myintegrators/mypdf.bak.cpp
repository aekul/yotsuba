
double myPdf(const RandomSequence<Vertex> &target, const RandomSequence<Vertex> &given,
		Node<bidir_mutation_t> &strategy,
		int delete_length,
		int delete_begin,
		int insert_length) {
	int delete_end = delete_begin + insert_length;
	for (int i = 0; i <= delete_begin; i++) {
		if (target[i].Value() != given[i].Value()) {
			//std::cerr << "myPdf:" << 0 << std::endl;
			return 0.;
		}
	}
	for (int i = delete_begin + delete_length; i < given.Size(); i++) {
		int offset = given.Size() - 1 - i;
		int targetIndex = target.Size() - 1 - offset;

		if (target[targetIndex].Value() != given[i].Value()) {
			//std::cerr << "myPdf:" << 0 << std::endl;
			return 0.;
		}
	}
	//int delete_begin = strategy.sampled_delete_begin;
	//int delete_length = strategy.forward_phase ? strategy.sampled_delete_length : strategy.sampled_insert_length;
	//int insert_length = strategy.forward_phase ? strategy.sampled_insert_length : strategy.sampled_delete_length;
	double desired_length_pmf = 0.;
	double desired_length_sum = 0.;
	const int path_length = given.Size() + 1;
	//std::cerr << "path_length:" << path_length << std::endl;
	BidirectionalMutation *mutation = strategy.mutation;
	for (int length = mutation->min_length; length <= mutation->max_length; length++) {
		const int offset = std::abs(length - path_length);
		const double w = pow(2.0, -offset);
		desired_length_sum += w;
		if (length == target.Size() + 1) {
			desired_length_pmf = w;
		}
	}
	desired_length_pmf /= desired_length_sum;
	const int target_length = target.Size() + 1;
	//std::cerr << "target_length:" << target_length << std::endl;
	//std::cerr << "target_length_pdf:" << desired_length_pmf << std::endl;
	double delete_length_pmf = 0.;
	double delete_length_sum = 0.;
	int min_deletion = std::max((path_length == target_length) ? 2 : 1, path_length-target_length+1);
	for (int length = min_deletion; length <= path_length; length++) {
		const int offset = std::abs(length - 2);
		const double w = pow(2.0, -offset);
		delete_length_sum += w;
		if (length == delete_length) {
			delete_length_pmf = w;
		}
	}
	delete_length_pmf /= delete_length_sum;

	//std::cerr << "delete_length:" << delete_length << std::endl;
	//std::cerr << "delete_length_pdf:" << delete_length_pmf << std::endl;

	//std::cerr << "delete_begin:" << delete_begin << std::endl;
	const Real specular_roughness_threshold = mutation->specular_roughness_threshold;
	int min_index = 0;
	int max_index = path_length - delete_length - 1;
	if (delete_length == 1 || insert_length == 1) {
		max_index--;
	}
	int delete_begin_count = 0;
	for (int index = min_index; index <= max_index; index++) {
		bool left_dirac = index > 0 && index < given.Size() - 1 && is_specular(given[index], specular_roughness_threshold);
		const int right_index = index + delete_length;
		bool right_dirac = right_index > 0 && right_index < given.Size() - 1 && is_specular(given[right_index], specular_roughness_threshold);
		if (!left_dirac && !right_dirac) {
			delete_begin_count++;
		}
	}
	//std::cerr << "delete_begin_count:" << delete_begin_count << std::endl;
	double delete_begin_pmf = 1.0 / delete_begin_count;

	SAssert(path_length - delete_length + insert_length == target_length);
	//std::cerr << "insert_length:" << insert_length << std::endl;
	const int insert_size = insert_length - 1;
	//std::cerr << "insert_size:" << insert_size << std::endl;
	const Scene *scene = mutation->scene;
	double pdfSum = 0.;
	//SourceDistributionHelper source_distribution_helper;
	for (int eye_insert = 0; eye_insert <= insert_size; eye_insert++) {
		//std::cerr << "eye_insert:" << eye_insert << std::endl;
		double pdf = 1.;
		for (int index = delete_begin; index < delete_begin + eye_insert; index++) {
			double thisPdf = 0.;
			//std::cerr << "[eye] index:" << index << std::endl;
			if (index == 0) {
				DirectionSamplingRecord dRec;
				dRec.measure = ESolidAngle;
				PositionSamplingRecord pRec;
				pRec.time = 0.5;
				dRec.d = to_vector((target[index + 1].Value() - target[index].Value()).normalized());
				const Sensor *sensor = scene->getSensor();
				const double sensorPdf = sensor->pdfDirection(dRec, pRec);
				//std::cerr << "sensorPdf:" << sensorPdf << std::endl;
				thisPdf = sensorPdf;
			} else if (index > 0) {
				Intersection its = target[index].Get(intersection_);
				its.wi = its.toLocal(to_vector((target[index - 1].Value() - target[index].Value()).normalized()));
				const Vector3 wo = its.toLocal(to_vector((target[index + 1].Value() - target[index].Value()).normalized()));
				const BSDF *bsdf = its.getBSDF();
				BSDFSamplingRecord bRec(its, wo, ERadiance);
				const double bsdfPdf = bsdf->pdf(bRec);
				//std::cerr << "bsdfPdf:" << bsdfPdf << std::endl;
				thisPdf = bsdfPdf;
			}
			if (index >= 0) {
				const Vector d = to_vector((target[index + 1].Value() - target[index].Value()).normalized());
				const double length = (target[index + 1].Value() - target[index].Value()).norm();
				thisPdf /= (length * length);
				const Intersection &its = target[index + 1].Get(intersection_);
				thisPdf *= fabs(dot(its.geoFrame.n, d));
			}

			//std::cerr << "thisPdf:" << thisPdf << std::endl;
			//std::cerr << "mcPdf:" << mcPdf << std::endl;
			pdf *= thisPdf;
		}

		for (int index = delete_end; index > delete_begin + eye_insert + 1; index--) {
			double thisPdf = 0.;
			//std::cerr << "[light] index:" << index << std::endl;
			if (index == target.Size()) {
				const Emitter *emitter = target[index - 1].Get(emitter_);
				PositionSamplingRecord pRec;
				const double pdfEmitterDiscrete = scene->pdfEmitterDiscrete(emitter);
				const double pdfPosition = emitter->pdfPosition(pRec);
				//std::cerr << "pdfEmitterDiscrete:" << pdfEmitterDiscrete << std::endl;
				//std::cerr << "pdfPosition:" << pdfPosition << std::endl;
				thisPdf = pdfEmitterDiscrete * pdfPosition;
			} else if (index == target.Size() - 1) {
				const Emitter *emitter = target[index].Get(emitter_);
				DirectionSamplingRecord dRec;
				dRec.measure = ESolidAngle;
				PositionSamplingRecord pRec;
				pRec.n = target[index].Get(intersection_).shFrame.n;
				dRec.d = to_vector((target[index - 1].Value() - target[index].Value()).normalized());
				const double pdfDirection = emitter->pdfDirection(dRec, pRec);
				//std::cerr << "pdfDirection:" << pdfDirection << std::endl;
				thisPdf = pdfDirection;
			} else {
				Intersection its = target[index].Get(intersection_);
				its.wi = its.toLocal(to_vector((target[index - 1].Value() - target[index].Value()).normalized()));
				const Vector3 wo = its.toLocal(to_vector((target[index + 1].Value() - target[index].Value()).normalized()));
				const BSDF *bsdf = its.getBSDF();
				BSDFSamplingRecord bRec(its, wo, EImportance);
				const double bsdfPdf = bsdf->pdf(bRec);
				//std::cerr << "bsdfPdf:" << bsdfPdf << std::endl;
				thisPdf = bsdfPdf;
			}
			if (index < target.Size()) {
				const Vector d = to_vector((target[index - 1].Value() - target[index].Value()).normalized());
				const double length = (target[index - 1].Value() - target[index].Value()).norm();
				thisPdf /= (length * length);
				const Intersection &its = target[index - 1].Get(intersection_);
				thisPdf *= fabs(dot(its.geoFrame.n, d));
			}
			//const int storeIndex = index - 1;
			//double mcPdf = target.stores[storeIndex].strategy->Pdf(
				//target, target[storeIndex], storeIndex, target.stores[storeIndex].reversed, source_distribution_helper);
			//std::cerr << "thisPdf:" << thisPdf << std::endl;
			//std::cerr << "mcPdf:" << mcPdf << std::endl;
			pdf *= thisPdf;
		}

		//std::cerr << "myPdf:" << pdf << std::endl;
		pdfSum += pdf * (1.0 / (insert_size + 1));
	}
	if (insert_size == 0) {
		pdfSum = 1.;
	}
	//std::cerr << "pdfSum:" << pdfSum << std::endl;

	return pdfSum * delete_length_pmf * desired_length_pmf * delete_begin_pmf;
}

double myPdf(const RandomSequence<Vertex> &target, const RandomSequence<Vertex> &given,
		Node<bidir_mutation_t> &strategy) {
	const int path_length = given.Size() + 1;
	const int target_length = target.Size() + 1;
	double pdfSum = 0.;
	int min_deletion = std::max((path_length == target_length) ? 2 : 1, path_length-target_length+1);
	for (int delete_length = min_deletion; delete_length <= path_length; delete_length++) {
		const int insert_length = target_length - path_length + delete_length;
		BidirectionalMutation *mutation = strategy.mutation;
		const Real specular_roughness_threshold = mutation->specular_roughness_threshold;
		int min_index = 0;
		int max_index = path_length - delete_length - 1;
		if (delete_length == 1 || insert_length == 1) {
			max_index--;
		}
		for (int delete_begin = min_index; delete_begin <= max_index; delete_begin++) {
			bool left_dirac = delete_begin > 0 && delete_begin < given.Size() - 1 && is_specular(given[delete_begin], specular_roughness_threshold);
			const int right_index = delete_begin + delete_length;
			bool right_dirac = right_index > 0 && right_index < given.Size() - 1 && is_specular(given[right_index], specular_roughness_threshold);
			if (!left_dirac && !right_dirac) {
				//std::cerr << "delete_begin:" << delete_begin << std::endl;
				//std::cerr << "delete_length:" << delete_length << std::endl;
				//std::cerr << "insert_length:" << insert_length << std::endl;
				double subPdf = myPdf(target, given, strategy, delete_length, delete_begin, insert_length);
				//std::cerr << "subPdf:" << subPdf << std::endl;
				pdfSum += subPdf;
			}
		}
	}
	return pdfSum;
}