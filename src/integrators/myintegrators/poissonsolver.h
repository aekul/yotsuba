#pragma once

#include <Eigen/Dense>
#include <numeric>
#include <vector>

namespace ps {

typedef Eigen::Matrix<double, 3, 1> Vector3;

struct Image3 {
	Image3() {
	}
	Image3(int w, int h) : pixelWidth(w), pixelHeight(h) {
		data.resize(w * h);
		std::fill(data.begin(), data.end(), Vector3::Zero());
	}
	Image3(int w, int h, const Vector3 &init) : pixelWidth(w), pixelHeight(h) {
		data.resize(w * h);
		std::fill(data.begin(), data.end(), init);
	}

	Vector3 &At(int x) {
		return data[x];
	}
	const Vector3 &At(int x) const {
		return data[x];
	}
	Vector3 &At(int x, int y) {
		return data[y * pixelWidth + x];
	}
	const Vector3 &At(int x, int y) const {
		return data[y * pixelWidth + x];
	}
	int NumPixels() const {
		return pixelWidth * pixelHeight;
	}
	void Clear() {
		std::fill(data.begin(), data.end(), Vector3::Zero());
	}
	Vector3 Sum() {
		Vector3 sum = Vector3::Zero();
		return std::accumulate(data.begin(), data.end(), sum);
	}

	int pixelWidth;
	int pixelHeight;
	std::vector< Vector3, Eigen::aligned_allocator<Vector3> > data;
};

void Solve(const Image3 &throughputImage,
		   const Image3 &dxImage,
		   const Image3 &dyImage,
		   const Image3 &visibleEmitterImage,
		   const double alpha,
		   const bool doL1,
		   Image3 &reconstructImage);

}
