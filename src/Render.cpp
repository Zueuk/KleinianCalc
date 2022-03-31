#include "Render.h"

#include <chrono>

#include <boost/math/constants/constants.hpp>

#include "pcg32.h"

static const real_t PI = boost::math::constants::pi<real_t>();

static const long batch = 100;

void Renderer::clear() {
	for (long i = 0; i < width * height; ++i)
		histogram[i] = 0;
	itersCounter = 0;
}

void Renderer:: init(const Kleinian& K, const complex_t& root, const Moebius<complex_t>& camera) {
	transforms.clear();

	Moebius<complex_t> Mv1{
		0.0, 1.0,
		1.0, complex_t(0, offsetN<real_t>(K.nV1))
	};
	transforms.emplace_back(Mv1);
	transforms.emplace_back(Mv1.inverse());

	Moebius<complex_t> Mv2{
		0.0, 1.0,
		1.0, complex_t(0, -offsetN<real_t>(K.nV2))
	};
	transforms.emplace_back(Mv2);
	transforms.emplace_back(Mv2.inverse());

	Moebius<complex_t> Ma = {
		0.0, 1.0,
		-1.0, root
	};
	transforms.emplace_back(Ma);
	transforms.emplace_back(Ma.inverse());

	cameraTransform = camera;
}

void Renderer::iterate(long iters) {
	pcg32 rng(std::chrono::system_clock::now().time_since_epoch().count());

	real_t camScale = std::min(width, height) / (real_t)2;
	real_t camX = width/2;
	real_t camY = height/2;

	for (long i = 0; i < iters; i += batch) {
		// point on circle -> vertical line
		real_t a = rng.nextFloat() * 2*PI;
		complex_t z(0, std::tan(a));

		for (long j = 0; j < batch; ++j) {
			int n = rng.nextUInt() % transforms.size();
			const auto& M = transforms[n];
			z = M.apply(z);

			complex_t zc = cameraTransform.apply(z);

			int x = int(zc.real() * camScale + camX);
			int y = int(zc.imag() * camScale + camY);

			if (x >= 0 && x < width && y >= 0 && y < height) {
				histogram[y * width + x] += real_t(1);
			}
		}
		itersCounter += batch;
	}
}

void Renderer::tonemap(QImage& image) {
	real_t k1 = 256 * log1p(1) / log1p(1 / 4.0);
	real_t k2 = real_t(width * height) / (4.0 * itersCounter);

	for (long y = 0; y < height; ++y) {
		auto line = image.scanLine(y);
		for (long x = 0; x < width; ++x) {
			float fc = histogram[y * width + x];
			float fl = k1 * log1pf(fc * k2);
			int c = std::clamp((int)fl, 0, 255);
			line[x] = c;
		}
	}
}
