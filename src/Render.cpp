#include "Render.h"

#include <chrono>

#include <boost/math/constants/constants.hpp>

#include "pcg32.h"

static const auto PI = boost::math::constants::pi<Renderer::real>();

static const long batch = 1000;

void Renderer::clear() {
	for (long i = 0; i < width * height; ++i)
		histogram[i] = 0;
	itersCounter = 0;
}

void Renderer:: init(const Kleinian& K, const complex& root, const Moebius<complex>& camera) {
	transforms.clear();

	Moebius<complex> Ma = {
		0.0, 1.0,
		-1.0, root
	};
	Moebius<complex> Mb = {
		0.0, 1.0,
		1.0, complex(0, offsetN<real>(K.v1))
	};
	Moebius<complex> Mc = {
		0.0, 1.0,
		1.0, complex(0, -offsetN<real>(K.v2))
	};

	transforms.emplace_back(Ma);
	transforms.emplace_back(Ma.inverse());
	transforms.emplace_back(Mb);
	transforms.emplace_back(Mb.inverse());
	transforms.emplace_back(Mc);
	transforms.emplace_back(Mc.inverse());

	cameraTransform = camera;
}

void Renderer::iterate(long iters) {
	pcg32 rng(std::chrono::system_clock::now().time_since_epoch().count());
	const auto bound = (uint32_t)transforms.size();

	real camScale = std::min(width, height) / 2.f;
	real camX = width / 2.f;
	real camY = height / 2.f;

	for (long i = 0; i < iters; i += batch) {
		// point on circle -> vertical line
		real a = rng.nextFloat() * 2*PI;
		complex z(0, std::tan(a));

		for (long j = 0; j < batch; ++j) {
			auto n = rng.nextUInt(bound);
			z = transforms[n].apply(z);

			complex zc = cameraTransform.apply(z);

			auto x = long(zc.real() * camScale + camX);
			auto y = long(zc.imag() * camScale + camY);

			if (x >= 0 && x < width && y >= 0 && y < height) {
				histogram[y * width + x] += 1;
			}
		}
		itersCounter += batch;
	}
}

void Renderer::tonemap(QImage& image) {
	float k1 = 256 * log1p(1.f) / log1p(1 / 4.f);
	float k2 = float(width * height) / (4.f * itersCounter);

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
