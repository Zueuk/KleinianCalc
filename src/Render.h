#pragma once

#include <memory>
#include <vector>
#include <complex>

#include "Kleinian.h"

#include <QImage>

using real_t = float;
using complex_t = std::complex<real_t>;

class Renderer
{
	long width, height;

	std::unique_ptr<float[]> histogram;

	std::vector<Moebius<complex_t>> transforms;
	Moebius<complex_t> cameraTransform;

	long itersCounter;

public:
	Renderer(long w, long h)
		: width(w)
		, height(h)
		, histogram(new float[w * h])
	{ }

	void clear();
	void init(const Kleinian& K, const complex_t& root, const Moebius<complex_t>& camera);
	void iterate(long iters);
	void tonemap(QImage& image);
};
