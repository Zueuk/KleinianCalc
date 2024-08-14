#pragma once

#include <memory>
#include <vector>
#include <complex>

#include "Kleinian.h"
#include "Moebius.hpp"

#include <QImage>

class Renderer
{
public:
	using real = float;
	using complex = std::complex<real>;

private:
	long width, height;

	std::unique_ptr<float[]> histogram;

	std::vector<Moebius<complex>> transforms;
	Moebius<complex> cameraTransform;

	long itersCounter;

public:
	Renderer(long w, long h)
		: width(w)
		, height(h)
		, histogram(new float[w * h])
	{ }

	void clear();
	void init(const Kleinian& K, const complex& root, const Moebius<complex>& camera);
	void iterate(long iters);
	void tonemap(QImage& image);
};
