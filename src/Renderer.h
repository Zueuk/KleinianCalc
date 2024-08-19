#pragma once

#include <complex>
#include <memory>
#include <vector>

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
	long itersCounter;

public:
	Renderer(long w, long h)
		: width(w)
		, height(h)
		, histogram(new float[w * h])
	{ }

	enum class ViewMode { Disc, Halfplane };

	void clear();
	void iterate(long iters, const std::vector<Moebius<complex>>& transforms, ViewMode calcMode, ViewMode viewMode);
	void tonemap(QImage& image);
};

extern const Moebius<Renderer::complex> viewMatrix[];
