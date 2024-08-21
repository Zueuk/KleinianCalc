#pragma once

#include <boost/multiprecision/cpp_complex.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/tools/polynomial.hpp>

#include "Moebius.hpp"

inline bool isInf(int n) { return n == INT_MAX; }

template <class realT>
constexpr realT offsetN(int n) {
	if (n == 2)
		return 0;
	if (isInf(n))
		return 2;
	return 2 * cos(boost::math::constants::pi<realT>() / n);
}

class Kleinian
{
	// Using boost::multiprecision to solve higher order polynomials
	using complex = boost::multiprecision::cpp_complex_oct;
	using real = complex::value_type;

	// Using boost::polynomial class for convenience
	using polynomial = boost::math::tools::polynomial<complex>;

	std::vector<complex> findRoots_Newton(polynomial& P, complex guess, bool findAllRoots) const;
	std::vector<complex> findRoots_SkowronGould(polynomial& P, complex guess, bool findAllRoots) const;

public:
	int a, b;
	int n, h;
	int v1, v2;

	Kleinian() = default;

	template <typename moebiusT>
	moebiusT createSolverMatrix(const moebiusT& Ma, const moebiusT& Mb, const moebiusT& Mc) const;

	std::vector<complex> solve(int method, bool findAllRoots) const;

	template <typename complexT>
	std::vector<Moebius<complexT>> createGenerators(complexT root, bool asDisc) const;
};

template <typename moebiusT>
moebiusT Kleinian::createSolverMatrix(const moebiusT& Ma, const moebiusT& Mb, const moebiusT& Mc) const {
	const moebiusT& MB = (v1 != 2) ? Mb : Mc;
	const moebiusT& MC = (v2 != 2) ? Mc : Mb;

	auto stepA = [&](moebiusT& m) {
		m = m * Ma;
	};
	bool toggle = false;
	auto stepB = [&](moebiusT& m) {
		m = m * ((toggle = !toggle) ? MB : MC);
	};

	moebiusT M;

	// This is basically the line drawing algorithm
	if (a >= b) {
		int x = 0;
		for (int i = 0; i < a; ++i) {
			stepA(M);
			x += b;
			if (x >= a) {
				stepB(M);
				x -= a;
			}
		}
	}
	else {
		int y = 0;
		for (int i = 0; i < b; ++i) {
			stepB(M);
			y += a;
			if (y >= b) {
				stepA(M);
				y -= b;
			}
		}
	}

	return M;
}

template <typename realT>
class TilngParameters
{
	static const inline realT PI = boost::math::constants::pi<realT>();

public:
	realT ph, qh, rh;
	realT A, B, C, B1, C1;

	// Hyperbolic tiling with p pairs of q- and r-gons
	TilngParameters(int p, int q, int r) {
		if (isInf(p)) {
			A = 0;
			B = isInf(q) ? 0 : PI/q;
			C = isInf(r) ? 0 : PI/r;
			B1 = 0;
			C1 = 0;

			ph = 1;
			qh = sin(PI/4 - B/2) / sin(PI/4 + B/2);
			rh = sin(PI/4 - C/2) / sin(PI/4 + C/2);
		}
		else {
			A = PI/p;
			B = isInf(q) ? 0 : PI/q;
			C = isInf(r) ? 0 : PI/r;
			B1 = atan( sin(A) / (cos(C)/cos(B) + cos(A)) );
			C1 = atan( sin(A) / (cos(B)/cos(C) + cos(A)) );

			realT t = cos(B) / sin(B1);
			ph = sqrt( (t - 1)/(t + 1) );

			qh = ph * tan(PI/4 + B1/2 - B/2);
			rh = ph * tan(PI/4 + C1/2 - C/2);
		}
	}
};

template <typename complexT>
void cleanMatrix(Moebius<complexT>& M) {
	auto cleanValue = [](complexT& z) {
		if (abs(z.real()) < FLT_EPSILON) z.real(0);
		if (abs(z.imag()) < FLT_EPSILON) z.imag(0);
	};
	cleanValue(M.a);
	cleanValue(M.b);
	cleanValue(M.c);
	cleanValue(M.d);
};

template <typename complexT>
std::vector<Moebius<complexT>> Kleinian::createGenerators(complexT root, bool asDisc) const {
	Moebius<complexT> Ma, Mb, Mc;
	using realT = typename complexT::value_type;

	static const Moebius<complexT> Disc2halfplane = {
		1, 1,
		-1, 1
	};
	static const Moebius<complexT> Halfplane2disc = {
		1, -1,
		1, 1
	};

	if (isInf(h)) {
		Ma = {
			0, 1,
			-1.0, root
		};
		Mb = {
			0, 1,
			1, complexT(0, offsetN<realT>(v1))
		};
		Mc = {
			0, 1,
			1, complexT(0, -offsetN<realT>(v2))
		};

		if (asDisc) {
			Ma = Halfplane2disc * Ma * Disc2halfplane;
			Mb = Halfplane2disc * Mb * Disc2halfplane;
			Mc = Halfplane2disc * Mc * Disc2halfplane;

			Ma.divideBy(Ma.c);
			Mb.divideBy(Mb.c);
			Mc.divideBy(Mc.c);

			cleanMatrix(Ma);
			cleanMatrix(Mb);
			cleanMatrix(Mc);
		}
	}
	else {
		TilngParameters<realT> tiling(h, v1, v2);

		realT pZ = -tiling.ph / (1 + tiling.ph * tiling.ph);
		Ma = {
			pZ + pZ * root, root,
			1, pZ + pZ * root
		};

		auto NormalForm2fp = [](complexT fp1, complexT fp2, complexT k) {
			if (isinf(fp2.real()) || isinf(fp2.imag()))
				return Moebius<complexT> {
					k, (complexT(1) - k)*fp1,
					0, 1
				};
			else
				return Moebius<complexT> {
					(fp1 - k*fp2), (k - complexT(1))*fp1*fp2,
					(complexT(1) - k), (k*fp1 - fp2)
				};
		};
		auto NormalForm1fp = [](complexT fp, complexT b) {
			return Moebius<complexT> {
				b + fp, -fp*fp,
				1, b - fp
			};
		};

		if (!isInf(v1))
			Mb = NormalForm2fp(complexT(0, -tiling.qh), complexT(0, -1/tiling.qh), exp(complexT(0, 2*tiling.B)));
		else
			Mb = NormalForm1fp(complexT(0, -1), -0.5*tiling.ph - 0.5/tiling.ph);

		if (!isInf(v2))
			Mc = NormalForm2fp(complexT(0, tiling.rh), complexT(0, 1/tiling.rh), exp(complexT(0, -2*tiling.C)));
		else
			Mc = NormalForm1fp(complexT(0, 1), -0.5*tiling.ph - 0.5/tiling.ph);

		if (!asDisc) {
			Ma = Disc2halfplane * Ma * Halfplane2disc;
			Mb = Disc2halfplane * Mb * Halfplane2disc;
			Mc = Disc2halfplane * Mc * Halfplane2disc;

			Ma.divideBy(-Ma.c);
		}
		Mb.divideBy(Mb.c);
		Mc.divideBy(Mc.c);

		cleanMatrix(Ma);
		cleanMatrix(Mb);
		cleanMatrix(Mc);
	}

	return { Ma, Mb, Mc };
}
