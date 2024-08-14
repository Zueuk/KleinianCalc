#include "Kleinian.h"

#include "Moebius.hpp"
#include "SkowronGould.hpp"

#include <boost/math/tools/roots.hpp>

std::vector<Kleinian::complex> Kleinian::findRoots_Newton(Kleinian::polynomial& P, complex guess, bool findAllRoots) const {
	std::vector<complex> roots;

	auto minDegree = findAllRoots ? 0 : P.degree() - 1;
	while (P.degree() > minDegree) {
		// Trivial case
		if (P.degree() == 1) {
			roots.push_back(-P[0] / P[1]);
			break;
		}

		auto functorTuple = [&](complex z) -> auto {
			// Calculate both P(z) and P'(z) at once
			complex p = P[P.degree()];
			complex dp = 0.0;
			for (long k = (long)P.degree() - 1; k >= 0; --k) {
				dp = p + dp * z;
				p = P[k] + p * z;
			}
			return std::pair<complex, complex>(p, dp);
		};

		// Find one root using Newton's method
		complex root = boost::math::tools::complex_newton(functorTuple, guess);

		roots.push_back(root);

		// Remove this root from the polynomial
		complex coef = P[P.degree()];
		for (long i = (long)P.degree() - 1; i >= 0; --i) {
			complex prev = P[i];
			P[i] = coef;
			coef = prev + root * coef;
		}
		P.data().pop_back();
	}

	return roots;
}

std::vector<Kleinian::complex> Kleinian::findRoots_SkowronGould(Kleinian::polynomial& P, complex guess, bool findAllRoots) const {
	if (!findAllRoots) {
		std::vector<complex> roots(1, guess);
		SkowronGould::solve_one(P, roots[0], true);
		return roots;
	}
	else {
		std::vector<complex> roots(P.degree(), guess);
		SkowronGould::solve(P, roots, true);
		return roots;
	}
}

std::vector<Kleinian::complex> Kleinian::solve(int method, bool findAllRoots) const {
	complex cv1, cv2;
	if (v1 == 2) {
		cv2 = complex(0, offsetN<real>(v2));
		cv1 = cv2;
	}
	else if (v2 == 2) {
		cv1 = complex(0, offsetN<real>(v1));
		cv2 = cv1;
	}
	else {
		cv1 = complex(0, offsetN<real>(v1));
		cv2 = complex(0, offsetN<real>(v2));
	}

	using PolyMoebius = Moebius<polynomial>;

	// Create matrices out of polynomials
	PolyMoebius Ma = {
		{ 0.0 }, { 1.0 },
		{ -1.0 }, { 0.0, 1.0 }
	};
	// Matrices used for solving differ from the ones used for generation
	PolyMoebius Mb = Ma * PolyMoebius{
		{ 1.0 }, { cv1 },
		{ 0.0 }, { 1.0 }
	};
	PolyMoebius Mc = Ma * PolyMoebius{
		{ 1.0 }, { cv2 },
		{ 0.0 }, { 1.0 }
	};

	PolyMoebius M;

	auto applyA = [&](PolyMoebius& m) {
		m = m * Ma;
	};
	auto applyB = [&](PolyMoebius& m) {
		static bool toggle = false;
		m = m * ((toggle = !toggle) ? Mb : Mc);
	};

	// This is basically the line drawing algorithm
	if (a >= b) {
		int x = 0;
		for (int i = 0; i < a; ++i) {
			applyA(M);
			x += b;
			if (x >= a) {
				applyB(M);
				x -= a;
			}
		}
	}
	else {
		int y = 0;
		for (int i = 0; i < b; ++i) {
			applyB(M);
			y += a;
			if (y >= b) {
				applyA(M);
				y -= b;
			}
		}
	}

	// We know that the generated matrix M is an elliptic Moebius transformation
	// with period = n, so its trace (M.a + M.d) must be equal to 2cos(pi/n)
	//
	// M.a + M.d - 2cos(pi/n) = 0
	//
	// Since members of M are polynomials, this equation
	// becomes a polynomial that we're going to solve

	polynomial P = M.a + M.d - offsetN<real>(n);

	// We are looking for the root closest to (2 + 0i),
	// a small imaginary value is added to improve root finding
	if (method == 0) {
		return findRoots_Newton(P, complex(2.0, 0.1), findAllRoots);
	}
	else {
		return findRoots_SkowronGould(P, complex(2.0, 0.1), findAllRoots);
	}
}
