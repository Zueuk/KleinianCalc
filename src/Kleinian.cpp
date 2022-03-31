#include "Kleinian.h"

#include "SkowronGould.hpp"

#include <boost/math/tools/roots.hpp>

std::vector<Kleinian::complex> Kleinian::findRoots_Newton(Kleinian::polynomial& P, complex guess, bool findAllRoots) {
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

std::vector<Kleinian::complex> Kleinian::findRoots_SkowronGould(Kleinian::polynomial& P, complex guess, bool findAllRoots) {
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

std::vector<Kleinian::complex> Kleinian::solve(int method, bool findAllRoots) {
	// Calculate stripe offsets
	complex v1, v2;
	if (nV1 == 2) {
		v2 = complex(0, offsetN<real_t>(nV2));
		v1 = v2;
	}
	else if (nV2 == 2) {
		v1 = complex(0, offsetN<real_t>(nV1));
		v2 = v1;
	}
	else {
		v1 = complex(0, offsetN<real_t>(nV1));
		v2 = complex(0, offsetN<real_t>(nV2));
	}

	using PolyMoebius = Moebius<polynomial>;

	// Create generator matrices out of polynomials
	PolyMoebius genA = {
		{ 0.0 }, { 1.0 },
		{ -1.0 }, { 0.0, 1.0 } // here is our variable
	};
	PolyMoebius genB1 = genA * PolyMoebius{
		{ 1.0 }, { v1 },
		{ 0.0 }, { 1.0 }
	};
	PolyMoebius genB2 = genA * PolyMoebius{
		{ 1.0 }, { v2 },
		{ 0.0 }, { 1.0 }
	};

	// Generate sequence matrix
	PolyMoebius M;

	auto applyA = [&](PolyMoebius& m) {
		m = m * genA;
	};
	auto applyB = [&](PolyMoebius& m) {
		static bool toggle = false;
		m = m * ((toggle = !toggle) ? genB1 : genB2);
	};

	// This is basically a line drawing algorithm
	if (nA >= nB) {
		int x = 0;
		for (int i = 0; i < nA; ++i) {
			applyA(M);
			x += nB;
			if (x >= nA) {
				applyB(M);
				x -= nA;
			}
		}
	}
	else {
		int y = 0;
		for (int i = 0; i < nB; ++i) {
			applyB(M);
			y += nA;
			if (y >= nB) {
				applyA(M);
				y -= nB;
			}
		}
	}

	// The generated sequence matrix also happens to be
	// an elliptic Moebius transformation with a period = n,
	// so its trace (M.a + M.d) must be equal to 2cos(pi/n).
	// This corresponds to the equation:
	//
	// tr(M) = 2cos(pi/n)
	//
	// and since members of M are polynomials, the whole equation
	// becomes a polynomial that we're going to solve:
	//
	// M.a + M.d - 2cos(pi/n) = 0

	polynomial P = M.a + M.d - offsetN<real_t>(nN);

	// We are looking for roots closest to (2 + 0i),
	// a small imaginary value is added to improve root finding
	std::vector<complex> roots;
	if (method == 0) {
		roots = findRoots_Newton(P, complex(2.0, 0.1), findAllRoots);
	}
	else {
		roots = findRoots_SkowronGould(P, complex(2.0, 0.1), findAllRoots);
	}

	for (auto& root : roots) {
		if (root.real() < 0)
			root = -root;
		if (fabs(root.imag()) < 1e-10)
			root.imag(0);
	}

	return roots;
}
