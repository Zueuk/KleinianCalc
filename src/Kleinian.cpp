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
		int n = (P.size() > 2) ? 2 : 1;
		std::vector<complex> roots(n, guess);
		for (auto& root : roots)
			SkowronGould::solve_one(P, root, true);
		return roots;
	}
	else {
		std::vector<complex> roots(P.degree(), guess);
		SkowronGould::solve(P, roots, true);
		return roots;
	}
}

std::vector<Kleinian::complex> Kleinian::solve(int method, bool findAllRoots) const {
	using PolyMoebius = Moebius<polynomial>;
	PolyMoebius Ma, Mb, Mc;

	if (isInf(h)) {
		Ma = {
			{ 0.0 }, { 1.0 },
			{ -1.0 }, { 0.0, 1.0 }
		};
		Mb = Ma * PolyMoebius{
			{ 1.0 }, { complex(0, offsetN<real>(v1)) },
			{ 0.0 }, { 1.0 }
		};
		Mc = Ma * PolyMoebius{
			{ 1.0 }, { complex(0, offsetN<real>(v2)) },
			{ 0.0 }, { 1.0 }
		};
	}
	else {
		TilngParameters<real> tiling(h, v1, v2);

		/* Building Ma using fixed points and poles, see
		https://en.wikipedia.org/wiki/M%C3%B6bius_transformation#Poles_of_the_transformation

		We are making a matrix with fp1 = -fp2, therefore Z_inf = -z_inf
		Then the matrix will look like
			[ Z_inf   fp²  ]
			[   1    Z_inf ]

		One of the tiling's centers of rotation is at the point = -ph, and
		we know that Ma applied to (1/ph) = -ph, that is
			Z_inf * 1/ph + fp²
			------------------ = -ph,
			   1/ph + Z_inf

		Solving this for Z_inf, we get
					-1 - fp²
			Z_inf = ---------,
					1/ph + ph

		Converted to polynomail with fp² as the variable
			Z_inf = -1/(1/ph + ph) + -1/(1/ph + ph)*fp²
		or,
			pZ = -1/(1/ph + ph) = -ph/(1 + ph²)
			Z_inf = pZ + pZ*fp²

		Now we can construct the matrix of polynomials to solve for fp² */

		real pZ = -tiling.ph/(1 + tiling.ph * tiling.ph);
		Ma = {
			{ pZ, pZ }, { 0.0, 1.0 },
			{ 1.0 },    { pZ, pZ }
		};

		static const PolyMoebius Disc2halfplane = {
			{ 1.0 }, { 1.0 },
			{-1.0 }, { 1.0 }
		};
		static const PolyMoebius Halfplane2disc = {
			{ 1.0 }, {-1.0 },
			{ 1.0 }, { 1.0 }
		};
		auto Rotation = [](real a) {
			return PolyMoebius {
				{ exp(complex(0, a)) }, { 0.0 },
				{ 0.0 }, { 1.0 }
			};
		};
		auto PoincareShift = [](complex p) {
			return PolyMoebius {
				{ 1.0 }, { p },
				{ conj(p) }, { 1.0 }
			};
		};

		Mb = Ma * PoincareShift(tiling.ph) * Rotation(-2*tiling.B1) * PoincareShift(-tiling.ph);
		Mc = Ma * PoincareShift(tiling.ph) * Rotation(-2*tiling.C1) * PoincareShift(-tiling.ph);
	}

	PolyMoebius M = createSolverMatrix(Ma, Mb, Mc);

	// We know that the generated matrix M is an elliptic Moebius transformation
	// with period = n, so its trace (M.a + M.d) must be equal to 2cos(pi/n)
	//
	// M.a + M.d = 2cos(pi/n)
	//
	// Since members of M are polynomials, this equation
	// becomes a polynomial that we're going to solve

	polynomial P;
	complex guess;
	bool (*filterPred)(const complex& z);
	bool (*sortPred)(const complex& a, const complex& b);

	if (isInf(h)) {
		P = M.a + M.d - offsetN<real>(n);
		// We are looking for the root closest to (2 + 0i),
		// a small imaginary value is added to improve root finding
		guess = complex(2.0, 0.01);
		// Roots must have their real value between 1 and 2
		filterPred = [](const complex& z) {
			auto re = abs(z.real());
			return re >= 0.0 && re <= 2.0;
		};
		// Sort by abs value, descending
		sortPred = [](const complex& a, const complex& b) {
			return norm(a) > norm(b);
		};
	}
	else {
		if (n == 2)
			P = M.a + M.d; // 2cos(pi/2) = 0
		else {
			// Here |det(M)| != 1, so the equation becomes
			//
			//  M.a + M.d
			// ------------ = 2cos(pi/n)
			// sqrt(det(M))
			//
			// To get a polynomial, we square it to get rid of sqrt(),
			// and then mulyiply by det(M) = (M.a * M.d - M.b * M.c)
			auto sqr = [](auto x){ return x * x; };
			P = sqr(M.a + M.d) - sqr(offsetN<real>(n)) * (M.a * M.d - M.b * M.c);
		}
		// We are looking for the root closest to zero,
		// a small imaginary value is added to improve root finding
		guess = complex(0.01, 0.01);
		// roots must have their abs value inside the unit disc
		filterPred = [](const complex& z) {
			return norm(z) <= 1.0;
		};
		// Sort by abs value, ascending
		sortPred = [](const complex& a, const complex& b) {
			return norm(a) < norm(b);
		};
	}

	if (P.size() == 0)
		return {};

	auto allRoots = (method == 0)
		? findRoots_Newton(P, guess, findAllRoots)
		: findRoots_SkowronGould(P, guess, findAllRoots);

	std::vector<complex> roots;
	std::copy_if(allRoots.begin(), allRoots.end(), std::back_inserter(roots), filterPred);

	std::sort(roots.begin(), roots.end(), sortPred);

	return roots;
}
