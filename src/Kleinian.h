#pragma once

#include <boost/multiprecision/cpp_complex.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/tools/polynomial.hpp>

#include "Moebius.hpp"

template <class realT>
constexpr realT offsetN(int n) { return (n == 2) ? 0 : 2*cos(boost::math::constants::pi<realT>() / n); }

class Kleinian
{
	// Using boost::multiprecision to solve higher order polynomials
	using complex = boost::multiprecision::cpp_complex_oct;
	using real_t = complex::value_type;

	// Using boost::polynomial class for convenience
	using polynomial = boost::math::tools::polynomial<complex>;

	std::vector<complex> findRoots_Newton(polynomial& P, complex guess, bool findAllRoots);
	std::vector<complex> findRoots_SkowronGould(polynomial& P, complex guess, bool findAllRoots);

	int nA, nB;
	int nN;
	int nV1, nV2;

public:
	Kleinian(int na, int nb, int nn, int nv1, int nv2)
		: nA(na), nB(nb)
		, nN(nn)
		, nV1(nv1), nV2(nv2)
	{ }

	std::vector<complex> solve(int method, bool findAllRoots);
};
