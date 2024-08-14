#pragma once

#include <boost/multiprecision/cpp_complex.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/tools/polynomial.hpp>

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
	int n;
	int v1, v2;

	Kleinian() = default;

	std::vector<complex> solve(int method, bool findAllRoots) const;
};
