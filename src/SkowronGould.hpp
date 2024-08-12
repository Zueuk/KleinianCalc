#pragma once

#include <boost/multiprecision/cpp_complex.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/tools/polynomial.hpp>

/// @brief Based on the algorithm published in the paper https://arxiv.org/abs/1203.1034
/// "General Complex Polynomial Root Solver and Its Further Optimization for Binary Microlenses"
/// by J. Skowron, A. Gould (Department of Astronomy, Ohio State University)
/// Submitted on 5 Mar 2012
///
/// Adapted for use with boost::multiprecision and polynomials by Peter Sdobnov

namespace SkowronGould {

/// @brief Solves quadratic equation.
///
/// @param[in]  poly   A vector of size = 3, holding the quadratic equation coefficients.
/// @param[out] out_x0 First root of the equation.
/// @param[out] out_x1 Second root of the equation.
template <typename complexT>
void solve_quadratic(const std::vector<complexT>& poly, complexT& out_x0, complexT& out_x1) {
	assert(poly.size() == 2 + 1);

	const complexT& a = poly[2];
	const complexT& b = poly[1];
	const complexT& c = poly[0];

	complexT sqrtD = sqrt(b*b - 4.0 * a*c);

	if (real(conj(b) * sqrtD) >= 0)
		out_x0 = -0.5 * (b + sqrtD);
	else
		out_x0 = -0.5 * (b - sqrtD);

	if (out_x0 == 0.0)
		out_x1 = 0.0;
	else { // Viete's formula
		out_x1 = c / out_x0;
		out_x0 = out_x0 / a;
	}
}

/// @brief Solves cubic equation using Lagrange's method.
///        See: http://en.wikipedia.org/wiki/Cubic_function Lagrange's method
///
/// @param[in]  poly   A vector of size = 4, holding the cubic equation coefficients.
/// @param[out] out_x0 First root of the equation.
/// @param[out] out_x1 Second root of the equation.
/// @param[out] out_x2 Third root of the equation.
template <typename complexT>
void solve_cubic(const std::vector<complexT>& poly, complexT& out_x0, complexT& out_x1, complexT& out_x2) {
	assert(poly.size() == 3 + 1);

	complexT E1 = -poly[2] / poly[3]; // x0 + x1 + x2
	complexT E2 =  poly[1] / poly[3]; // x0*x1 + x1*x2 + x2*x0
	complexT E3 = -poly[0] / poly[3]; // x0*x1*x2

	complexT s0 = E1;
	complexT E12 = E1 * E1;
	complexT A = 2.0 * E1 * E12 - 9.0 * E1 * E2 + 27.0 * E3;
	complexT B = E12 - 3.0 * E2;

	// quadratic equation z^2 - A * z + B^3 where roots are equal to s1^3 and s2^3
	complexT A2 = A * A;
	complexT sqrtD = sqrt(A2 - 4.0 * (B * B * B));

	complexT s1, s2;
	if (real(conj(A) * sqrtD) >= 0.0) // scalar product to decide the sign yielding bigger magnitude
		s1 = cbrt(0.5 * (A + sqrtD));
	else
		s1 = cbrt(0.5 * (A - sqrtD));

	if (s1 == 0.0)
		s2 = 0.0;
	else
		s2 = B / s1;

	using realT = typename complexT::valut_type;
	constexpr auto one_third = boost::math::constants::third<realT>();
	constexpr auto zeta1 = complex(-0.5, 0.5 * boost::math::constants::root_three<realT>());
	constexpr auto zeta2 = complex(-0.5, -0.5 * boost::math::constants::root_three<realT>());

	out_x0 = one_third * (s0 + s1 + s2);
	out_x1 = one_third * (s0 + s1 * zeta2 + s2 * zeta1);
	out_x2 = one_third * (s0 + s1 * zeta1 + s2 * zeta2);
}

namespace detail {

constexpr int MaxIters = 80;

enum class SolvingMethod {
	Newton = 0,
	SecondOrderGeneral = 1,
	Laguerre = 2,
};

constexpr double Frac_Jumps[] = { // some arbitrary numbers for 'random' jumps
	0.64109297, 0.91577881, 0.25921289, 0.50487203, 0.08177045,
	0.13653241, 0.30616200, 0.37794326, 0.04618805, 0.75132137
};
constexpr int Frac_Jumps_Count = (int)std::size(Frac_Jumps);
constexpr int Frac_Jump_Every = 10;
constexpr double Frac_Err = 2e-15;

const inline struct ComplexJumpsInitStruct { // to statically initialize 'random' jump points
	static inline std::complex<double> Jumps[Frac_Jumps_Count];

	ComplexJumpsInitStruct() {
		constexpr auto TwoPi = boost::math::constants::two_pi<double>();
		for (int i = 0; i < Frac_Jumps_Count; ++i)
			Jumps[i] = std::exp(std::complex<double>(0.0, Frac_Jumps[i] * TwoPi));
	}
	inline std::complex<double> operator [](int n) const { return Jumps[n]; }
} Complex_Jumps;

/// @brief Finds one root of a complex polynomial using Newton's method.
///
/// @details Calculates simplified Adams' stopping criterion for the value of
/// the polynomial once per 10 iterations after initial iteration.
/// This is done to speed up calculations when polishing roots that are known
/// preety well, and stopping criterion does significantly change in their
/// neighborhood.
///
/// Remember to initialize 'root' to some initial guess.
/// Do not initilize 'root' to point (0,0) if the polynomial
/// coefficients are strictly real, because it will make going
/// to imaginary roots impossible.
///
/// @param[in]     poly Vector of polynomial coefs.
/// @param[in,out] root input: Initial guess for the root value.
///                    output: Computed root value.
template <typename complexT>
void solve_1_newton_spec(const std::vector<complexT>& poly, complexT& root) {
	auto degree = (long)poly.size() - 1;

	using realT = typename complexT::value_type;

	bool good_to_go = false;

	for (int i = 1; i <= MaxIters; i++) {
		// prepare stoping criterion
		realT stopping_crit2; // will be initialized on the first iteration

		// calculate value of polynomial and its first two derivatives
		complexT p = poly[degree];
		complexT dp = 0.0;
		if (i % 10 == 1) { // calculate stopping criterion every tenth iteration
			realT e_k = abs(poly[degree]);
			realT absroot = abs(root);
			for (long k = degree - 1; k >= 0; k--) {
				dp = p + dp * root;
				p = poly[k] + p * root;
				// b_k
				// Adams, Duane A., 1967, "A stopping criterion for polynomial root finding",
				// Communications of ACM, Volume 10 Issue 10, Oct. 1967, p. 655
				// ftp://reports.stanford.edu/pub/cstr/reports/cs/tr/67/55/CS-TR-67-55.pdf
				// Eq. 8
				e_k = absroot * e_k + abs(p);
			}
			stopping_crit2 = pow(Frac_Err * e_k, 2);
		}
		else {
			// calculate value of the polynomial and its derivative using Horner's method
			// see for eg. Numerical Recipes Sec. 5.3 how to evaluate polynomials and derivatives
			for (long k = degree - 1; k >= 0; k--) {
				dp = p + dp * root;
				p = poly[k] + p * root;
			}
		}

		realT abs2p = norm(p); // re^2 + im^2
		if (abs2p == 0.0) // value of the polynomial == 0 means that we found the root
			return;

		if (abs2p < stopping_crit2) { // simplified a little Eq. 10 of Adams 1967
			if (dp == 0.0)
				return; // if we have problem with zero, but we are close to the root, just accept

			// do additional iteration if we are less than 10x from stopping criterion
			if (abs2p < 0.01 * stopping_crit2)
				return; // return immediatley because we are at very good place
			else
				good_to_go = true; // do one iteration more
		}
		else
			good_to_go = false; // reset if we are outside the zone of the root

		complexT dx = (dp != 0.0) ? // avoid division by zero
			p / dp :
			dx = (abs(root) + 1.0) * (complexT)Complex_Jumps[i % Frac_Jumps_Count]; // make a 'random' jump

		complexT newroot = root - dx;
		if (newroot == root) // no change, return
			return;
		if (good_to_go) { // this was jump already after stopping criterion was met
			root = newroot;
			return;
		}

		if (i % Frac_Jump_Every == 0) { // decide whether to do a jump of modified length (to break cycles)
			double faq = Frac_Jumps[(i / Frac_Jump_Every - 1) % Frac_Jumps_Count];
			newroot = root - faq * dx;
		}
		root = newroot;
	}
}

/// @brief Find one root of a complex polynomial using Laguerre's, Second-order
/// general and Newton's method - depending on the value of function F, which
/// is a combination of second derivative, first derivative and
/// value of polynomial [F = -(p"*p)/(p'p')].
///
/// @details The function has 3 modes of operation. It starts with Laguerre's method,
/// and continues until F < 0.50, at which point, it switches to SG method (see
/// the paper). While in the first two modes, the stopping criterion is
/// calculated once per every iteration.
/// Switch to the last mode - Newton's method - happens when F becomes < 0.05.
/// In this mode, stopping criterion is calculated only once at the beginning,
/// under an assumption that we are already very close to the root.
///
/// If there are more than 10 iterations in Newton's mode, it means that in
/// fact we were far from the root, so we go back to Laguerre's method.
///
/// Remember to initialize 'root' to some initial guess or to (0, 0)
///
/// @param[in]     poly Vector of polynomial coefs.
/// @param[in,out] root input: Initial guess for the value of a root.
///                    output: Computed root value.
/// @param mode This should be by default = 2. However if you choose to start
///             with SG method put 1 instead. Zero will cause the routine to
///             start with Newton for first 10 iterations, and then go back to
///             mode 2.
/// @returns true when a root is found,
///          false if reaches maximum number of iterations
template <typename complexT>
bool solve_1_laguerre2newton(const std::vector<complexT>& poly, complexT& root, SolvingMethod mode) {
	using realT = typename complexT::value_type;

	auto degree = poly.size() - 1;

	int iter = 0; // this iteration counter looks unused

	int i;
	int j = 1;
	bool good_to_go = false;

	realT stopping_crit2; // will be initialized on the first iteration

	for (;;) { // infinite loop to cycle through solving methods

		//
		// Laguerre method
		//

		if (mode >= SolvingMethod::Laguerre) {
			realT one_nth = realT(1.0) / degree;
			realT n_1_nth = (degree - 1) * one_nth;
			realT two_n_div_n_1 = 2.0 / n_1_nth;

			for (i = 1; i <= MaxIters; i++) {
				// prepare the stoping criterion
				realT e_k = abs(poly[degree]);
				realT absroot = abs(root);
				// calculate value of polynomial and its first two derivatives
				complexT p = poly[degree];
				complexT dp = 0.0;
				complexT d2p_half = 0.0;
				// Horner's method
				for (auto k = degree; k >= 1; k--) {
					d2p_half = dp + d2p_half * root;
					dp = p + dp * root;
					p = poly[k - 1] + p * root;
					// b_k
					// Adams, Duane A., 1967, "A stopping criterion for polynomial root finding",
					// Communications of the ACM, Volume 10 Issue 10, Oct. 1967, p. 655
					// ftp://reports.stanford.edu/pub/cstr/reports/cs/tr/67/55/CS-TR-67-55.pdf
					// Eq. 8
					e_k = absroot * e_k + abs(p);
				}
				realT abs2p = norm(p);
				iter += 1;
				if (abs2p == 0) // value of the polynomial == 0 means that we found the root
					return true;

				stopping_crit2 = pow(Frac_Err * e_k, 2);
				if (abs2p < stopping_crit2) {
					// (simplified a little Eq. 10 of Adams 1967)
					// do additional iteration if we are less than 10x from stopping criterion
					if (abs2p < 0.01*stopping_crit2)
						return true; // ten times better than stopping criterion
									 // return immediately, because we are at very good place
					else
						good_to_go = true; //do one iteration more
				}
				else
					good_to_go = false; //reset if we are outside the zone of the root

				complexT denom = 0.0;
				complexT fac_netwon = 0.0;
				if (dp != 0.0) {
					fac_netwon = p / dp;
					complexT fac_extra = d2p_half / dp;
					complexT F_half = fac_netwon * fac_extra;

					double abs2_F_half = (double)norm(F_half);
					if (abs2_F_half <= 0.0625) {
						// F < 0.50, F/2 < 0.25
						// switch to a different method
						if (abs2_F_half <= 0.000625) {
							// F < 0.05, F/2 < 0.02
							// switch to Newton's method
							mode = SolvingMethod::Newton;
						}
						else {
							// switch to Second-order General method
							mode = SolvingMethod::SecondOrderGeneral;
						}
					}

					complexT denom_sqrt = sqrt(1.0 - two_n_div_n_1 * F_half);
					// Probably unnecessary check, the real part is usually non-negative here
					if (real(denom_sqrt) >= 0.0)
						denom = one_nth + n_1_nth * denom_sqrt;
					else
						denom = one_nth - n_1_nth * denom_sqrt;
				}
				complexT dx = (denom != 0.0) ? // avoid division by zero
					fac_netwon / denom :
					(abs(root) + 1.0) + (complexT)Complex_Jumps[i % Frac_Jumps_Count]; // make a 'random' jump

				complexT newroot = root - dx;
				if (newroot == root) // no change, return
					return true;
				if (good_to_go) { // this was jump already after stopping criterion was met
					root = newroot;
					return true;
				}

				if (mode != SolvingMethod::Laguerre) {
					root = newroot;
					j = i + 1; //remember iteration index
					break; // go to Newton's or SG
				}

				if ((i % Frac_Jump_Every) == 0) { // decide whether to do a jump of modified length to break cycles
					double faq = Frac_Jumps[((i / Frac_Jump_Every - 1) % Frac_Jumps_Count)];
					newroot = root - faq * dx; // do jump of some semi-random length (0 < faq < 1)
				}
				root = newroot;
			}

			if (i >= MaxIters)
				return false;
		}

		//
		// Second-order general method
		//

		if (mode == SolvingMethod::SecondOrderGeneral) {

			for (i = j; i <= MaxIters; i++) {
				// calculate value of polynomial and its first two derivatives
				complexT p = poly[degree];
				complexT dp = 0.0;
				complexT d2p_half = 0.0;
				if ((i - j) % 10 == 0) {
					// prepare stopping criterion
					realT e_k = abs(poly[degree]);
					realT absroot = abs(root);
					// Horner's method
					for (auto k = degree; k >= 1; k--) {
						d2p_half = dp + d2p_half * root;
						dp = p + dp * root;
						p = poly[k - 1] + p * root;
						// b_k
						// Adams, Duane A., 1967, "A stopping criterion for polynomial root finding",
						// Communications of the ACM, Volume 10 Issue 10, Oct. 1967, p. 655
						// ftp://reports.stanford.edu/pub/cstr/reports/cs/tr/67/55/CS-TR-67-55.pdf
						// Eq. 8
						e_k = absroot * e_k + abs(p);
					}
					stopping_crit2 = pow(Frac_Err * e_k, 2);
				}
				else {
					// Horner's method
					for (auto k = degree; k >= 1; k--) {
						d2p_half = dp + d2p_half * root;
						dp = p + dp * root;
						p = poly[k - 1] + p * root; // b_k
					}
				}
				realT abs2p = norm(p);
				iter += 1;
				if (abs2p == 0.0) // value of the polynomial == 0 means that we found the root
					return true;

				if (abs2p < stopping_crit2) { // (simplified a little Eq. 10 of Adams 1967)
					if (dp == 0.0)
						return true;
					// do additional iteration if we are less than 10x from stopping criterion
					if (abs2p < 0.01*stopping_crit2)
						return true; // ten times better than stopping criterion
					else
						good_to_go = true; // do one iteration more
				}
				else
					good_to_go = false; // reset if we are outside the zone of the root

				complexT dx;
				if (dp == 0.0) // avoid division by zero
					dx = (abs(root) + 1.0) * (complexT)Complex_Jumps[i % Frac_Jumps_Count]; // make a 'random' jump
				else {
					complexT fac_netwon = p / dp;
					complexT fac_extra = d2p_half / dp;
					complexT F_half = fac_netwon * fac_extra;

					double abs2_F_half = (double)norm(F_half);
					if (abs2_F_half <= 0.000625) {
						// F < 0.05, F/2 < 0.025
						mode = SolvingMethod::Newton; // switch to Newton's method after the jump
					}
					dx = fac_netwon * (1.0 + F_half); // SG
				}

				complexT newroot = root - dx;
				if (newroot == root)
					return true; // no change, return
				if (good_to_go) {
					root = newroot; // this was jump already after stopping criterion was met
					return true;
				}

				if (mode != SolvingMethod::SecondOrderGeneral) {
					root = newroot;
					j = i + 1; // remember iteration number
					break; // switch to Newton's method
				}
				if ((i % Frac_Jump_Every) == 0) {
					// make a jump to break cycle
					newroot = root - dx * Frac_Jumps[(i / Frac_Jump_Every - 1) % Frac_Jumps_Count];
				}
				root = newroot;
			}
			if (i >= MaxIters)
				return false;
		}

		//
		// Newton's method
		//

		if (mode == SolvingMethod::Newton) {
			// Do only 10 iterations the most then go back to Laguerre
			for (i = j; i <= j + 10; i++) {
				// calculate value of the polynomial and its derivative
				complexT p = poly[degree];
				complexT dp = 0.0;
				if (i == j) { // Calculating stopping criterion only at the beginning
					realT e_k = abs(poly[degree]);
					realT absroot = abs(root);
					for (auto k = degree; k >= 1; k--) {
						dp = p + dp * root;
						p = poly[k - 1] + p * root;
						e_k = absroot * e_k + abs(p);
					}
					stopping_crit2 = pow(Frac_Err * e_k, 2);
				}
				else {
					for (auto k = degree; k >= 1; k--) {
						dp = p + dp * root;
						p = poly[k - 1] + p * root;
					}
				}
				realT abs2p = norm(p);
				iter += 1;
				if (abs2p == 0.0) // value of the polynomial == 0 means that we found the root
					return true;

				if (abs2p < stopping_crit2) {
					if (dp == 0.0)
						return true;
					// do additional iteration if we are less than 10x from stopping criterion
					if (abs2p < 0.01 * stopping_crit2)
						return true; // return immediately since we are close enough
					else
						good_to_go = true; // do one more iteration
				}
				else
					good_to_go = false;

				complexT dx = (dp != 0.0) ? // avoid division by zero
					dx = p / dp :
					dx = (abs(root) + 1.0) * (complexT)Complex_Jumps[i % Frac_Jumps_Count]; // make a 'random' jump
				complexT newroot = root - dx;
				if (newroot == root)
					return true;
				if (good_to_go) {
					root = newroot;
					return true;
				}
				root = newroot;
			}
			if (iter >= MaxIters)
				return false;

			mode = SolvingMethod::Laguerre; // Could not converge with 10 steps of Newton, switch back to Laguerre's method
		}

	} // end of infinite loop
}

/// @brief Finds one root of a complex polynomial using Laguerre's method.
/// In every loop it calculates simplified Adams' stopping criterion for
/// the value of the polynomial.
///
/// @param[in]     poly Vector of polynomial coefs.
/// @param[in,out] root input: Initial guess for the root value.
///                    output: Computed root value.
template <typename complexT>
void solve_1_laguerre(const std::vector<complexT>& poly, complexT& root) {
	using realT = typename complexT::value_type;

	auto degree = (long)poly.size() - 1;

	bool good_to_go = false;
	realT one_nth = realT(1.0) / degree;
	realT n_1_nth = (degree - 1.0) * one_nth;
	realT two_n_div_n_1 = 2.0 / n_1_nth;

	for (int i = 1; i <= MaxIters; i++) {
		realT e_k = abs(poly[degree]); // Preparing stopping criterion
		realT absroot = abs(root);
		// Calculate values of the polynomial and its first and second derivatives
		complexT p = poly[degree];
		complexT dp = 0.0;
		complexT d2p_half = 0.0;
		for (long k = degree - 1; k >= 0; k--) {
			d2p_half = dp + d2p_half * root;
			dp = p + dp * root;
			p = poly[k] + p * root;
			// b_k
			// Adams, Duane A., 1967, "A stopping criterion for polynomial root finding",
			// Communications of the ACM, Volume 10 Issue 10, Oct. 1967, p. 655
			// ftp://reports.stanford.edu/pub/cstr/reports/cs/tr/67/55/CS-TR-67-55.pdf
			// Eq. 8
			e_k = absroot * e_k + abs(p);
		}

		realT abs2p = norm(p);
		if (abs2p == 0) // value of the polynomial == 0 means that we found the root
			return;

		realT stopping_crit2 = pow(Frac_Err * e_k, 2);
		if (abs2p < stopping_crit2) {
			// (simplified a little Eq. 10 of Adams 1967)
			// do additional iteration if we are less than 10x from stopping criterion
			if (abs2p < 0.01 * stopping_crit2) {
				return; // we are close enough
			}
			else
				good_to_go = true;
		}
		else
			good_to_go = false;

		complexT denom = 0.0;
		complexT fac_newton = 0.0;

		if (dp != 0.0) {
			fac_newton = p / dp;
			complexT fac_extra = d2p_half / dp;
			complexT F_half = fac_newton * fac_extra;
			complexT denom_sqrt = sqrt(1.0 - two_n_div_n_1 * F_half);
			// Probably unnecessary check, the real part is usually non-negative here
			if (real(denom_sqrt) >= 0.0)
				denom = one_nth + n_1_nth * denom_sqrt;
			else
				denom = one_nth - n_1_nth * denom_sqrt;
		}

		complexT dx = (denom != 0.0) ? // avoid division by zero
			fac_newton / denom :
			(absroot + 1.0) * (complexT)Complex_Jumps[i % Frac_Jumps_Count]; // make a 'random' jump

		complexT newroot = root - dx;
		if (newroot == root)
			return; // no change, return
		if (good_to_go) {
			root = newroot;
			return;
		}
		if (i % Frac_Jump_Every == 0) {
			// make a jump to break cycles
			newroot = root - dx * Frac_Jumps[(i / Frac_Jump_Every - 1) % Frac_Jumps_Count];
		}
		root = newroot;
	}

	// failed to find the root, perhaps we should return false or something
}

/// @brief Divides the polynomial by (x - root)
template<typename complexT>
void divide_poly_by_root(/* in, out */ std::vector<complexT>& poly, const complexT& root) {
	auto n = (long)poly.size() - 1;
	complexT coef = poly[n];
	for (long i = n - 1; i >= 0; i--) {
		complexT prev = poly[i];
		poly[i] = coef;
		coef = prev + root * coef;
	}
	poly.pop_back();
}

} // namespace detail

/// @brief Finds all roots of the complex polynomial.
///
/// @param[in]     poly Vector of polynomial coefs.
///                     Highest order coef (last vector item) must be nonzero.
/// @param[in,out] root input: Initial guess for the root value.
///                    output: Computed root value.
/// @param       polish Perform an extra step for more precise results.
template <typename complexT>
void solve(const std::vector<complexT>& poly, std::vector<complexT>& roots, bool polish = true) {
	auto degree = (long)poly.size() - 1;

	// If `roots` has some values, we will use them as initial guess
	// the rest will be zero-initialized
	roots.resize(degree);

	// Check trivial cases first
	if (degree <= 2) {
		if (degree == 2)
			solve_quadratic(poly, roots[0], roots[1]);
		else if (degree == 1)
			roots[0] = -poly[0] / poly[1];
		return;
	}

	auto workPoly = poly; // work on a copy, keep the original for 'polishing'

	for (auto n = degree; n >= 3; n--) {
		auto& currentRoot = roots[n - 1];

		bool success = detail::solve_1_laguerre2newton<complexT>(workPoly, currentRoot, detail::SolvingMethod::Laguerre);
		if (!success) {
			currentRoot = 0.0;
			detail::solve_1_laguerre(workPoly, currentRoot);
		}

		// divide the polynomial by (x - root)
		complexT coef = workPoly[n];
		for (long i = n - 1; i >= 0; i--) {
			complexT prev = workPoly[i];
			workPoly[i] = coef;
			coef = prev + currentRoot * coef;
		}
		workPoly.pop_back(); // not really necessary
	}

	// At this point, degree of the polynomial == 2,
	// we can solve it as a quadratic eq to find the last two roots
	solve_quadratic(workPoly, roots[0], roots[1]);

	if (polish) {
		// 'polish' roots with full polynomial
		for (size_t n = 0; n < roots.size(); n++)
			detail::solve_1_newton_spec(poly, roots[n]);
	}
}

template <typename complexT>
using polynomial = boost::math::tools::polynomial<complexT>;

template <typename complexT>
inline void solve(const polynomial<complexT>& poly, std::vector<complexT>& roots, bool polish_roots_after = true) {
	solve(poly.data(), roots, polish_roots_after);
}

template <typename complexT>
inline void solve_one(polynomial<complexT>& poly, /* in, out */ complexT& root, bool polish_roots_after = true) {
	auto degree = poly.size() - 1;
	if (degree == 2) {
		complexT root1, root2;
		solve_quadratic(poly.data(), root1, root2);
		if (norm(root1) >= norm(root2))
			root = root1;
		else
			root = root2;
	}
	else if (degree == 1) {
		root = -poly[0] / poly[1];
	}
	else {
		detail::solve_1_laguerre2newton<complexT>(poly.data(), root, detail::SolvingMethod::Laguerre);
		if (polish_roots_after)
			detail::solve_1_newton_spec<complexT>(poly.data(), root);
	}
	detail::divide_poly_by_root(poly.data(), root);
}

} // namespace SkowronGould
