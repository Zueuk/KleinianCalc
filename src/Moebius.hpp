#pragma once

template <typename complexT>
class Moebius
{
public:
	complexT a, b, c, d;

	Moebius() : a(1), b(0), c(0), d(1) { }

	Moebius(complexT a, complexT b, complexT c, complexT d) : a(a), b(b), c(c), d(d) { }

	Moebius operator *(const Moebius& other) const {
		return Moebius{
			a * other.a + b * other.c,
			a * other.b + b * other.d,
			c * other.a + d * other.c,
			c * other.b + d * other.d
		};
	}

	Moebius inverse() const { return Moebius{ d, -b, -c, a }; }

	complexT apply(complexT z) const { return (a * z + b) / (c * z + d); }
};
