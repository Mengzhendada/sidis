#include "sidis/math.hpp"

#include <cmath>
#include <limits>

using namespace sidis;
using namespace sidis::math;

// Compute `log(|x|) log(|1 - x|)`
template<typename T>
static T log_log1m(T x) {
	if (!std::isfinite(x)) {
		if (std::isinf(x)) {
			return std::numeric_limits<T>::infinity();
		} else {
			return x;
		}
	} else if (x > 1.) {
		return std::log(x) * std::log(x - 1.);
	} else if (x == 1.) {
		return 0.;
	} else if (x > 0.5) {
		return std::log1p(x - 1.) * std::log(1. - x);
	} else if (x > 0.) {
		return std::log(x) * std::log1p(-x);
	} else {
		return std::log(-x) * std::log1p(-x);
	}
}

template<typename T>
static T dilog_impl(T x) {
	// Ensure we have as many digits of pi as can fit the `T` type.
	T pi = 3.141592653589793238462643383279502797479068098137295573004504331874296718662975536062731407582759857177734375;

	// To compute the dilogarithm, first bring the argument into the range
	// (-0.5, 0.5).
	T a = 1.;
	T b = 0.;
	if (!std::isfinite(x)) {
		if (std::isinf(x)) {
			return -std::numeric_limits<T>::infinity();
		} else {
			return x;
		}
	} else if (x > 2.) {
		a = -1.;
		b = sq(pi) / 3. - 0.5 * sq(std::log(x));
		x = 1. / x;
	} else if (x > 1.) {
		a = 1.;
		//b = sq(pi) / 6. + 0.5 * std::log(x) * (std::log(x / sq(1. - x));
		b = sq(pi) / 6. + 0.5 * sq(std::log(x)) - log_log1m(x);
		x = (x - 1.) / x;
	} else if (x > 0.5) {
		a = -1.;
		//b = sq(pi) / 6. - std::log(x) * std::log(1. - x);
		b = sq(pi) / 6. - log_log1m(x);
		x = 1. - x;
	} else if (x > -0.5) {
		a = 1.;
		b = 0.;
		x = x;
	} else if (x > -1.) {
		a = -1.;
		b = -0.5 * sq(std::log1p(-x));
		x = x / (x - 1.);
	} else {
		a = 1.;
		//b = -sq(pi) / 6. + 0.5 * std::log(1. - x) * std::log((1. - x) / sq(x));
		b = -sq(pi) / 6. + 0.5 * sq(std::log1p(-x)) - log_log1m(x);
		x = 1. / (1. - x);
	}

	// Size of the mantissa in base 2.
	double d = std::numeric_limits<T>::digits
		* std::log2(std::numeric_limits<T>::radix);
	// Number of terms needed to reach desired precision.
	unsigned n_max = (unsigned) std::ceil(
		(1.4 * d + 6. * (1. - std::log(d))) * std::abs(x) + 0.3 * d) + 1;
	T result = 0.;
	T numerator = x;
	for (unsigned n = 1; n < n_max; ++n) {
		result += numerator / (n * n);
		numerator *= x;
	}

	return a * result + b;
}

template<typename T>
static T trapezoid_impl(T (*f)(T), T a, T b, unsigned n) {
	T delta = (b - a) / n;
	T result = 0.5 * (f(a) + f(b));
	for (unsigned i = 1; i <= n - 1; ++i) {
		T x = i / n * (b - a) + a;
		result += f(x);
	}
	return result * delta;
}

float math::dilog(float x) {
	return dilog_impl<float>(x);
}
double math::dilog(double x) {
	return dilog_impl<double>(x);
}
long double math::dilog(long double x) {
	return dilog_impl<long double>(x);
}

float math::trapezoid(float (*f)(float), float a, float b, unsigned n) {
	return trapezoid_impl<float>(f, a, b, n);
}
double math::trapezoid(double (*f)(double), double a, double b, unsigned n) {
	return trapezoid_impl<double>(f, a, b, n);
}
long double math::trapezoid(
		long double (*f)(long double),
		long double a,
		long double b,
		unsigned n) {
	return trapezoid_impl<long double>(f, a, b, n);
}

