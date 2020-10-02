#ifndef SIDIS_MATH_HPP
#define SIDIS_MATH_HPP

namespace sidis {
namespace math {

inline float sq(float x) {
	return x * x;
}
inline double sq(double x) {
	return x * x;
}
inline long double sq(long double x) {
	return x * x;
}

float sqrt1p_1m(float x);
double sqrt1p_1m(double x);
long double sqrt1p_1m(long double x);

float dilog(float x);
double dilog(double x);
long double dilog(long double x);

}
}

#endif
