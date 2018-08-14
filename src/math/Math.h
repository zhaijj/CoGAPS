#ifndef __COGAPS_MATH_H__
#define __COGAPS_MATH_H__

#include <stdint.h>
#include <string>
#include <sstream>

namespace gaps
{
    const float epsilon = 1.0e-10f;
    const float pi = 3.14159265358979323846264f;
    const double pi_double = 3.14159265358979323846264;

    float min(float a, float b);
    uint32_t min(uint32_t a, uint32_t b);
    uint64_t min(uint64_t a, uint64_t b);

    float max(float a, float b);
    uint32_t max(uint32_t a, uint32_t b);
    uint64_t max(uint64_t a, uint64_t b);

    double lgamma(double x);

    float dgamma(float d, float shape, float scale); // pdf
    float pgamma(float p, float shape, float scale); // cdf
    float qgamma(float q, float shape, float scale); // quantile/inverse cdf

    float dnorm(float d, float mean, float sd); // pdf
    float pnorm(float p, float mean, float sd); // cdf
    float qnorm(float q, float mean, float sd); // quantile/inverse cdf

    float dgamma_fast(float d, float shape, float scale);
    float pgamma_fast(float p, float shape, float scale);
    float qgamma_fast(float q, float shape, float scale);

    float dnorm_fast(float d, float mean, float sd);
    float pnorm_fast(float p, float mean, float sd);
    float qnorm_fast(float q, float mean, float sd);

    template <class T>
    std::string to_string(T a);
}

template <class T>
std::string gaps::to_string(T a)
{
    std::stringstream ss;
    ss << a;
    return ss.str();
}

#endif