#include "Math.h"

#define Q_GAMMA_THRESHOLD 0.000001f

#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/exponential.hpp>
#include <boost/math/distributions/gamma.hpp>
#include <boost/math/special_functions/gamma.hpp>

float gaps::min(float a, float b)
{
    return a < b ? a : b;
}

uint32_t gaps::min(uint32_t a, uint32_t b)
{
    return a < b ? a : b;
}

uint64_t gaps::min(uint64_t a, uint64_t b)
{
    return a < b ? a : b;
}

float gaps::max(float a, float b)
{
    return a < b ? b : a;
}

uint32_t gaps::max(uint32_t a, uint32_t b)
{
    return a < b ? b : a;
}

uint64_t gaps::max(uint64_t a, uint64_t b)
{
    return a < b ? b : a;
}

double gaps::lgamma(double x)
{
    return boost::math::lgamma(x);
}

float gaps::d_gamma(float d, float shape, float scale)
{
    boost::math::gamma_distribution<> gam(shape, scale);
    return pdf(gam, d);
}

float gaps::p_gamma(float p, float shape, float scale)
{
    boost::math::gamma_distribution<> gam(shape, scale);
    return cdf(gam, p);
}

float gaps::q_gamma(float q, float shape, float scale)
{
    if (q < Q_GAMMA_THRESHOLD)
    {
        return 0.f;
    }
    else
    {
        boost::math::gamma_distribution<> gam(shape, scale);
        return quantile(gam, q);
    }
}

float gaps::d_norm(float d, float mean, float sd)
{
    boost::math::normal_distribution<> norm(mean, sd);
    return pdf(norm, d);
}

float gaps::p_norm(float p, float mean, float sd)
{
    boost::math::normal_distribution<> norm(mean, sd);
    return cdf(norm, p);
}

float gaps::q_norm(float q, float mean, float sd)
{
    boost::math::normal_distribution<> norm(mean, sd);
    return quantile(norm, q);
}

float gaps::d_norm_fast(float d, float mean, float sd)
{
    return std::exp((d - mean) * (d - mean) / (-2.f * sd * sd))
        / std::sqrt(2.f * gaps::pi * sd * sd);
}
