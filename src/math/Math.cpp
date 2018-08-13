#include "Math.h"

#define Q_GAMMA_THRESHOLD 0.000001f

float gaps::min(float a, float b)
{
    return a < b ? a : b;
}

unsigned gaps::min(unsigned a, unsigned b)
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

unsigned gaps::max(unsigned a, unsigned b)
{
    return a < b ? b : a;
}

uint64_t gaps::max(uint64_t a, uint64_t b)
{
    return a < b ? b : a;
}

float gaps::dgamma(float d, float shape, float scale)
{
    boost::math::gamma_distribution<> gam(shape, scale);
    return pdf(gam, d);
}

float gaps::pgamma(float p, float shape, float scale)
{
    boost::math::gamma_distribution<> gam(shape, scale);
    return cdf(gam, p);
}

float gaps::qgamma(float q, float shape, float scale)
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

float gaps::dnorm(float d, float mean, float sd)
{
    boost::math::normal_distribution<> norm(mean, sd);
    return pdf(norm, d);
}

float gaps::pnorm(float p, float mean, float sd)
{
    boost::math::normal_distribution<> norm(mean, sd);
    return cdf(norm, p);
}

float gaps::qnorm(float q, float mean, float sd)
{
    boost::math::normal_distribution<> norm(mean, sd);
    return quantile(norm, q);
}

float gaps::dnorm_fast(float d, float mean, float sd)
{
    return std::exp((d - mean) * (d - mean) / (-2.f * sd * sd))
        / std::sqrt(2.f * gaps::pi * sd * sd);
}
