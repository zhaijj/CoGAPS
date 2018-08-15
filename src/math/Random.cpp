#include "Math.h"
#include "Random.h"
#include "../GapsAssert.h"

#include <stdint.h>

const float maxU32AsFloat = static_cast<float>(std::numeric_limits<uint32_t>::max());
const double maxU32AsDouble = static_cast<double>(std::numeric_limits<uint32_t>::max());

static Xoroshiro128plus seeder;

static uint64_t rotl(const uint64_t x, int k)
{
    return (x << k) | (x >> (64 - k));
}

void GapsRng::setSeed(uint64_t seed)
{
    seeder.seed(seed);
}

void GapsRng::save(Archive &ar)
{
    ar << seeder;
}

void GapsRng::load(Archive &ar)
{
    ar >> seeder;
}

void Xoroshiro128plus::seed(uint64_t seed)
{
    mState[0] = seed|1;
    mState[1] = seed|1;
    warmup();
}

uint64_t Xoroshiro128plus::next()
{
    uint64_t result = 0;
    #pragma omp critical(RngCreation)
    {
        const uint64_t s0 = mState[0];
        uint64_t s1 = mState[1];
        result = s0 + s1;
        s1 ^= s0;
        mState[0] = rotl(s0, 24) ^ s1 ^ (s1 << 16); // a, b
        mState[1] = rotl(s1, 37); // c
    }
    return result;
}

void Xoroshiro128plus::warmup()
{
    for (unsigned i = 0; i < 5000; ++i)
    {
        next();
    }
}

Archive& operator<<(Archive &ar, Xoroshiro128plus &gen)
{
    ar << gen.mState[0] << gen.mState[1];
    return ar;
}

Archive& operator>>(Archive &ar, Xoroshiro128plus &gen)
{
    ar >> gen.mState[0] >> gen.mState[1];
    return ar;
}

GapsRng::GapsRng() : mState(seeder.next()) {}

uint32_t GapsRng::next()
{
    advance();
    return get();
}

void GapsRng::advance()
{
    mState = mState * 6364136223846793005ull + (54u|1);
}

uint32_t GapsRng::get() const
{
    uint32_t xorshifted = ((mState >> 18u) ^ mState) >> 27u;
    uint32_t rot = mState >> 59u;
    return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
}

double GapsRng::uniformd()
{
    return static_cast<double>(uniform32()) / maxU32AsDouble;
}

float GapsRng::uniform()
{
    return static_cast<float>(uniform32()) / maxU32AsFloat;
}

float GapsRng::uniform(float a, float b)
{
    return uniform() * (b - a) + a;
}

uint32_t GapsRng::uniform32()
{
    return next();
}

uint32_t GapsRng::uniform32(uint32_t a, uint32_t b)
{
    uint32_t range = b - a;
    if (range == 0)
    {
        return a;
    }

    uint32_t x = uniform32();
    uint32_t iPart = std::numeric_limits<uint32_t>::max() / range;
    while (x >= range * iPart)
    {
        x = uniform32();
    }
    return x / iPart + a;
}

uint64_t GapsRng::uniform64()
{
    uint64_t high = (static_cast<uint64_t>(uniform32()) << 32) & 0xFFFFFFFF00000000ull;
    uint64_t low = uniform32();
    return high | low;
}

uint64_t GapsRng::uniform64(uint64_t a, uint64_t b)
{
    uint64_t range = b - a;
    if (range == 0)
    {
        return a;
    }

    uint64_t x = uniform64();
    uint64_t iPart = std::numeric_limits<uint64_t>::max() / range;
    while (x >= range * iPart)
    {
        x = uniform64();
    }
    return x / iPart + a;
}

int GapsRng::poisson(double lambda)
{
    if (lambda <= 5.0)
    {
        return poissonSmall(lambda);
    }
    else
    {
        return poissonLarge(lambda);
    }
}

// lambda <= 5
int GapsRng::poissonSmall(double lambda)
{
    int x = 0;
    double p = uniformd();
    double cutoff = std::exp(-lambda);
    while (p >= cutoff)
    {
        p *= uniformd();
        ++x;
    }
    return x;
}

// lambda > 5
int GapsRng::poissonLarge(double lambda)
{
    double c = 0.767 - 3.36 / lambda;
    double beta = gaps::pi_double / sqrt(3.0 * lambda);
    double alpha = beta * lambda;
    double k = std::log(c) - lambda - std::log(beta);

    while(true)
    {
        double u = uniformd();
        double x = (alpha - std::log((1.0 - u) / u)) / beta;
        double n = floor(x + 0.5);
        if (n < 0.0)
        {
            continue;
        }

        double v = uniformd();
        double y = alpha - beta * x;
        double w = 1.0 + std::exp(y);
        double lhs = y + std::log(v / (w * w));
        double rhs = k + n * std::log(lambda) - gaps::lgamma(n+1);
        if (lhs <= rhs)
        {
            return n;
        }
    }
}

float GapsRng::exponential(float lambda)
{
    return -1.f * std::log(uniform()) / lambda;
}


OptionalFloat GapsRng::truncNormal(float a, float b, float mean, float sd)
{
    float lower = gaps::p_norm(a, mean, sd);
    float upper = gaps::p_norm(b, mean, sd);

    if (lower > 0.95f || upper < 0.05f) // too far in tail
    {
        return OptionalFloat();
    }

    float u = uniform(lower, upper);
    while (u == 0.f || u == 1.f)
    {
        u = uniform(lower, upper);
    }
    float ret = gaps::q_norm(u, mean, sd);
    GAPS_ASSERT(ret >= a);
    GAPS_ASSERT(ret <= b);
    return ret;
}

float GapsRng::truncGammaUpper(float b, float shape, float scale)
{
    float upper = gaps::p_gamma(b, shape, scale);
    GAPS_ASSERT(upper > 0.f);
    float u = uniform(0.f, upper);
    while (u == 0.f || u == 1.f)
    {
        u = uniform(0.f, upper);
    }
    float ret = gaps::q_gamma(u, shape, scale);
    GAPS_ASSERT(ret <= b);
    return ret;
}