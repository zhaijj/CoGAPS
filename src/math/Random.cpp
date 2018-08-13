#include "Random.h"

#include <stdint.h>

static Xoroshiro128plus seeder;

static uint64_t rotl(const uint64_t x, int k)
{
    return (x << k) | (x >> (64 - k));
}

void gaps::random::seed(uint64_t seed)
{
    seeder.seed(seed);
}

void gaps::random::save(Archive &ar)
{
    ar << seeder;
}

void gaps::random::load(Archive &ar)
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
    const uint64_t s0 = mState[0];
    uint64_t s1 = mState[1];
    const uint64_t result = s0 + s1;
    s1 ^= s0;
    s[0] = rotl(s0, 24) ^ s1 ^ (s1 << 16); // a, b
    s[1] = rotl(s1, 37); // c
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
    mState = mState * 6364136223846793005ULL + (54u|1);
}

uint32_t GapsRng::get() const
{
    uint32_t xorshifted = ((mState >> 18u) ^ mState) >> 27u;
    uint32_t rot = mState >> 59u;
    return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
}

float GapsRng::uniform()
{
    return uniform32() / std::numeric_limits<uint32_t>::max();
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
    uint32_t x = uniform32();
    uint32_t range = b - a;
    uint32_t iPart = std::numeric_limits<uint32_t>::max() / range;
    while (x >= range * iPart)
    {
        x = uniform32();
    }
    return x / iPart + a;
}

uint64_t GapsRng::uniform64()
{
    uint64_t high = uniform32() << 32 & 0xFFFFFFFF000000000ul;
    uint64_t low = uniform32();
    return high | low;
}

uint64_t GapsRng::uniform64(uint64_t a, uint64_t b)
{
    uint64_t x = uniform32();
    uint64_t range = b - a;
    uint64_t iPart = std::numeric_limits<uint64_t>::max() / range;
    while (x >= range * iPart)
    {
        x = uniform64();
    }
    return x / iPart + a;
}

int GapsRng::poisson(float lambda)
{
    int x = 0;
    float p = uniform();
    while (p >= std::exp(-lambda))
    {
        p *= uniform();
        ++x;
    }
    return x;        
}

float GapsRng::exponential(float lambda)
{
    return -1.f * std::log(uniform()) / lambda;
}


float GapsRng::inverseNormSample(float a, float b, float mean, float sd)
{
    float u = uniform(a, b);
    while (u == 0.f || u == 1.f)
    {
        u = uniform(a, b);
    }
    return gaps::qnorm(u, mean, sd);
}

float GapsRng::inverseGammaSample(float a, float b, float mean, float sd)
{
    float u = uniform(a, b);
    while (u == 0.f || u == 1.f)
    {
        u = uniform(a, b);
    }
    return gaps::qgamma(u, mean, sd);
}