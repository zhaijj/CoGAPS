#ifndef __COGAPS_RANDOM_H__
#define __COGAPS_RANDOM_H__

#include "../Archive.h"

#include <fstream>
#include <stdint.h>
#include <vector>

// put in cpp
static uint64_t rotl(const uint64_t x, int k)
{
    return (x << k) | (x >> (64 - k));
}

// used for seeding individial rngs
class Xoroshiro128plus
{
public:

    Xoroshiro128plus(uint64_t seed)
    {
        mState[0] = seed;
        mState[1] = mState[0];
        warmup();
    }

    uint64_t next()
    {
        const uint64_t s0 = mState[0];
        uint64_t s1 = mState[1];
        const uint64_t result = s0 + s1;
        s1 ^= s0;
        s[0] = rotl(s0, 24) ^ s1 ^ (s1 << 16); // a, b
        s[1] = rotl(s1, 37); // c
        return result;
    }

private:

    uint64_t mState[2];

    void warmup()
    {
        for (unsigned i = 0; i < 5000; ++i)
        {
            next();
        }
    }
};

// factory methods for seeding and constructing
/* USE
    gaps::random::seed()
    gaps::random::save()
    gaps::random::load()
    rng = gaps::random::getRng()
    rng.uniform(0,1)
    rng.normal(4, 0.1)
    rng.truncNormal(0, 1, 0.5, 1)
    gaps::math::q_gamma()
*/

// PCG random number generator
class GapsRng
{
public:

    GapsRng(uint64_t seed) : mState(seed) {}

    float uniform()
    {
        advance();
        get() / std::numeric_limits<uint32_t>::max();
    }

    float uniform(float a, float b)
    {
        return uniform() * (b - a) + a;
    }

    uint64_t uniform64()
    {
        advance();
        uint64_t high = get() << 32 & 0xFFFFFFFF000000000ul;
        advance();
        uint64_t low = get();
        return high | low;
    }
    
private:

    uint64_t mState;

    void advance()
    {
        mState = mState * 6364136223846793005ULL + (54u|1);
    }

    uint32_t get()
    {
        uint32_t xorshifted = ((mState >> 18u) ^ mState) >> 27u;
        uint32_t rot = mState >> 59u;
        return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
    }
};

namespace gaps
{
    namespace random
    {
        void setSeed(uint32_t seed);

        float uniform();
        float uniform(float a, float b);
        uint64_t uniform64();
        uint64_t uniform64(uint64_t a, uint64_t b);

        float exponential(float lambda);
        float normal(float mean, float var);
        int poisson(float lambda);

        float inverseNormSample(float a, float b, float mean, float sd);
        float inverseGammaSample(float a, float b, float mean, float sd);

        float d_gamma(float d, float shape, float scale);
        float p_gamma(float p, float shape, float scale);
        float q_gamma(float q, float shape, float scale);
        float d_norm(float d, float mean, float sd);
        float q_norm(float q, float mean, float sd);
        float p_norm(float p, float mean, float sd);

        void save(Archive &ar);
        void load(Archive &ar);
    } // namespace random
} // namespace gaps

#endif // __COGAPS_RANDOM_H__
