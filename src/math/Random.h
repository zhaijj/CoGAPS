#ifndef __COGAPS_RANDOM_H__
#define __COGAPS_RANDOM_H__

#include "../Archive.h"

#include <fstream>
#include <stdint.h>
#include <vector>

// used for seeding individual rngs
class Xoroshiro128plus
{
public:

    void seed(uint64_t);
    uint64_t next();

private:

    uint64_t mState[2];
    void warmup();

    friend Archive& operator<<(Archive &ar, Xoroshiro128plus &gen);
    friend Archive& operator>>(Archive &ar, Xoroshiro128plus &gen);
};

// PCG random number generator
class GapsRng
{
public:

    GapsRng();

    float uniform();
    float uniform(float a, float b);

    uint32_t uniform32();
    uint32_t uniform32(uint32_t a, uint32_t b);

    uint64_t uniform64();
    uint64_t uniform64(uint64_t a, uint64_t b);

    int poisson(float lambda);
    float exponential(float lambda);

    float inverseNormSample(float a, float b, float mean, float sd);
    float inverseGammaSample(float a, float b, float mean, float sd);
   
private:

    uint64_t mState;

    uint32_t next();
    void advance();
    uint32_t get() const;
};

namespace gaps
{
    namespace random
    {
        void setSeed(uint32_t seed);
        GapsRng getRng();

        void save();
        void load();
    }
}

#endif // __COGAPS_RANDOM_H__
