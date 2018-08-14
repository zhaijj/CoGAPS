#ifndef __COGAPS_RANDOM_H__
#define __COGAPS_RANDOM_H__

#include "../Archive.h"
#include "Math.h"

#include <fstream>
#include <stdint.h>
#include <vector>

struct OptionalFloat
{
    bool hasValue;
    float value;

    OptionalFloat() : hasValue(false), value(0.f) {}
    OptionalFloat(float f) : hasValue(true), value(f) {}
};

// used for seeding individual rngs
class Xoroshiro128plus
{
public:

    void seed(uint64_t seed);
    uint64_t next();

private:

    uint64_t mState[2];
    void warmup();

    friend Archive& operator<<(Archive &ar, Xoroshiro128plus &gen);
    friend Archive& operator>>(Archive &ar, Xoroshiro128plus &gen);
};

// PCG random number generator
// This is constructed with a seed pulled from the global seeder
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

    int poisson(double lambda);
    float exponential(float lambda);

    OptionalFloat truncNormal(float a, float b, float mean, float sd);
    float truncGammaUpper(float b, float shape, float scale);

    // these functions control the global seed
    static void setSeed(uint64_t seed);
    static void save(Archive &ar);
    static void load(Archive &ar);
   
private:

    uint64_t mState;

    uint32_t next();
    void advance();
    uint32_t get() const;

    double uniformd();
    int poissonSmall(double lambda);
    int poissonLarge(double lambda);
};

#endif // __COGAPS_RANDOM_H__
