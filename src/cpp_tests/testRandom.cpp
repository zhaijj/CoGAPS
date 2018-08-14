#include "catch.h"
#include "../math/Math.h"
#include "../math/Random.h"

#define TEST_SQ(x) ((x) * (x))

template <typename T>
void testMeanVar(const std::string &type, const std::vector<T> &vals,
float mean, float var)
{
    printf("\nTesting %s\n", type.c_str());

    float sum = 0.f;
    for (unsigned i = 0; i < vals.size(); ++i)
    {
        sum += vals[i];
    }
    float sampleMean = sum / vals.size();
    REQUIRE(sampleMean == Approx(mean).epsilon(0.01f));

    float var_sum = 0.f;
    for (unsigned i = 0; i < vals.size(); ++i)
    {
        var_sum += TEST_SQ(vals[i] - sampleMean);
    }
    float sampleVar = var_sum / (vals.size() - 1);
    REQUIRE(sampleVar == Approx(var).epsilon(5.f));
}

TEST_CASE("Test Random.h - RNG Interface")
{
    GapsRng::setSeed(0);
    GapsRng rng;

    SECTION("Make sure uniform01 is working")
    {
        REQUIRE(rng.uniform() != rng.uniform());
    }

    SECTION("Test bounds")
    {
        for (unsigned i = 1; i < 10000; ++i)
        {
            // Uniform[0,1]
            REQUIRE(rng.uniform() >= 0.f);
            REQUIRE(rng.uniform() <= 1.f);

            // poisson
            REQUIRE(rng.poisson(i) >= 0);

            // Exponential
            float exp = rng.exponential(i);
            REQUIRE(exp >= 0.f);

            // Uniform[a,b]
            REQUIRE(rng.uniform(exp, 2.f * exp) >= exp);
            REQUIRE(rng.uniform(exp, 2.f * exp) <= 2.f * exp);

            // uniform32[a,b]
            uint32_t a32 = rng.uniform32();
            uint32_t b32 = rng.uniform32();
            uint32_t lower32 = gaps::min(a32,b32);
            uint32_t upper32 = gaps::max(a32,b32);
            REQUIRE(rng.uniform32(lower32, upper32) >= lower32);
            REQUIRE(rng.uniform32(lower32, upper32) <= upper32);

            // uniform64[a,b]
            uint64_t a64 = rng.uniform64();
            uint64_t b64 = rng.uniform64();
            uint64_t lower64 = gaps::min(a64,b64);
            uint64_t upper64 = gaps::max(a64,b64);
            REQUIRE(rng.uniform64(lower64, upper64) >= lower64);
            REQUIRE(rng.uniform64(lower64, upper64) <= upper64);

            // truncNormal
            float mean = rng.uniform(-250.f, 250.f);
            float sd = rng.uniform(0.1f, 10.f);
            float lowerNormal = rng.uniform(-500.f, 500.f);
            float upperNormal = rng.uniform(lowerNormal, 501.f);
            OptionalFloat f = rng.truncNormal(lowerNormal, upperNormal, mean, sd);
            if (f.hasValue)
            {
                REQUIRE(f.value >= lowerNormal);
                REQUIRE(f.value <= upperNormal);
            }

            // truncGammaUpper
            float upperGamma = rng.uniform(0.f, 1000.f);
            float shape = 2.f, scale = gaps::max(rng.uniform(), gaps::epsilon);
            REQUIRE(rng.truncGammaUpper(upperGamma, shape, scale) <= upperGamma);
        }
    }

    SECTION("Test Mean/Variance")
    {
        std::vector<float> uab;
        float uab_A = 2.3f, uab_B = 152.1f;

        std::vector<uint32_t> u32;
        uint32_t u32_A = 2341, u32_B = 1235123;

        std::vector<uint64_t> u64;
        uint64_t u64_A = 14147423647, u64_B = 14148423647;

        std::vector<int> pois;
        unsigned pois_L = 1230;

        std::vector<float> exp;
        float exp_L = 53.f;

        for (unsigned i = 0; i < 50000; ++i)
        {
            uab.push_back(rng.uniform(uab_A, uab_B));
            u32.push_back(rng.uniform32(u32_A, u32_B));
            u64.push_back(rng.uniform64(u64_A, u64_B));
            pois.push_back(rng.poisson(pois_L));
            exp.push_back(rng.exponential(exp_L));
        }

        testMeanVar("Uniform[a,b]", uab, (uab_A + uab_B) / 2.f, TEST_SQ(uab_B - uab_A) / 12.f);
        testMeanVar("Uniform32[a,b]", u32, (u32_A + u32_B) / 2.f, (TEST_SQ(u32_B - u32_A + 1) - 1) / 12.f);
        testMeanVar("Uniform64[a,b]", u64, (u64_A + u64_B) / 2.f, (TEST_SQ(u64_B - u64_A + 1) - 1) / 12.f);
        testMeanVar("Poisson", pois, pois_L, pois_L);
        testMeanVar("Exponential", exp, 1.f / exp_L, 1.f / TEST_SQ(exp_L));
    }
}

TEST_CASE("Test Random.h - RNG Quality")
{

}