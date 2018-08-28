#ifndef __COGAPS_SIMD_H__
#define __COGAPS_SIMD_H__

#if defined( __AVX2__ ) || defined( __AVX__ )

    #define __GAPS_AVX__
    #include <immintrin.h>
    typedef __m256 gaps_packed_t;

#elif defined( __SSE4_2__ ) || defined ( __SSE4_1__ )

    #define __GAPS_SSE__
    #include <nmmintrin.h>
    typedef __m128 gaps_packed_t;

#else

    typedef float gaps_packed_t;

#endif

// forward declarations
namespace gaps { namespace simd { class Index; } }
const float* operator+(const float *ptr, gaps::simd::Index ndx);
float* operator+(float *ptr, gaps::simd::Index ndx);

namespace gaps
{
namespace simd
{

class Index
{
public:

    explicit Index(unsigned i);
    Index& operator=(unsigned val);
    bool operator<(unsigned comp) const;
    bool operator<=(unsigned comp) const;
    void operator++();
    unsigned value() const;
    unsigned increment() const;

    friend const float* ::operator+(const float *ptr, Index ndx);
    friend float* ::operator+(float *ptr, Index ndx);

private:

    unsigned index;
};

class PackedFloat
{
public:

    PackedFloat();
    explicit PackedFloat(float val);

#if defined( __GAPS_SSE__ ) || defined( __GAPS_AVX__ )
    explicit PackedFloat(gaps_packed_t val) : mData(val) {}
#endif

    PackedFloat operator+(PackedFloat b) const;
    PackedFloat operator-(PackedFloat b) const;
    PackedFloat operator*(PackedFloat b) const;
    PackedFloat operator/(PackedFloat b) const;

    void operator+=(PackedFloat val);
    void load(const float *ptr);
    void store(float *ptr);

    float scalar();

private:

    gaps_packed_t mData;
};

} // namespace simd
} // namespace gaps

#endif // __COGAPS_SIMD_H__

