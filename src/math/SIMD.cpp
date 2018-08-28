#include "SIMD.h"

namespace gs = gaps::simd;

/////////////////////////////// Index //////////////////////////////////////////

#if defined( __GAPS_AVX__ )
    static const unsigned index_increment = 8;
#elif defined( __GAPS_SSE__ )
    static const unsigned index_increment = 4;
#else
    static const unsigned index_increment = 1;
#endif

gs::Index::Index(unsigned i) : index(i) {}

gs::Index& gs::Index::operator=(unsigned val)
{
    index = val;
    return *this;
}

bool gs::Index::operator<(unsigned comp) const
{
    return index < comp;
}

bool gs::Index::operator<=(unsigned comp) const
{
    return index <= comp;
}

void gs::Index::operator++()
{
    index += index_increment;
}

unsigned gs::Index::value() const
{
    return index;
}

unsigned gs::Index::increment() const
{
    return index_increment;
}

const float* operator+(const float *ptr, gs::Index ndx)
{
    return ptr + ndx.index;
}

float* operator+(float *ptr, gs::Index ndx)
{
    return ptr + ndx.index;
}

///////////////////////////// PackedFloat //////////////////////////////////////

#if defined( __GAPS_AVX__ )
    #define SET_SCALAR(x) _mm256_set1_ps(x)
    #define LOAD_PACKED(x) _mm256_load_ps(x)
    #define STORE_PACKED(p,x) _mm256_store_ps(p,x)
    #define ADD_PACKED(a,b) _mm256_add_ps(a,b)
    #define SUB_PACKED(a,b) _mm256_sub_ps(a,b)
    #define MUL_PACKED(a,b) _mm256_mul_ps(a,b)
    #define DIV_PACKED(a,b) _mm256_div_ps(a,b)
#elif defined( __GAPS_SSE__ )
    #define SET_SCALAR(x) _mm_set1_ps(x)
    #define LOAD_PACKED(x) _mm_load_ps(x)
    #define STORE_PACKED(p,x) _mm_store_ps(p,x)
    #define ADD_PACKED(a,b) _mm_add_ps(a,b)
    #define SUB_PACKED(a,b) _mm_sub_ps(a,b)
    #define MUL_PACKED(a,b) _mm_mul_ps(a,b)
    #define DIV_PACKED(a,b) _mm_div_ps(a,b)
#else
    #define SET_SCALAR(x) x
    #define LOAD_PACKED(x) *(x)
    #define STORE_PACKED(p,x) *(p) = (x)
    #define ADD_PACKED(a,b) ((a)+(b))
    #define SUB_PACKED(a,b) ((a)-(b))
    #define MUL_PACKED(a,b) ((a)*(b))
    #define DIV_PACKED(a,b) ((a)/(b))
#endif

gs::PackedFloat::PackedFloat() : mData() {}

gs::PackedFloat::PackedFloat(float val) : mData(SET_SCALAR(val)) {}

gs::PackedFloat gs::PackedFloat::operator+(gs::PackedFloat b) const
{
    return PackedFloat(ADD_PACKED(mData, b.mData));
}

gs::PackedFloat gs::PackedFloat::operator-(gs::PackedFloat b) const
{
    return PackedFloat(SUB_PACKED(mData, b.mData));
}

gs::PackedFloat gs::PackedFloat::operator*(gs::PackedFloat b) const
{
    return PackedFloat(MUL_PACKED(mData, b.mData));
}

gs::PackedFloat gs::PackedFloat::operator/(gs::PackedFloat b) const
{
    return PackedFloat(DIV_PACKED(mData, b.mData));
}

void gs::PackedFloat::operator+=(PackedFloat val)
{
    mData = ADD_PACKED(mData, val.mData);
}

void gs::PackedFloat::load(const float *ptr)
{
    mData = LOAD_PACKED(ptr);
}

void gs::PackedFloat::store(float *ptr)
{
    STORE_PACKED(ptr, mData);
}

float gs::PackedFloat::scalar()
{
#if defined( __GAPS_AVX__ )

    float* ra = reinterpret_cast<float*>(&mData); // NOLINT
    mData = _mm256_hadd_ps(mData, mData);
    mData = _mm256_hadd_ps(mData, mData);
    return ra[0] + ra[4];

#elif defined( __GAPS_SSE__ )

    float* ra = reinterpret_cast<float*>(&mData); // NOLINT
    return ra[0] + ra[1] + ra[2] + ra[3];

#else

    return mData;

#endif
}
