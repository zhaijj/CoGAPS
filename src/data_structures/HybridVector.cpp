#include "HybridVector.h"
#include "../math/Math.h"

HybridVector::HybridVector(unsigned size)
    : mIndexBitFlags(size / 64 + 1, 0), mData(size, 0.f)
{}

HybridVector::HybridVector(const std::vector<float> &v)
    :
mIndexBitFlags(v.size() / 64 + 1, 0),
mData(v.size(), 0.f)
{
    for (unsigned i = 0; i < v.size(); ++i)
    {
        mData[i] = v[i];
        if (v[i] > 0.f)
        {
            mIndexBitFlags[i / 64] ^= (1ull << (i % 64));
        }
    }
}

bool HybridVector::empty() const
{
    for (unsigned i = 0; i < mIndexBitFlags.size(); ++i)
    {
        if (mIndexBitFlags[i] != 0)
        {
            return false;
        }
    }
    return true;
}

unsigned HybridVector::size() const
{
    return mData.size();
}

bool HybridVector::add(unsigned i, float v)
{
    if (mData[i] + v < gaps::epsilon)
    {
        mIndexBitFlags[i / 64] ^= (1ull << (i % 64));
        mData[i] = 0.f;
        return true;
    }
    else
    {
        mIndexBitFlags[i / 64] |= (1ull << (i % 64));
        mData[i] += v;
        return false;
    }
}

float HybridVector::operator[](unsigned i) const
{
    return mData[i];
}

const float* HybridVector::densePtr() const
{
    return &(mData[0]);
}

Archive& operator<<(Archive &ar, HybridVector &vec)
{
    ar << vec.mData.size();
    for (unsigned i = 0; i < vec.mIndexBitFlags.size(); ++i)
    {
        ar << vec.mIndexBitFlags[i];
    }

    for (unsigned i = 0; i < vec.mData.size(); ++i)
    {
        ar << vec.mData[i];
    }
    return ar;
}

Archive& operator>>(Archive &ar, HybridVector &vec)
{
    unsigned sz = 0;
    ar >> sz;
    GAPS_ASSERT(sz == vec.mData.size());

    for (unsigned i = 0; i < vec.mIndexBitFlags.size(); ++i)
    {
        ar >> vec.mIndexBitFlags[i];
    }

    for (unsigned i = 0; i < vec.mData.size(); ++i)
    {
        ar >> vec.mData[i];
    }
    return ar;
}


