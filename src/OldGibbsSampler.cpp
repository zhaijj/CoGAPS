#include "GibbsSampler.h"
#include "math/SIMD.h"

void AmplitudeGibbsSampler::updateAPMatrix(unsigned row, unsigned col, float delta)
{
    const float *other = mOtherMatrix->rowPtr(col);
    float *ap = mAPMatrix.rowPtr(row);
    unsigned size = mAPMatrix.nCol();

    gaps::simd::packedFloat pOther, pAP;
    gaps::simd::Index i(0);
    gaps::simd::packedFloat pDelta(delta);
    for (; i <= size - i.increment(); ++i)
    {
        pOther.load(other + i);
        pAP.load(ap + i);
        pAP += pDelta * pOther;
        pAP.store(ap + i);
    }

    for (unsigned j = i.value(); j < size; ++j)
    {
        ap[j] += delta * other[j];
    }
}

void PatternGibbsSampler::updateAPMatrix(unsigned row, unsigned col, float delta)
{
    const float *other = mOtherMatrix->colPtr(row);
    float *ap = mAPMatrix.colPtr(col);
    unsigned size = mAPMatrix.nRow();

    gaps::simd::packedFloat pOther, pAP;
    gaps::simd::Index i(0);
    gaps::simd::packedFloat pDelta(delta);
    for (; i <= size - i.increment(); ++i)
    {
        pOther.load(other + i);
        pAP.load(ap + i);
        pAP += pDelta * pOther;
        pAP.store(ap + i);
    }

    for (unsigned j = i.value(); j < size; ++j)
    {
        ap[j] += delta * other[j];
    }
}

AlphaParameters AmplitudeGibbsSampler::alphaParameters(unsigned r1, unsigned c1,
unsigned r2, unsigned c2)
{
    if (r1 == r2)
    {
        return gaps::algo::alphaParameters(mDMatrix.nCol(), mDMatrix.rowPtr(r1),
            mSMatrix.rowPtr(r1), mAPMatrix.rowPtr(r1), mOtherMatrix->rowPtr(c1),
            mOtherMatrix->rowPtr(c2));
    }
    return alphaParameters(r1, c1) + alphaParameters(r2, c2);
}

AlphaParameters PatternGibbsSampler::alphaParameters(unsigned r1, unsigned c1,
unsigned r2, unsigned c2)
{
    if (c1 == c2)
    {
        return gaps::algo::alphaParameters(mDMatrix.nRow(), mDMatrix.colPtr(c1),
            mSMatrix.colPtr(c1), mAPMatrix.colPtr(c1), mOtherMatrix->colPtr(r1),
            mOtherMatrix->colPtr(r2));
    }
    return alphaParameters(r1, c1) + alphaParameters(r2, c2);
}

float AmplitudeGibbsSampler::computeDeltaLL(unsigned r1, unsigned c1, float m1,
unsigned r2, unsigned c2, float m2)
{
    if (r1 == r2)
    {
        return gaps::algo::deltaLL(mDMatrix.nCol(), mDMatrix.rowPtr(r1),
            mSMatrix.rowPtr(r1), mAPMatrix.rowPtr(r1), mOtherMatrix->rowPtr(c1),
            m1, mOtherMatrix->rowPtr(c2), m2);
    }
    return computeDeltaLL(r1, c1, m1) + computeDeltaLL(r2, c2, m2);
}

float PatternGibbsSampler::computeDeltaLL(unsigned r1, unsigned c1, float m1,
unsigned r2, unsigned c2, float m2)
{
    if (c1 == c2)
    {
        return gaps::algo::deltaLL(mDMatrix.nRow(), mDMatrix.colPtr(c1),
            mSMatrix.colPtr(c1), mAPMatrix.colPtr(c1), mOtherMatrix->colPtr(r1),
            m1, mOtherMatrix->colPtr(r2), m2);
    }
    return computeDeltaLL(r1, c1, m1) + computeDeltaLL(r2, c2, m2);
}
