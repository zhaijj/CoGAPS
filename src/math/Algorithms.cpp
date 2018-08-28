#include "Algorithms.h"
#include "../data_structures/Matrix.h"
#include "../utils/GapsAssert.h"
#include "SIMD.h"

#include <algorithm>

float gaps::algo::sum(const ColMatrix &mat)
{
    float sum = 0.f;
    for (unsigned i = 0; i < mat.nRow(); ++i)
    {
        for (unsigned j = 0; j < mat.nCol(); ++j)
        {
            sum += mat(i,j);
        }
    }
    return sum;
}

float gaps::algo::mean(const ColMatrix &mat)
{
    return gaps::algo::sum(mat) / (mat.nRow() * mat.nCol());
}

float gaps::algo::nonZeroMean(const ColMatrix &mat)
{
    float sum = 0.f;
    unsigned nNonZeros = 0;
    for (unsigned i = 0; i < mat.nRow(); ++i)
    {
        for (unsigned j = 0; j < mat.nCol(); ++j)
        {
            if (mat(i,j) != 0.f)
            {
                nNonZeros++;
                sum += mat(i,j);
            }
        }
    }
    return sum / static_cast<float>(nNonZeros);
}

ColMatrix gaps::algo::computeStdDev(const ColMatrix &stdMat,
const ColMatrix &meanMat, unsigned nUpdates)
{
    GAPS_ASSERT(nUpdates > 1);
    ColMatrix retMat(stdMat.nRow(), stdMat.nCol());
    for (unsigned r = 0; r < retMat.nRow(); ++r)
    {
        for (unsigned c = 0; c < retMat.nCol(); ++c)
        {
            float meanTerm = meanMat(r,c) * meanMat(r,c) / static_cast<float>(nUpdates);
            float numer = gaps::max(0.f, stdMat(r,c) - meanTerm);
            retMat(r,c) = std::sqrt(numer / (static_cast<float>(nUpdates) - 1.f));
        }
    }
    return retMat;
}

float gaps::algo::loglikelihood(const ColMatrix &D, const ColMatrix &S,
const ColMatrix &AP)
{
    GAPS_ASSERT_MSG(D.nRow() == S.nRow(), D.nRow() << " " << S.nRow());
    GAPS_ASSERT_MSG(S.nRow() == AP.nRow(), S.nRow() << " " << AP.nRow());
    GAPS_ASSERT_MSG(D.nCol() == S.nCol(), D.nCol() << " " << S.nCol());
    GAPS_ASSERT_MSG(S.nCol() == AP.nCol(), S.nCol() << " " << AP.nCol());

    float chi2 = 0.f;
    for (unsigned i = 0; i < D.nRow(); ++i)
    {
        for (unsigned j = 0; j < D.nCol(); ++j)
        {
            GAPS_ASSERT(AP(i,j) > 0.f);
            chi2 += GAPS_SQ(D(i,j) - AP(i,j)) / GAPS_SQ(S(i,j));
        }
    }
    return chi2 / 2.f;
}

float gaps::algo::sum(const Vector &vec)
{
    float sum = 0.f;
    for (unsigned i = 0; i < vec.size(); ++i)
    {
        sum += vec[i];
    }
    return sum;
}

float gaps::algo::min(const Vector &vec)
{
    float min = vec[0];
    for (unsigned i = 0; i < vec.size(); ++i)
    {
        min = gaps::min(min, vec[i]);
    }
    return min;
}

float gaps::algo::max(const Vector &vec)
{
    float max = vec[0];
    for (unsigned i = 0; i < vec.size(); ++i)
    {
        max = gaps::max(max, vec[i]);
    }
    return max;
}

float gaps::algo::dot(const Vector &A, const Vector &B)
{
    float dotProd = 0.f;
    for (unsigned i = 0; i < A.size(); ++i)
    {
        dotProd += A[i] * B[i];
    }
    return dotProd;
}

unsigned gaps::algo::whichMin(const Vector &vec)
{
    float min = vec[0];
    unsigned minNdx = 0;
    for (unsigned i = 1; i < vec.size(); ++i)
    {
        if (vec[i] < min)
        {
            min = vec[i];
            minNdx = i;
        }
    }
    return minNdx;
}

Vector gaps::algo::rank(Vector vec)
{
    std::vector< std::pair<float, float> > sortVec(vec.size());
    for (unsigned i = 0; i < vec.size(); ++i)
    {
        sortVec[i] = std::pair<float, float>(vec[i], i);
    }
    
    std::sort(sortVec.begin(), sortVec.end());
    Vector ranks(vec.size());
    for (unsigned i = 0; i < vec.size(); ++i)
    {
        ranks[i] = sortVec[i].second;
    }
    return ranks;
}

Vector gaps::algo::elementSq(Vector vec)
{
    for (unsigned i = 0; i < vec.size(); ++i)
    {
        vec[i] *= vec[i];
    }
    return vec;
}

bool gaps::algo::isVectorZero(const float *vec, unsigned size)
{
    for (unsigned i = 0; i < size; ++i)
    {
        if (vec[i] != 0.f)
        {
            return false;
        }
    }
    return true;
}
    
ColMatrix gaps::algo::pmax(const ColMatrix &mat, float factor)
{
    ColMatrix temp(mat.nRow(), mat.nCol());
    for (unsigned i = 0; i < mat.nRow(); ++i)
    {
        for (unsigned j = 0; j < mat.nCol(); ++j)
        {
            temp(i,j) = gaps::max(mat(i,j) * factor, factor);
        }
    }
    return temp;
}

void gaps::algo::copyTranspose(ColMatrix *dest, const ColMatrix &src)
{
    for (unsigned j = 0; j < dest->nCol(); ++j)
    {
        for (unsigned i = 0; i < dest->nRow(); ++i)
        {
            dest->operator()(i,j) = src(j,i); // TODO test which order is better
        }
    }
}

AlphaParameters gaps::algo::alphaParameters(unsigned size, const float *D,
const float *S, const float *AP, const float *mat)
{
    gaps::simd::PackedFloat pMat, pD, pAP, pS;
    gaps::simd::PackedFloat partialS(0.f), partialSU(0.f);
    gaps::simd::Index i(0);
    for (; i <= size - i.increment(); ++i)
    {   
        pMat.load(mat + i);
        pD.load(D + i);
        pAP.load(AP + i);
        pS.load(S + i);
        gaps::simd::PackedFloat ratio(pMat / pS);
        partialS += ratio * ratio;
        partialSU += (ratio * (pD - pAP)) / pS;
    }
    float s = partialS.scalar(), su = partialSU.scalar();
    for (unsigned j = i.value(); j < size; ++j)
    {
        float term = mat[j] / (S[j] * S[j]);
        s += mat[j] * term;
        su += term * (D[j] - AP[j]);
    }
    return AlphaParameters(s,su);
}

//
AlphaParameters gaps::algo::alphaParameters(unsigned size, const float *D,
const float *S, const float *AP, const float *mat1, const float *mat2)
{
    gaps::simd::PackedFloat ratio, pMat1, pMat2, pD, pAP, pS;
    gaps::simd::PackedFloat partialS(0.f), partialSU(0.f);
    gaps::simd::Index i(0);
    for (; i <= size - i.increment(); ++i)
    {   
        pMat1.load(mat1 + i);
        pMat2.load(mat2 + i);
        pD.load(D + i);
        pAP.load(AP + i);
        pS.load(S + i);
        ratio = (pMat1 - pMat2) / pS;
        partialS += ratio * ratio;
        partialSU += ratio * (pD - pAP) / pS;
    }

    float s = partialS.scalar(), su = partialSU.scalar();
    for (unsigned j = i.value(); j < size; ++j)
    {
        float term = (mat1[j] - mat2[j]) / (S[j] * S[j]);
        s += (mat1[j] - mat2[j]) * term;
        su += term * (D[j] - AP[j]);
    }
    return AlphaParameters(s,su);
}

float gaps::algo::deltaLL(unsigned size, const float *D, const float *S,
const float *AP, const float *mat, float delta)
{
    const gaps::simd::PackedFloat pDelta(delta), two(2.f);
    gaps::simd::PackedFloat d, pMat, pD, pAP, pS, partialSum(0.f);
    gaps::simd::Index i(0);
    for (; i <= size - i.increment(); ++i)
    {   
        pMat.load(mat + i);
        pD.load(D + i);
        pAP.load(AP + i);
        pS.load(S + i);
        d = pDelta * pMat;
        partialSum += (d * (two * (pD - pAP) - d)) / (two * pS * pS);
    }
    float delLL = partialSum.scalar();
    for (unsigned j = i.value(); j < size; ++j)
    {
        float fd = delta * mat[j];
        delLL += (fd * (2.f * (D[j] - AP[j]) - fd)) / (2.f * S[j] * S[j]);
    }
    return delLL;
}

float gaps::algo::deltaLL(unsigned size, const float *D, const float *S,
const float *AP, const float *mat1, float delta1, const float *mat2,
float delta2)
{
    const gaps::simd::PackedFloat pDelta1(delta1), pDelta2(delta2), two(2.f);
    gaps::simd::PackedFloat d, pMat1, pMat2, pD, pAP, pS, partialSum(0.f);
    gaps::simd::Index i(0);
    for (; i <= size - i.increment(); ++i)
    {   
        pMat1.load(mat1 + i);
        pMat2.load(mat2 + i);
        pD.load(D + i);
        pAP.load(AP + i);
        pS.load(S + i);
        d = pDelta1 * pMat1 + pDelta2 * pMat2;
        partialSum += (d * (two * (pD - pAP) - d)) / (two * pS * pS);
    }
    float delLL = partialSum.scalar();
    for (unsigned j = i.value(); j < size; ++j)
    {
        float fd = delta1 * mat1[j] + delta2 * mat2[j];
        delLL += (fd * (2.f * (D[j] - AP[j]) - fd)) / (2.f * S[j] * S[j]);
    }
    return delLL;
}

// multiply A by B transpose - super slow, don't call often
ColMatrix gaps::algo::matrixMultiplication(const ColMatrix &A, const ColMatrix &BT)
{
    GAPS_ASSERT_MSG(A.nCol() == BT.nCol(), A.nCol() << " " << BT.nCol());
    ColMatrix temp(A.nRow(), BT.nRow());
    for (unsigned i = 0; i < A.nRow(); ++i)
    {
        for (unsigned j = 0; j < BT.nRow(); ++j)
        {
            float sum = 0.0;
            for (unsigned k = 0; k < A.nCol(); ++k)
            {
                sum += A(i,k) * BT(j,k);
            }
            temp(i,j) = sum;
        }
    }
    return temp;    
}