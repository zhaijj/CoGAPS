#ifndef __COGAPS_ALGORITHMS_H__
#define __COGAPS_ALGORITHMS_H__

#include "../data_structures/Matrix.h"
#include "Math.h"

#include <cmath>

#define GAPS_SQ(x) ((x) * (x))

struct AlphaParameters
{
    float s;
    float su;
    
    AlphaParameters(float inS, float inSU)
        : s(inS), su(inSU)
    {}

    AlphaParameters operator+(const AlphaParameters &other) const
    {
        float rs = s + other.s;
        float rsu = su - other.su; // weird
        return AlphaParameters(rs, rsu);
    }
};

namespace gaps
{
namespace algo
{
    bool isVectorZero(const float *vec, unsigned size);
    
    // vector algorithms    
    unsigned whichMin(const Vector &vec);
    float sum(const Vector &vec);
    float min(const Vector &vec);
    float max(const Vector &vec);
    float dot(const Vector &A, const Vector &B);
    Vector rank(Vector vec);
    Vector elementSq(Vector vec);

    ColMatrix pmax(const ColMatrix &mat, float factor);
    
    // generic matrix algorithms
    float sum(const ColMatrix &mat);
    float mean(const ColMatrix &mat);
    float nonZeroMean(const ColMatrix &mat);
    ColMatrix computeStdDev(const ColMatrix &stdMat,
        const ColMatrix &meanMat, unsigned nUpdates);

    // specific matrix algorithms
    ColMatrix matrixMultiplication(const ColMatrix &A, const ColMatrix &BT);

    void copyTranspose(ColMatrix *dest, const ColMatrix &src);

    // chiSq / 2
    float loglikelihood(const ColMatrix &D, const ColMatrix &S,
        const ColMatrix &AP);

    AlphaParameters alphaParameters(unsigned size, const float *D,
        const float *S, const float *AP, const float *mat);

    AlphaParameters alphaParameters(unsigned size, const float *D,
        const float *S, const float *AP, const float *mat1, const float *mat2);

    float deltaLL(unsigned size, const float *D, const float *S,
        const float *AP, const float *mat, float delta);

    float deltaLL(unsigned size, const float *D, const float *S,
        const float *AP, const float *mat1, float delta1, const float *mat2,
        float delta2);

} // namespace algo
} // namespace gaps

#endif
