#ifndef __COGAPS_GIBBS_SAMPLER_H__
#define __COGAPS_GIBBS_SAMPLER_H__

#include "Archive.h"
#include "AtomicDomain.h"
#include "GapsAssert.h"
#include "ProposalQueue.h"
#include "data_structures/Matrix.h"
#include "math/Algorithms.h"
#include "math/Random.h"

#include <algorithm>

// forward declarations needed for friend classes/functions

class AmplitudeGibbsSampler;
class PatternGibbsSampler;
class GapsStatistics;

template <class T, class MatA, class MatB>
class GibbsSampler;

template <class T, class MatA, class MatB>
inline Archive& operator<<(Archive &ar, GibbsSampler<T, MatA, MatB> &samp);

template <class T, class MatA, class MatB>
inline Archive& operator>>(Archive &ar, GibbsSampler<T, MatA, MatB> &samp);

/*************************** GIBBS SAMPLER INTERFACE **************************/

template <class T, class MatA, class MatB>
class GibbsSampler
{
public:

    template <class DataType>
    GibbsSampler(const DataType &data, bool transpose, unsigned nPatterns,
        bool amp, bool partitionRows, const std::vector<unsigned> &indices);

    template <class DataType>
    void setUncertainty(const DataType &unc, bool transposeData,
        bool partitionRows, const std::vector<unsigned> &indices);
    
    void setSparsity(float alpha, bool singleCell);
    void setMaxGibbsMass(float max);
    void setAnnealingTemp(float temp);

    void setMatrix(const Matrix &mat);

    void update(unsigned nSteps, unsigned nCores);

    unsigned dataRows() const;
    unsigned dataCols() const;

    float chi2() const;
    uint64_t nAtoms() const;

    #ifdef GAPS_DEBUG
    float getAvgQueue() const;
    bool internallyConsistent();
    #endif

protected:

    MatB mDMatrix;
    MatB mSMatrix;
    MatB mAPMatrix;

    MatA mMatrix;
    MatB* mOtherMatrix;

    ProposalQueue mQueue;
    AtomicDomain mDomain;

    float mLambda;
    float mMaxGibbsMass;
    float mAnnealingTemp;
    
    unsigned mNumRows;
    unsigned mNumCols;
    uint64_t mBinSize;

    float mAvgQueue;
    float mNumQueues;

    T* impl();

    void processProposal(const AtomicProposal &prop);

    void addMass(uint64_t pos, float mass, unsigned row, unsigned col);
    void removeMass(uint64_t pos, float mass, unsigned row, unsigned col);

    void birth(const AtomicProposal &proposal);
    void death(const AtomicProposal &proposal);
    void move(const AtomicProposal &proposal);
    void exchange(const AtomicProposal &proposal);
    void exchangeMH(const AtomicProposal &proposal, unsigned r1, unsigned c1,
        unsigned r2, unsigned c2);

    bool updateAtomMass(uint64_t pos, float mass, float delta);
    void acceptExchange(const AtomicProposal &proposal, float delta1,
        unsigned r1, unsigned c1, unsigned r2, unsigned c2);

private:

    friend class GapsStatistics;

    // serialization
    friend Archive& operator<< <T, MatA, MatB> (Archive &ar, GibbsSampler &sampler);
    friend Archive& operator>> <T, MatA, MatB> (Archive &ar, GibbsSampler &sampler);
};

class AmplitudeGibbsSampler : public GibbsSampler<AmplitudeGibbsSampler, ColMatrix, RowMatrix>
{
private:

    friend class GibbsSampler;
    friend class PatternGibbsSampler;

    unsigned getRow(uint64_t pos) const;
    unsigned getCol(uint64_t pos) const;
    bool canUseGibbs(unsigned row, unsigned col) const;
    bool canUseGibbs(unsigned r1, unsigned c1, unsigned r2, unsigned c2) const;
    void updateAPMatrix(unsigned row, unsigned col, float delta);

    AlphaParameters alphaParameters(unsigned row, unsigned col);
    AlphaParameters alphaParameters(unsigned r1, unsigned c1, unsigned r2,
        unsigned c2);

    float computeDeltaLL(unsigned row, unsigned col, float mass);
    float computeDeltaLL(unsigned r1, unsigned c1, float m1, unsigned r2,
        unsigned c2, float m2);

public:

    template <class DataType>
    AmplitudeGibbsSampler(const DataType &data, bool transposeData,
        unsigned nPatterns, bool partitionRows,
        const std::vector<unsigned> &indices);

    void sync(PatternGibbsSampler &sampler);
    void recalculateAPMatrix();
};

class PatternGibbsSampler : public GibbsSampler<PatternGibbsSampler, RowMatrix, ColMatrix>
{
private:

    friend class GibbsSampler;
    friend class AmplitudeGibbsSampler;

    unsigned getRow(uint64_t pos) const;
    unsigned getCol(uint64_t pos) const;
    bool canUseGibbs(unsigned row, unsigned col) const;
    bool canUseGibbs(unsigned r1, unsigned c1, unsigned r2, unsigned c2) const;
    void updateAPMatrix(unsigned row, unsigned col, float delta);

    AlphaParameters alphaParameters(unsigned row, unsigned col);
    AlphaParameters alphaParameters(unsigned r1, unsigned c1, unsigned r2,
        unsigned c2);

    float computeDeltaLL(unsigned row, unsigned col, float mass);
    float computeDeltaLL(unsigned r1, unsigned c1, float m1, unsigned r2,
        unsigned c2, float m2);

public:

    template <class DataType>
    PatternGibbsSampler(const DataType &data, bool transposeData,
        unsigned nPatterns, bool partitionRows,
        const std::vector<unsigned> &indices);

    void sync(AmplitudeGibbsSampler &sampler);
    void recalculateAPMatrix();
};

/******************** IMPLEMENTATION OF TEMPLATED FUNCTIONS *******************/

template <class DataType>
AmplitudeGibbsSampler::AmplitudeGibbsSampler(const DataType &data,
bool transposeData, unsigned nPatterns, bool partitionRows,
const std::vector<unsigned> &indices)
    :
GibbsSampler(data, transposeData, nPatterns, true, partitionRows, indices)
{
    mQueue.setDimensionSize(mBinSize, mNumCols);
}

template <class DataType>
PatternGibbsSampler::PatternGibbsSampler(const DataType &data,
bool transposeData, unsigned nPatterns, bool partitionRows,
const std::vector<unsigned> &indices)
    :
GibbsSampler(data, transposeData, nPatterns, false, partitionRows, indices)
{
    mQueue.setDimensionSize(mBinSize, mNumRows);
}

template <class T, class MatA, class MatB>
template <class DataType>
GibbsSampler<T, MatA, MatB>::GibbsSampler(const DataType &data,
bool transposeData, unsigned nPatterns, bool amp, bool partitionRows,
const std::vector<unsigned> &indices)
    :
mDMatrix(data, transposeData, partitionRows, indices),
mSMatrix(mDMatrix.pmax(0.1f, 0.1f)),
mAPMatrix(mDMatrix.nRow(), mDMatrix.nCol()),
mMatrix(amp ? mDMatrix.nRow() : nPatterns, amp ? nPatterns : mDMatrix.nCol()),
mOtherMatrix(NULL), mQueue(mMatrix.nRow() * mMatrix.nCol()), mLambda(0.f),
mMaxGibbsMass(100.f), mAnnealingTemp(1.f), mNumRows(mMatrix.nRow()),
mNumCols(mMatrix.nCol()), mAvgQueue(0.f), mNumQueues(0.f)
{
    // calculate atomic domain size
    mBinSize = std::numeric_limits<uint64_t>::max()
        / static_cast<uint64_t>(mNumRows * mNumCols);
    mQueue.setDomainSize(mBinSize * mNumRows * mNumCols);
    mDomain.setDomainSize(mBinSize * mNumRows * mNumCols);

    // default sparsity parameters
    setSparsity(0.01, false);
}

template <class T, class MatA, class MatB>
template <class DataType>
void GibbsSampler<T, MatA, MatB>::setUncertainty(const DataType &unc,
bool transpose, bool partitionRows, const std::vector<unsigned> &indices)
{
    mSMatrix = MatB(unc, transpose, partitionRows, indices);
}

template <class T, class MatA, class MatB>
void GibbsSampler<T, MatA, MatB>::setSparsity(float alpha, bool singleCell)
{
    mQueue.setAlpha(alpha);

    float meanD = singleCell ? gaps::algo::nonZeroMean(mDMatrix) :
        gaps::algo::mean(mDMatrix);

    unsigned nPatterns = mDMatrix.nRow() == mMatrix.nRow() ? mMatrix.nCol() :
        mMatrix.nRow();

    mLambda = alpha * std::sqrt(nPatterns / meanD);
}

template <class T, class MatA, class MatB>
void GibbsSampler<T, MatA, MatB>::setMaxGibbsMass(float max)
{
    mMaxGibbsMass = max;
}

template <class T, class MatA, class MatB>
void GibbsSampler<T, MatA, MatB>::setAnnealingTemp(float temp)
{
    mAnnealingTemp = temp;
}

template <class T, class MatA, class MatB>
void GibbsSampler<T, MatA, MatB>::setMatrix(const Matrix &mat)
{   
    mMatrix = mat;
}

template <class T, class MatA, class MatB>
void GibbsSampler<T, MatA, MatB>::update(unsigned nSteps, unsigned nCores)
{
    unsigned n = 0;
    while (n < nSteps)
    {
        // populate queue, prepare domain for this queue
        mQueue.populate(mDomain, nSteps - n);
        mDomain.resetCache(mQueue.size());
        n += mQueue.size();
        
        // update average queue count
        #ifdef GAPS_DEBUG
        mNumQueues += 1.f;
        mAvgQueue *= (mNumQueues - 1.f) / mNumQueues;
        mAvgQueue += mQueue.size() / mNumQueues;
        #endif

        // process all proposed updates
        #pragma omp parallel for num_threads(nCores)
        for (unsigned i = 0; i < mQueue.size(); ++i)
        {
            processProposal(mQueue[i]);
        }
        mDomain.flushCache();
        mQueue.clear();
    }
}

template <class T, class MatA, class MatB>
unsigned GibbsSampler<T, MatA, MatB>::dataRows() const
{
    return mDMatrix.nRow();
}

template <class T, class MatA, class MatB>
unsigned GibbsSampler<T, MatA, MatB>::dataCols() const
{
    return mDMatrix.nCol();
}

template <class T, class MatA, class MatB>
float GibbsSampler<T, MatA, MatB>::chi2() const
{   
    return 2.f * gaps::algo::loglikelihood(mDMatrix, mSMatrix, mAPMatrix);
}
  
template <class T, class MatA, class MatB>
uint64_t GibbsSampler<T, MatA, MatB>::nAtoms() const
{   
    return mDomain.size();
}

template <class T, class MatA, class MatB>
Archive& operator<<(Archive &ar, GibbsSampler<T, MatA, MatB> &sampler)
{
    ar << sampler.mMatrix << sampler.mAPMatrix << sampler.mQueue <<
        sampler.mDomain << sampler.mLambda << sampler.mMaxGibbsMass <<
        sampler.mAnnealingTemp << sampler.mNumRows << sampler.mNumCols <<
        sampler.mBinSize << sampler.mAvgQueue << sampler.mNumQueues;
    return ar;
}

template <class T, class MatA, class MatB>
Archive& operator>>(Archive &ar, GibbsSampler<T, MatA, MatB> &sampler)
{
    ar >> sampler.mMatrix >> sampler.mAPMatrix >> sampler.mQueue >>
        sampler.mDomain >> sampler.mLambda >> sampler.mMaxGibbsMass >>
        sampler.mAnnealingTemp >> sampler.mNumRows >> sampler.mNumCols >>
        sampler.mBinSize >> sampler.mAvgQueue >> sampler.mNumQueues;
    return ar;
}

static inline OptionalFloat gibbsMass(const AtomicProposal &proposal,
AlphaParameters alpha, float lambdaTerm, float lower, float upper)
{
    if (alpha.s > gaps::epsilon)
    {
        float mean = (alpha.su - lambdaTerm) / alpha.s;
        float sd = 1.f / std::sqrt(alpha.s);
        return proposal.rng.truncNormal(lower, upper, mean, sd);
    }
    return OptionalFloat();
}

template <class T, class MatA, class MatB>
T* GibbsSampler<T, MatA, MatB>::impl()
{
    return static_cast<T*>(this);
}

template <class T, class MatA, class MatB>
void GibbsSampler<T, MatA, MatB>::processProposal(const AtomicProposal &prop)
{
    switch (prop.type)
    {
        case 'B':
            birth(prop);
            break;
        case 'D':
            death(prop);
            break;
        case 'M':
            move(prop);
            break;
        case 'E':
            exchange(prop);
            break;
    }
}

template <class T, class MatA, class MatB>
void GibbsSampler<T, MatA, MatB>::addMass(uint64_t pos, float mass,
unsigned row, unsigned col)
{
    mDomain.cacheInsert(pos, mass);
    mMatrix(row, col) += mass;
    impl()->updateAPMatrix(row, col, mass);
}

template <class T, class MatA, class MatB>
void GibbsSampler<T, MatA, MatB>::removeMass(uint64_t pos, float mass,
unsigned row, unsigned col)
{
    mDomain.cacheErase(pos);
    float newValue = gaps::max(mMatrix(row, col) - mass, 0.f);
    impl()->updateAPMatrix(row, col, newValue - mMatrix(row, col));
    mMatrix(row, col) = newValue;
}

// add an atom at pos, calculate mass either with an exponential distribution
// or with the gibbs mass calculation
template <class T, class MatA, class MatB>
void GibbsSampler<T, MatA, MatB>::birth(const AtomicProposal &proposal)
{
    unsigned row = impl()->getRow(proposal.pos1);
    unsigned col = impl()->getCol(proposal.pos1);

    // calculate proposed mass
    float mass = 0.f;
    if (impl()->canUseGibbs(row, col))
    {
        AlphaParameters alpha = impl()->alphaParameters(row, col);
        mass = gibbsMass(proposal, alpha * mAnnealingTemp, mLambda, 0.f,
            mMaxGibbsMass / mLambda).value; // if it fails this stays 0
    }
    else
    {
        mass = proposal.rng.exponential(mLambda);
    }
    GAPS_ASSERT_MSG(mass < mMaxGibbsMass / mLambda, mass << " " << mLambda);

    // accept mass as long as it's non-zero
    if (mass >= gaps::epsilon)
    {
        addMass(proposal.pos1, mass, row, col);
        mQueue.acceptBirth();
    }
    else
    {
        mQueue.rejectBirth();
    }
}

// automatically accept death, attempt a rebirth at the same position, using
// the original mass or the gibbs mass calculation
template <class T, class MatA, class MatB>
void GibbsSampler<T, MatA, MatB>::death(const AtomicProposal &proposal)
{
    unsigned row = impl()->getRow(proposal.pos1);
    unsigned col = impl()->getCol(proposal.pos1);

    // kill off atom
    float newValue = gaps::max(mMatrix(row, col) - proposal.mass1, 0.f);
    impl()->updateAPMatrix(row, col, newValue - mMatrix(row, col));
    mMatrix(row, col) = newValue;

    // calculate rebirth mass
    float rebirthMass = proposal.mass1;
    AlphaParameters alpha = impl()->alphaParameters(row, col);
    if (impl()->canUseGibbs(row, col))
    {
        OptionalFloat gMass = gibbsMass(proposal, alpha * mAnnealingTemp,
            mLambda, 0.f, mMaxGibbsMass / mLambda);
        if (gMass.hasValue)
        {
            rebirthMass = gMass.value;
            GAPS_ASSERT_MSG(rebirthMass < mMaxGibbsMass / mLambda, rebirthMass);
        }
    }

    // accept/reject rebirth
    float deltaLL = rebirthMass * (alpha.su - alpha.s * rebirthMass / 2.f);
    if (deltaLL * mAnnealingTemp >= std::log(proposal.rng.uniform()))
    {
        mDomain.updateMass(proposal.pos1, rebirthMass);
        mMatrix(row, col) += rebirthMass;
        impl()->updateAPMatrix(row, col, rebirthMass);
        mQueue.rejectDeath();
    }
    else
    {
        mDomain.cacheErase(proposal.pos1);
        mQueue.acceptDeath();
    }
}

template <class T, class MatA, class MatB>
void GibbsSampler<T, MatA, MatB>::move(const AtomicProposal &proposal)
{
    unsigned r1 = impl()->getRow(proposal.pos1);
    unsigned c1 = impl()->getCol(proposal.pos1);
    unsigned r2 = impl()->getRow(proposal.pos2);
    unsigned c2 = impl()->getCol(proposal.pos2);

    if (r1 != r2 || c1 != c2) // automatically reject if change in same bin
    {
        float deltaLL = impl()->computeDeltaLL(r1, c1, -proposal.mass1,
            r2, c2, proposal.mass1);
        if (deltaLL * mAnnealingTemp > std::log(proposal.rng.uniform()))
        {
            removeMass(proposal.pos1, proposal.mass1, r1, c1);
            addMass(proposal.pos2, proposal.mass1, r2, c2);
        }
    }
}

// exchange some amount of mass between two positions, note it is possible
// for one of the atoms to be deleted if it's mass becomes too small
template <class T, class MatA, class MatB>
void GibbsSampler<T, MatA, MatB>::exchange(const AtomicProposal &proposal)
{
    unsigned r1 = impl()->getRow(proposal.pos1);
    unsigned c1 = impl()->getCol(proposal.pos1);
    unsigned r2 = impl()->getRow(proposal.pos2);
    unsigned c2 = impl()->getCol(proposal.pos2);

    if (r1 != r2 || c1 != c2) // automatically reject if change in same bin
    {
        if (impl()->canUseGibbs(r1, c1, r2, c2)) // attempt gibbs change
        {
            AlphaParameters alpha = impl()->alphaParameters(r1, c1, r2, c2);
            OptionalFloat gMass = gibbsMass(proposal, alpha * mAnnealingTemp,
                0.f, -proposal.mass1, proposal.mass2);
            if (gMass.hasValue) // false if too far in tail of distribution to sample
            {
                acceptExchange(proposal, gMass.value - proposal.mass1, r1, c1, r2, c2);
                return;
            }
        }
        // use metropolis-hastings if gibbs fails
        return exchangeMH(proposal, r1, c1, r2, c2);
    }
    mQueue.rejectDeath();
}

// metropolis-hastings version of exchange
template <class T, class MatA, class MatB>
void GibbsSampler<T, MatA, MatB>::exchangeMH(const AtomicProposal &proposal,
unsigned r1, unsigned c1, unsigned r2, unsigned c2)
{
    float totalMass = proposal.mass1 + proposal.mass2;
    GAPS_ASSERT(totalMass > 0.f);
    float newMass = proposal.rng.truncGammaUpper(totalMass, 2.f, 1.f / mLambda);

    float biggerMass = gaps::max(proposal.mass1, proposal.mass2);
    float smallerMass = gaps::min(proposal.mass1, proposal.mass2);

    float delta1 = proposal.mass1 > proposal.mass2 ? newMass - proposal.mass1
        : proposal.mass2 - newMass;
    float oldMass = 2.f * newMass > totalMass ? biggerMass : smallerMass;

    float pNew = gaps::d_gamma(newMass, 2.f, 1.f / mLambda);
    float pOld = gaps::d_gamma(oldMass, 2.f, 1.f / mLambda);

    if (pNew != 0.f || pOld == 0.f)
    {
        if (pNew > 0.f && pOld == 0.f)
        {
            acceptExchange(proposal, delta1, r1, c1, r2, c2);
            return;
        }
        float priorLL = (pOld == 0.f) ? 1.f : pOld / pNew;
        float deltaLL = impl()->computeDeltaLL(r1, c1, delta1, r2, c2, -delta1);
        if (std::log(proposal.rng.uniform() * priorLL) < deltaLL * mAnnealingTemp)
        {
            acceptExchange(proposal, delta1, r1, c1, r2, c2);
            return;
        }
    }
    mQueue.rejectDeath();
}

// helper function for exchange step, updates the atomic domain, matrix, and
// A*P matrix
template <class T, class MatA, class MatB>
void GibbsSampler<T, MatA, MatB>::acceptExchange(const AtomicProposal &proposal,
float delta1, unsigned r1, unsigned c1, unsigned r2, unsigned c2)
{
    float delta2 = -delta1;
    bool death1 = updateAtomMass(proposal.pos1, proposal.mass1, delta1);
    bool death2 = updateAtomMass(proposal.pos2, proposal.mass2, delta2);
    GAPS_ASSERT(!death1 || !death2);
    
    // delete entire atom if resize would make it too small
    if (death1) { delta1 = -proposal.mass1; }
    if (death2) { delta2 = -proposal.mass2; }

    if (!(death1 || death2))
    {
        mQueue.rejectDeath();
    }

    mMatrix(r1, c1) = gaps::max(mMatrix(r1, c1) + delta1, 0.f);
    mMatrix(r2, c2) = gaps::max(mMatrix(r2, c2) + delta2, 0.f);
    impl()->updateAPMatrix(r1, c1, delta1);
    impl()->updateAPMatrix(r2, c2, delta2);
}

// helper function for acceptExchange, used to conditionally update the mass
// at a single position, deleting the atom if the new mass is too small
template <class T, class MatA, class MatB>
bool GibbsSampler<T, MatA, MatB>::updateAtomMass(uint64_t pos, float mass,
float delta)
{
    if (mass + delta < gaps::epsilon)
    {
        mDomain.cacheErase(pos);
        mQueue.acceptDeath();
        return true;
    }
    mDomain.updateMass(pos, mass + delta);
    return false;
}

#ifdef GAPS_DEBUG
template <class T, class MatA, class MatB>
float GibbsSampler<T, MatA, MatB>::getAvgQueue() const
{
    return mAvgQueue;
}

template <class T, class MatA, class MatB>
bool GibbsSampler<T, MatA, MatB>::internallyConsistent()
{
    if (mDomain.size() > 0)
    {
        Atom a = mDomain.front();
        float current = a.mass;
        uint64_t row = impl()->getRow(a.pos);
        uint64_t col = impl()->getCol(a.pos);

        while (mDomain.hasRight(a))
        {
            a = mDomain.right(a);
            if (row != impl()->getRow(a.pos) || col != impl()->getCol(a.pos))
            {
                float matVal = mMatrix(row, col);
                if (std::abs(current - matVal) > 0.1f)
                {
                    gaps_printf("mass difference detected at row %lu, column %lu: %f %f\n",
                        row, col, current, matVal); 
                    return false;
                }
                
                row = impl()->getRow(a.pos);
                col = impl()->getCol(a.pos);
                current = a.mass;
            }
            else
            {
                current += a.mass;
            }
        }
        return true;
    }
    else
    {
        return gaps::algo::sum(mMatrix) == 0.f;
    }
}
#endif // GAPS_DEBUG

#endif // __COGAPS_GIBBS_SAMPLER_H__
