#ifndef __COGAPS_PROPOSAL_LOCKS_H__
#define __COGAPS_PROPOSAL_LOCKS_H__

#include "../AtomicDomain.h"
#include "../math/Random.h"

#include <vector>
#include <stdint.h>

struct AtomicProposal
{
    GapsRng rng; // used for consistency no matter number of threads 
 
    uint64_t pos; // used for birth/move
    Atom *atom1; // used for death/move/exchange
    Atom *atom2; // used for exchange

    char type; // birth (B), death (D), move (M), exchange (E)

    AtomicProposal(uint64_t seed, uint64_t id);
};

class ProposalTypeLock
{
public:

    ProposalTypeLock(unsigned nrow, unsigned npat);

    void setAlpha(double alpha);
    void setNumAtoms(unsigned n);

    AtomicProposal createProposal(uint64_t seed, uint64_t id);

    void rejectBirth();
    void acceptBirth();
    void rejectDeath();
    void acceptDeath();

    #ifdef GAPS_DEBUG
    bool consistent() const;
    #endif

private:

    double mDomainLength;
    double mNumBins;

    unsigned mMinAtoms;
    unsigned mMaxAtoms;
    double mAlpha;

    float deathProb(double nAtoms) const;

    void possibleBirth();
    void possibleDeath();
};

class ProposalLocationLock
{
public:

    ProposalLocationLock(unsigned nrow, unsigned npat);

    void fillProposal(AtomicProposal *prop, AtomicDomain *domain, unsigned id);

private:
    
    bool sameBin(uint64_t p1, uint64_t p2);

    uint64_t mBinSize;
    uint64_t mNumPatterns;
    uint64_t mDomainLength;
};

class IntFixedHashSet
{
public:

    IntFixedHashSet() : mCurrentKey(1) {}

    void setDimensionSize(unsigned size) {mSet.resize(size, 0);}
    void clear() {++mCurrentKey;}
    bool contains(unsigned n) {return mSet[n] == mCurrentKey;}
    void insert(unsigned n) {mSet[n] = mCurrentKey;}

private:

    std::vector<uint64_t> mSet;
    uint64_t mCurrentKey;
};

// TODO have sorted vector with at least some % of holes
// even distribute entries along it
// when shift happens, should be minimal
class IntDenseOrderedSet
{
public:

    void insert(uint64_t p) {mVec.push_back(p);}
    void clear() {mVec.clear();}

    // inclusive of a and b, TODO improve performance
    bool isEmptyInterval(uint64_t a, uint64_t b)
    {
        for (unsigned i = 0; i < mVec.size(); ++i)
        {
            if (mVec[i] >= a && mVec[i] <= b)
            {
                return false;
            }
        }
        return true;
    }

private:

    std::vector<uint64_t> mVec;
};

#endif // __COGAPS_PROPOSAL_LOCKS_H__