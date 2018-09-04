#ifndef __COGAPS_PROPOSAL_LOCKS_H__
#define __COGAPS_PROPOSAL_LOCKS_H__

#include "../AtomicDomain.h"
#include "../math/Random.h"
#include "../data_structures/HashSets.h"

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

    AtomicProposal createProposal(uint64_t seed, uint64_t id);

    void rejectBirth();
    void acceptBirth();
    void rejectDeath();
    void acceptDeath();

    void reset();

private:

    uint64_t mLastProcessed;

    double mDomainLength;
    double mNumBins;
    double mAlpha;

    unsigned mMinAtoms;
    unsigned mMaxAtoms;

    float deathProb(double nAtoms) const;

    void possibleBirth();
    void possibleDeath();
};

class ProposalLocationLock
{
public:

    ProposalLocationLock(unsigned nrow, unsigned npat);

    void fillProposal(AtomicProposal *prop, AtomicDomain *domain, unsigned id);

    void waitOnIndex(uint64_t pos);
    void release(const AtomicProposal &prop);
    void reset();

private:
    
    bool sameBin(uint64_t p1, uint64_t p2);

    FixedHashSetU32 mUsedIndices;
    HashSetU64 mDeathPositions;
    HashSetU64 mMovePositions;

    uint64_t mLastProcessed;

    uint64_t mBinSize;
    uint64_t mNumPatterns;
    uint64_t mDomainLength;
    uint64_t mRowLength;
};

#endif // __COGAPS_PROPOSAL_LOCKS_H__