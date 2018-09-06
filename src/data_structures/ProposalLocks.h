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

    void reset();

    void rejectBirth();
    void acceptBirth();
    void rejectDeath();
    void acceptDeath();

private:

    uint64_t mLastProcessed;

    double mDomainLength;
    double mNumBins;
    double mAlpha;

    unsigned mMinAtoms;
    unsigned mMaxAtoms;

    void waitForTurn(unsigned id);
    void finishTurn();

    float deathProb(double nAtoms) const;

    void possibleBirth();
    void possibleDeath();
};

class ProposalLocationLock
{
public:

    ProposalLocationLock(unsigned nrow, unsigned npat, unsigned nThreads);

    void fillProposal(AtomicProposal *prop, AtomicDomain *domain, unsigned id);

    void rejectDeath(uint64_t pos);
    void acceptDeath(uint64_t pos);
    void finishMove(uint64_t pos);

    void releaseIndex(unsigned index);
    void reset();

private:

    void waitForTurn(unsigned id);
    void finishTurn();

    void waitOnMatrixIndex(uint64_t pos);
    void waitOnMatrixIndex(uint64_t p1, uint64_t p2);
    void waitUntilMoveResolved(uint64_t pos);
    void watchForDeath(uint64_t pos);
    bool hasDied(uint64_t pos);
    void ignoreDeath(uint64_t pos);

    void fillBirth(AtomicProposal *prop, AtomicDomain *domain);
    void fillDeath(AtomicProposal *prop, AtomicDomain *domain);
    void fillMove(AtomicProposal *prop, AtomicDomain *domain);
    void fillExchange(AtomicProposal *prop, AtomicDomain *domain);

    FixedHashSetU32 mUsedMatrixIndices;
    SmallHashSetU64 mPotentialDeaths;
    SmallHashSetU64 mPotentialMoves;
    CanaryHashSetU64 mDeathCanaries;

    uint64_t mLastProcessed;

    uint64_t mBinSize;
    uint64_t mNumPatterns;
    uint64_t mDomainLength;
    uint64_t mRowLength;
};

#endif // __COGAPS_PROPOSAL_LOCKS_H__