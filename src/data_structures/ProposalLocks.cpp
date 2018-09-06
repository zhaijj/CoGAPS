#include "ProposalLocks.h"

////////////////////////////// AtomicProposal //////////////////////////////////

AtomicProposal::AtomicProposal(uint64_t seed, uint64_t id)
    : rng(seed, id), pos(0), atom1(NULL), atom2(NULL), type(0)
{}

///////////////////////////// ProposalTypeLock /////////////////////////////////

ProposalTypeLock::ProposalTypeLock(unsigned nrow, unsigned npat)
    :
mMinAtoms(0), mMaxAtoms(0), mAlpha(0.0), mLastProcessed(0)
{
    uint64_t binSize = std::numeric_limits<uint64_t>::max() / 
        static_cast<uint64_t>(nrow * npat);
    mDomainLength = static_cast<double>(binSize * nrow * npat);
    mNumBins = static_cast<double>(nrow * npat);        
}

void ProposalTypeLock::setAlpha(double alpha)
{
    mAlpha = alpha;
}

AtomicProposal ProposalTypeLock::createProposal(uint64_t seed, uint64_t id)
{   
    AtomicProposal prop(seed, id);
    float u1 = prop.rng.uniform();
    float u2 = prop.rng.uniform();

    waitForTurn(id);

    // special indeterminate case
    while (mMinAtoms < 2 && mMaxAtoms >= 2)
    {
        #pragma omp taskyield
    }

    if (mMaxAtoms < 2) // always birth with 0 or 1 atom
    {
        prop.type = 'B';
        possibleBirth();
    }
    else if (u1 < 0.5f) // decide between birth/death
    {
        // wait until we can determine birth/death
        float lowerBound = deathProb(static_cast<double>(mMinAtoms));
        float upperBound = deathProb(static_cast<double>(mMaxAtoms));
        while (lowerBound < u2 && u2 < upperBound)
        {
            #pragma omp taskyield
            lowerBound = deathProb(static_cast<double>(mMinAtoms));
            upperBound = deathProb(static_cast<double>(mMaxAtoms));
        }

        if (u2 < lowerBound)
        {
            prop.type = 'D';
            possibleDeath();
        }
        else
        {
            prop.type = 'B';
            possibleBirth();
        }
    }
    else // decide between move/exchange
    {
        prop.type = (u1 < 0.75f) ? 'M' : 'E';
    }

    finishTurn();
    return prop;
}

void ProposalTypeLock::reset()
{
    GAPS_ASSERT(mMinAtoms == mMaxAtoms);

    mLastProcessed = 0;
}

void ProposalTypeLock::rejectBirth()
{
    #pragma omp atomic
    --mMaxAtoms;
}

void ProposalTypeLock::acceptBirth()
{
    #pragma omp atomic
    ++mMinAtoms;
}

void ProposalTypeLock::rejectDeath()
{
    #pragma omp atomic
    ++mMinAtoms;
}

void ProposalTypeLock::acceptDeath()
{
    #pragma omp atomic
    --mMaxAtoms;
}

void ProposalTypeLock::waitForTurn(unsigned id)
{
    while (mLastProcessed < id - 1)
    {
        #pragma omp taskyield
    }
}

void ProposalTypeLock::finishTurn()
{
    #pragma omp atomic
    ++mLastProcessed;
}

float ProposalTypeLock::deathProb(double nAtoms) const
{
    double numer = nAtoms * mDomainLength;
    return numer / (numer + mAlpha * mNumBins * (mDomainLength - nAtoms));
}

void ProposalTypeLock::possibleBirth()
{
    ++mMaxAtoms;
}

void ProposalTypeLock::possibleDeath()
{
    --mMinAtoms;
}

/////////////////////////// ProposalLocationLock ///////////////////////////////

ProposalLocationLock::ProposalLocationLock(unsigned nrow, unsigned npat,
unsigned nThreads)
    :
mUsedMatrixIndices(nrow), mPotentialDeaths(nThreads),
mPotentialMoves(nThreads), mDeathCanaries(nThreads)
{
    mBinSize = std::numeric_limits<uint64_t>::max() / static_cast<uint64_t>(nrow * npat);
    mNumPatterns = npat;
    mDomainLength = mBinSize * static_cast<uint64_t>(nrow * npat);
    mRowLength = mBinSize * mNumPatterns;
}

void ProposalLocationLock::fillProposal(AtomicProposal *prop,
AtomicDomain *domain, unsigned id)
{
    waitForTurn(id); // proposals are created in order
    switch (prop->type)
    {
        case 'B':
            fillBirth(prop, domain);
            break;
        case 'D':
            fillDeath(prop, domain);
            break;
        case 'M':
            fillMove(prop, domain);
            break;
        case 'E':
            fillExchange(prop, domain);
            break;
    }
    finishTurn();
}

void ProposalLocationLock::waitForTurn(unsigned id)
{
    while (mLastProcessed < id - 1)
    {
        #pragma omp taskyield
    }
}

void ProposalLocationLock::finishTurn()
{
    #pragma omp atomic
    ++mLastProcessed;
}

void ProposalLocationLock::rejectDeath(uint64_t pos)
{
    GAPS_ASSERT(mPotentialDeaths.contains(pos));
    mPotentialDeaths.erase(pos);
}

void ProposalLocationLock::acceptDeath(uint64_t pos)
{
    GAPS_ASSERT(mPotentialDeaths.contains(pos));
    mPotentialDeaths.erase(pos);
    mDeathCanaries.signal(pos);
}

void ProposalLocationLock::finishMove(uint64_t pos)
{
    mPotentialMoves.erase(pos);
}

void ProposalLocationLock::releaseIndex(unsigned index)
{
    mUsedMatrixIndices.erase(index);
}

void ProposalLocationLock::reset()
{
    GAPS_ASSERT(mUsedMatrixIndices.isEmpty());
    GAPS_ASSERT(mPotentialMoves.isEmpty());
    GAPS_ASSERT(mPotentialDeaths.isEmpty());
    GAPS_ASSERT(mDeathCanaries.isEmpty());

    mLastProcessed = 0;
}

// data mutex on matrix rows/cols
void ProposalLocationLock::waitOnMatrixIndex(uint64_t pos)
{
    unsigned index = pos / mRowLength;
    while (mUsedMatrixIndices.contains(index))
    {
        #pragma omp taskyield
    }
    mUsedMatrixIndices.insert(index);
}

// data mutex on matrix rows/cols
void ProposalLocationLock::waitOnMatrixIndex(uint64_t p1, uint64_t p2)
{
    unsigned ndx1 = p1 / mRowLength;
    unsigned ndx2 = p2 / mRowLength;
    while (mUsedMatrixIndices.contains(ndx1) || mUsedMatrixIndices.contains(ndx2))
    {
        #pragma omp taskyield
    }
    mUsedMatrixIndices.insert(ndx1);
    mUsedMatrixIndices.insert(ndx2);
}

void ProposalLocationLock::waitUntilMoveResolved(uint64_t pos)
{
    while (mPotentialMoves.contains(pos))
    {
        #pragma omp taskyield
    }
}

void ProposalLocationLock::watchForDeath(uint64_t pos)
{
    mDeathCanaries.listen(pos);
}

bool ProposalLocationLock::hasDied(uint64_t pos)
{
    while (mPotentialDeaths.contains(pos))
    {
        #pragma omp taskyield
    }
    return mDeathCanaries.pop(pos);
}

void ProposalLocationLock::ignoreDeath(uint64_t pos)
{
    mDeathCanaries.pop(pos);
}

void ProposalLocationLock::fillBirth(AtomicProposal *prop, AtomicDomain *domain)
{
    // add an empty atom to a random position
    uint64_t pos = 0;
    #pragma omp critical(AtomicInsertOrErase)
    {
        pos = domain->randomFreePosition(&(prop->rng));
    }
    prop->atom1 = domain->insert(pos, 0.f);
    mPotentialDeaths.insert(prop->atom1->pos);

    // make sure the needed row/col of the matrix are available
    waitOnMatrixIndex(prop->atom1->pos);
}

void ProposalLocationLock::fillDeath(AtomicProposal *prop, AtomicDomain *domain)
{
    uint64_t pos = 0;
    do
    {
        #pragma omp critical(AtomicInsertOrErase)
        {
            prop->atom1 = domain->randomAtom(&(prop->rng));
            pos = prop->atom1->pos;
            watchForDeath(pos);
        }
        #pragma omp taskyield
    } while (hasDied(pos));
    mPotentialDeaths.insert(prop->atom1->pos);

    // make sure the needed row/col of the matrix are available
    waitOnMatrixIndex(prop->atom1->pos);
}

void ProposalLocationLock::fillMove(AtomicProposal *prop, AtomicDomain *domain)
{
    AtomNeighborhood hood;
    uint64_t pos1 = 0, pos2 = 0, pos3 = 0;
    do
    {
        ignoreDeath(pos2);
        ignoreDeath(pos3);
    
        #pragma omp critical(AtomicInsertOrErase)
        {
            hood = domain->randomAtomWithNeighbors(&(prop->rng));
            pos1 = hood.center->pos;
            pos2 = hood.hasLeft() ? hood.left->pos : 0;
            pos3 = hood.hasRight() ? hood.right->pos : mDomainLength;
            watchForDeath(pos1);
            watchForDeath(pos2);
            watchForDeath(pos3);
        }
        #pragma omp taskyield
    } while (hasDied(pos1));
    mPotentialMoves.insert(pos1);

    // make sure neighbors haven't been killed off by a previous proposal
    while (hasDied(pos2))
    {
        #pragma omp critical(AtomicInsertOrErase)
        {
            hood.left = domain->getLeftNeighbor(hood.center->pos);
            pos2 = hood.hasLeft() ? hood.left->pos : 0;
        }
        watchForDeath(pos2);
        #pragma omp taskyield
    }
    while (hasDied(pos3))
    {
        #pragma omp critical(AtomicInsertOrErase)
        {
            hood.right = domain->getRightNeighbor(hood.center->pos);
            pos3 = hood.hasRight() ? hood.right->pos : mDomainLength;
        }
        watchForDeath(pos3);
        #pragma omp taskyield
    }

    // wait for neighbors to be in final position before selecting location
    waitUntilMoveResolved(pos2);
    waitUntilMoveResolved(pos3);
    pos2 = hood.hasLeft() ? hood.left->pos : 0;
    pos3 = hood.hasRight() ? hood.right->pos : mDomainLength;
    prop->pos = prop->rng.uniform64(pos2 + 1, pos3 - 1);
    prop->atom1 = hood.center;

    // make sure the needed row/col of the matrix are available
    waitOnMatrixIndex(prop->atom1->pos, prop->pos);
}

void ProposalLocationLock::fillExchange(AtomicProposal *prop, AtomicDomain *domain)
{
    uint64_t pos1 = 0, pos2 = 0;
    do
    {
        ignoreDeath(pos2);

        #pragma omp critical(AtomicInsertOrErase)
        {
            AtomNeighborhood hood = domain->randomAtomWithRightNeighbor(&(prop->rng));
            prop->atom1 = hood.center;
            prop->atom2 = hood.hasRight() ? hood.right : domain->front();
            pos1 = prop->atom1->pos;
            pos2 = prop->atom2->pos;
            watchForDeath(pos1);
            watchForDeath(pos2);
        }
        #pragma omp taskyield
    } while (hasDied(pos1));

    // make sure right neighbor hasn't been killed off by a previous proposal
    while (hasDied(pos2))
    {
        #pragma omp critical(AtomicInsertOrErase)
        {
            prop->atom2 = domain->getRightNeighbor(prop->atom1->pos);
            if (prop->atom2 == NULL)
            {
                prop->atom2 = domain->front();
            }
            pos2 = prop->atom2->pos;
        }
        watchForDeath(pos2);
        #pragma omp taskyield
    }

    // make sure the needed row/col of the matrix are available
    waitOnMatrixIndex(prop->atom1->pos, prop->atom2->pos);
}

