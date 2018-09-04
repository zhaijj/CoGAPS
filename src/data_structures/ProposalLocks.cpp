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

// TODO use boost conditional variables to wait
static void busyWait()
{
    volatile unsigned a;
    for (unsigned i = 0; i < 10; ++i)
    {
        a = i;
    }
}

AtomicProposal ProposalTypeLock::createProposal(uint64_t seed, uint64_t id)
{   
    GAPS_ASSERT(mMinAtoms == mMaxAtoms);
    AtomicProposal prop(seed, id);

    // wait for turn - proposals are created in order
    //while (mLastProcessed != id - 1)
    //{
    //    busyWait();
    //}

    // special indeterminate case
    //while (mMinAtoms < 2 && mMaxAtoms >= 2)
    //{
    //    busyWait();
    //}

    // check special condition
    if (mMaxAtoms < 2)
    {
        prop.type = 'B';
        possibleBirth();
        ++mLastProcessed;
        return prop;
    }

    // decide between birth/death
    float u1 = prop.rng.uniform();
    if (u1 < 0.5f)
    {
        float lowerBound = deathProb(static_cast<double>(mMinAtoms));
        float upperBound = deathProb(static_cast<double>(mMaxAtoms));

        float u2 = prop.rng.uniform();
        //while (u2 > lowerBound && u2 < upperBound)
        //{
        //    busyWait();
        //}

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
        ++mLastProcessed;
        return prop;
    }

    // decide between move/exchange
    prop.type = (u1 < 0.75f) ? 'M' : 'E';
    ++mLastProcessed;
    return prop;
}

void ProposalTypeLock::reset()
{
    mLastProcessed = 0;
    GAPS_ASSERT(mMinAtoms == mMaxAtoms);
}

void ProposalTypeLock::rejectBirth()
{
    --mMaxAtoms;
}

void ProposalTypeLock::acceptBirth()
{
    ++mMinAtoms;
}

void ProposalTypeLock::rejectDeath()
{
    ++mMinAtoms;
}

void ProposalTypeLock::acceptDeath()
{
    --mMaxAtoms;
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

ProposalLocationLock::ProposalLocationLock(unsigned nrow, unsigned npat)
    :
mUsedIndices(nrow)
{
    mBinSize = std::numeric_limits<uint64_t>::max() / 
        static_cast<uint64_t>(nrow * npat);
    mNumPatterns = npat;
    mDomainLength = mBinSize * static_cast<uint64_t>(nrow * npat);
    mRowLength = mBinSize * mNumPatterns;
}

void ProposalLocationLock::reset()
{
    mLastProcessed = 0;
}

// TODO block on positions as well

void ProposalLocationLock::waitOnIndex(uint64_t pos)
{
    unsigned ndx = pos / mRowLength; // row length = bin size * nPatterns
    //while (mUsedIndices.contains(ndx))
    //{
    //    busyWait();
    //}
    mUsedIndices.insert(ndx);
}

void ProposalLocationLock::waitOnPosition(uint64_t pos)
{
    while (mUsedPositions.contains(pos)) { busyWait(); }
}

void ProposalLocationLock::release(const AtomicProposal &prop)
{
    mUsedIndices.erase(prop.atom1->pos / mRowLength);
    switch (prop.type)
    {
        case 'M':
            mUsedIndices.erase(prop.pos / mRowLength);
            break;
        case 'E':
            mUsedIndices.erase(prop.atom2->pos / mRowLength);
            break;
        default:
            break;
    }            
}

bool ProposalLocationLock::hasDied(uint64_t pos)
{
    while (mDeathPositions.contains(pos1))
    {
        busyWait();
    }
    return mDeaths.pop(pos); // gets and erases
}

void ProposalLocationLock::waitUntilMoveResolved
{
    while (mDeathPositions.contains(pos1))
    {
        busyWait();
    }
}

void ProposalLocationLock::fillProposal(AtomicProposal *prop,
AtomicDomain *domain, unsigned id)
{
    waitForTurn(); // proposals are created in order
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

void ProposalLocationLock::fillBirth(AtomicProposal *prop, AtomicDomain *dom)
{
    // add an empty atom to a random position
    pos = domain->randomFreePosition(&(prop->rng));
    GAPS_ASSERT(pos > 0);
    prop->atom1 = domain->insert(pos, 0.f);

    // make sure the needed row/col of the matrix are available
    waitOnIndex(prop->atom1->pos);
}

void ProposalLocationLock::fillDeath(AtomicProposal *prop, AtomicDomain *dom)
{
    uint64_t pos = 0;
    do
    {
        #pragma omp critical(AtomicInsertOrErase)
        {
            prop->atom1 = domain->randomAtom(&(prop->rng));
            pos = prop->atom1->pos;
        }
    } while (hadDied(pos));
    mDeathPositions.insert(prop->atom1->pos);

    // make sure the needed row/col of the matrix are available
    waitOnIndex(prop->atom1->pos);
}

void ProposalLocationLock::fillMove(AtomicProposal *prop, AtomicDomain *dom)
{
    uint64_t pos1 = 0, pos2 = 0, pos3 = 0;
    do
    {
        #pragma omp critical(AtomicInsertOrErase)
        {
            hood = domain->randomAtomWithNeighbors(&(prop->rng));
            prop->atom1 = hood.center;
            pos1 = prop->atom1->pos;
            pos2 = hood.hasLeft() ? hood.left->pos : 0;
            pos3 = hood.hasRight() ? hood.right->pos : mDomainLength;
        }
    } while (hasDied(pos1));
    mMovePositions.insert(prop->atom1->pos);

    // make sure neighbors haven't been killed off by a previous proposal
    while (hasDied(pos2))
    {
        pos2 = domain->getLeftNeighbor(prop->atom1->pos)->pos
    }
    while (hasDied(pos3))
    {
        pos3 = domain->getRightNeighbor(prop->atom1->pos)->pos
    }

    // wait for neighbors to be in final position before selecting location
    waitUntilMoveResolved(pos2);
    waitUntilMoveResolved(pos3);
    prop->pos = prop->rng.uniform64(pos2 + 1, pos3 - 1);

    // make sure the needed row/col of the matrix are available
    waitOnIndex(prop->atom1->pos);
    waitOnIndex(prop->pos);
}

void ProposalLocationLock::fillExchange(AtomicProposal *prop, AtomicDomain *dom)
{
    uint64_t pos1 = 0, pos2 = 0;
    do
    {
        #pragma omp critical(AtomicInsertOrErase)
        {
            hood = domain->randomAtomWithRightNeighbor(&(prop->rng));
            prop->atom1 = hood.center;
            prop->atom2 = hood.hasRight() ? hood.right : domain->front();
            pos1 = prop->atom1->pos;
            pos2 = prop->atom2->pos;
        }
    } while (hasDied(pos1));

    // make sure right neighbor hasn't been killed off by a previous proposal
    while (hasDied(pos2))
    {
        pos2 = domain->getRightNeighbor(prop->atom1->pos)->pos
    }

    // make sure the needed row/col of the matrix are available
    waitOnIndex(prop->atom1->pos);
    waitOnIndex(prop->atom2->pos);
}

