#include "ProposalLocks.h"

////////////////////////////// AtomicProposal //////////////////////////////////

AtomicProposal::AtomicProposal(uint64_t seed, uint64_t id)
    : rng(seed, id), pos(0), atom1(NULL), atom2(NULL), type(0)
{}

///////////////////////////// ProposalTypeLock /////////////////////////////////

ProposalTypeLock::ProposalTypeLock(unsigned nrow, unsigned npat)
    :
mMinAtoms(0), mMaxAtoms(0), mAlpha(0.0)
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

void ProposalTypeLock::setNumAtoms(unsigned n)
{
    mMinAtoms = n;
    mMaxAtoms = n;
}

AtomicProposal ProposalTypeLock::createProposal(uint64_t seed, uint64_t id)
{
    GAPS_ASSERT(mMinAtoms == mMaxAtoms);
    AtomicProposal prop(seed, id);

    // check special condition
    if (mMinAtoms < 2)
    {
        prop.type = 'B';
        possibleBirth();
        return prop;
    }

    // decide between birth/death
    float u1 = prop.rng.uniform();
    if (u1 < 1.f)
    {
        if (prop.rng.uniform() < deathProb(static_cast<double>(mMinAtoms)))
        {
            prop.type = 'D';
            possibleDeath();
        }
        else
        {
            prop.type = 'B';
            possibleBirth();
        }
        return prop;
    }

    // decide between move/exchange
    if (u1 < 0.75f)
    {
        prop.type = 'M';
    }
    else
    {
        prop.type = 'E';
        possibleDeath();
    }
    return prop;
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

#ifdef GAPS_DEBUG
bool ProposalTypeLock::consistent() const
{
    return mMinAtoms == mMaxAtoms;
}
#endif

/////////////////////////// ProposalLocationLock ///////////////////////////////

ProposalLocationLock::ProposalLocationLock(unsigned nrow, unsigned npat)
{
    mBinSize = std::numeric_limits<uint64_t>::max() / 
        static_cast<uint64_t>(nrow * npat);
    mNumPatterns = npat;
    mDomainLength = mBinSize * static_cast<uint64_t>(nrow * npat);
}

void ProposalLocationLock::fillProposal(AtomicProposal *prop,
AtomicDomain *domain, unsigned id)
{
    AtomNeighborhood hood;
    uint64_t lbound, rbound;
    switch(prop->type)
    {
        case 'B': // birth

            prop->pos = domain->randomFreePosition(&(prop->rng));
            break;

        case 'D': // death

            prop->atom1 = domain->randomAtom(&(prop->rng));
            break;

        case 'M': // move

            hood = domain->randomAtomWithNeighbors(&(prop->rng));
            lbound = hood.hasLeft() ? hood.left->pos : 0;
            rbound = hood.hasRight() ? hood.right->pos : mDomainLength;
            prop->pos = prop->rng.uniform64(lbound + 1, rbound - 1);
            prop->atom1 = hood.center;
            break;

        case 'E': // exchange

            hood = domain->randomAtomWithRightNeighbor(&(prop->rng));
            prop->atom2 = hood.hasRight() ? hood.right : domain->front();
            prop->atom1 = hood.center;
            break;
    }
}
