#include "ProposalLocks.h"

////////////////////////////// AtomicProposal //////////////////////////////////

AtomicProposal::AtomicProposal(uint64_t seed, uint64_t id)
    : rng(seed, id), pos(0), atom1(NULL), atom2(NULL), type(0)
{}

///////////////////////////// ProposalTypeLock /////////////////////////////////

ProposalTypeLock::ProposalTypeLock()
    : mMinAtoms(0), mMaxAtoms(0)
{}

float ProposalTypeLock::deathProb(double nAtoms) const
{
    double numer = nAtoms * mDomainLength;
    return numer / (numer + mAlpha * mNumBins * (mDomainLength - nAtoms));
}

AtomicProposal createProposal(uint64_t seed, uint64_t id)
{
    GAPS_ASSERT(mMinAtoms == mMaxAtoms);
    AtomicProposal prop(seed, id);

    if (mMinAtoms < 2)
    {
        //#pragma omp atomic
        ++mMaxAtoms;
        prop->type = 'B';
        return prop;
    }

    float u1 = prop.rng.uniform();
    if (u1 < 0.5f)
    {
        if (prop.rng.uniform() < deathProb(static_cast<double>(mMinAtoms)))
        {
            //#pragma omp atomic
            --mMinAtoms;
            prop->type = 'D';
        }
        else
        {
            prop->type = 'B';
        }
        return prop;
    }

    if (u1 < 0.75f)
    {
        prop->type = 'M';
    }
    else
    {
        //#pragma omp atomic
        --mMinAtoms;
        prop->type = 'E'
    }
    return prop;
}

void setNumAtoms(unsigned n)
{
    mMinAtoms = n;
    mMaxAtoms = n;
}

/////////////////////////// ProposalLocationLock ///////////////////////////////

ProposalLocationLock::ProposalLocationLock(unsigned nrow, unsigned npat)
{
    mBinSize = std::numeric_limits<uint64_t>::max() / 
        static_cast<uint64_t>(nrow * npat);
    mNumPatterns = npat;
}

void ProposalLocationLock::fillProposal(AtomicProposal *prop,
AtomicDomain *domain, unsigned id)
{
    switch(prop->type)
    {
        case 'B': // birth

            prop->pos = domain->randomFreePosition(&(prop->rng));
            break;

        case 'D': // death

            prop->atom1 = domain->randomAtom(&(prop->rng));
            break;

        case 'M': // move

            AtomNeighborhood hood = domain->randomAtomWithNeighbors(&(prop->rng));
            uint64_t lbound = hood.hasLeft() ? hood.left->pos : 0;
            uint64_t rbound = hood.hasRight() ? hood.right->pos : mDomainLength;
            prop->pos = prop->rng.uniform64(lbound + 1, rbound - 1);

            if (sameBin(hood.center, prop->pos))
            {
                hood.center->pos = prop->pos;
                prop->type = 'R'; // resolved
            }
            else
            {
                prop->atom1 = hood.center->pos;
            }
            break;

        case 'E': // exchange

            AtomNeighborhood hood = domain->randomAtomWithRightNeighbor(&(prop->rng));
            prop->atom2 = hood.hasRight() ? hood.right : domain.front();

            if (sameBin(hood.center, prop->atom2))
            {
                prop->type = 'R'; // resolved
            }
            else
            {
                prop->atom1 = hood.center->pos;
            }
            break;
    }
}

bool ProposalLocationLock::sameBin(Atom *a1, Atom *a2)
{
    return a1->pos / (mBinSize * mNumPatterns) == a2->pos / (mBinSize * mNumPatterns)
    && (a1->pos / mBinSize) % mNumPatterns == (a2->pos / mBinSize) % mNumPatterns;
}
