#include "GibbsSampler.h"

// TODO don't use this
unsigned GibbsSampler::dataRows() const
{
    return mDMatrix.nRow();
}

// TODO don't use this
unsigned GibbsSampler::dataCols() const
{
    return mDMatrix.nCol();
}

void GibbsSampler::setSparsity(float alpha, bool singleCell)
{
    float meanD = singleCell ? gaps::algo::nonZeroMean(mDMatrix) :
        gaps::algo::mean(mDMatrix);

    mAlpha = alpha;
    mLambda = alpha * std::sqrt(mMatrix.nCol() / meanD); // nPatterns / meanD
}

void GibbsSampler::setMaxGibbsMass(float max)
{
    mMaxGibbsMass = max;
}

void GibbsSampler::setAnnealingTemp(float temp)
{
    mAnnealingTemp = temp;
}

void GibbsSampler::setMatrix(const Matrix &mat)
{   
    mMatrix = mat;
}

float GibbsSampler::chi2() const
{   
    return 2.f * gaps::algo::loglikelihood(mDMatrix, mSMatrix, mAPMatrix);
}
  
uint64_t GibbsSampler::nAtoms() const
{   
    return mDomain.size();
}

void GibbsSampler::update(unsigned nSteps, unsigned nCores)
{
    for (unsigned n = 0; n < nSteps; ++n)
    {
        makeAndProcessProposal(); // use n for rng stream
    }
}

void GibbsSampler::makeAndProcessProposal()
{
    GapsRng rng; // need to make this independent of nThreads (use streams?)
    if (mDomain.size() < 2) // always birth with 0 or 1 atoms
    {
        return birth();
    }

    float u1 = rng.uniform();
    if (u1 <= 0.5f)
    {
        return rng.uniform() < deathProb(mDomain.size()) ? death(&rng) : birth(&rng);
    }
    return (u1 < 0.75f) ? move(&rng) : exchange(&rng);
}

float GibbsSampler::deathProb(uint64_t nAtoms) const
{
    double size = static_cast<double>(mDomainLength);
    double term1 = (size - static_cast<double>(nAtoms)) / size;
    double term2 = mAlpha * static_cast<double>(mNumBins) * term1;
    return static_cast<double>(nAtoms) / (static_cast<double>(nAtoms) + term2);
}

// add an atom at a random position, calculate mass either with an
// exponential distribution or with the gibbs mass calculation
void GibbsSampler::birth(GapsRng *rng)
{
    // get random position and matrix indices 
    uint64_t pos = mDomain.randomFreePosition(rng);
    unsigned mainNdx = mainIndex(pos);
    unsigned patNdx = patternIndex(pos);

    // calculate proposed mass
    float mass = canUseGibbs(patNdx)
        ? gibbsMass(alphaParameters(mainNdx, patNdx), rng).value()
        : rng->exponential(mLambda);
    
    // accept mass as long as it's non-zero
    if (mass >= gaps::epsilon)
    {
        mDomain.insert(pos, mass);
        mMatrix(mainNdx, patNdx) += mass;
        updateAPMatrix(mainNdx, patNdx, mass);
    }
}

// automatically accept death, attempt a rebirth at the same position, using
// the original mass or the gibbs mass calculation
void GibbsSampler::death(GapsRng *rng)
{
    // calculate matrix indices for this atom
    Atom* atom = mDomain.randomAtom(rng);
    unsigned mainNdx = mainIndex(pos);
    unsigned pattNdx = patternIndex(pos);

    // kill off atom
    safelyChangeMatrixValue(mainNdx, pattNdx, -1.f * atom->mass);

    // calculate rebirth mass
    AlphaParameters alpha = alphaParameters(mainNdx, pattNdx);
    float rebirthMass = atom->mass;
    if (canUseGibbs(pattNdx))
    {
        OptionalFloat mass = gibbsMass(alpha, rng);
        if (mass.hasValue())
        {
            rebirthMass = mass.value();
        }
    }

    // accept/reject rebirth
    float deltaLL = rebirthMass * (alpha.su - alpha.s * rebirthMass / 2.f);
    if (std::log(rng->uniform()) < deltaLL * mAnnealingTemp)
    {
        atom->mass = rebirthMass;
        mMatrix(mainNdx, pattNdx) += rebirthMass;
        updateAPMatrix(mainNdx, pattNdx, rebirthMass);
    }
    else
    {
        mDomain.erase(atom->pos);
    }
}

// move mass from src to dest in the atomic domain
void GibbsSampler::move(GapsRng *rng)
{
    // get atom and find the distance between it's neighbors
    AtomNeighborhood hood = mDomain.randomAtomWithNeighbors(rng);
    uint64_t lbound = hood.hasLeft() ? hood.left->pos : 0;
    uint64_t rbound = hood.hasRight() ? hood.right->pos : mDomainLength;

    // choose a new location to move the atom to, don't move past neighbors
    uint64_t newLocation = mRng.uniform64(lbound + 1, rbound - 1);

    // get matrix indices for both positions
    unsigned mainNdx1 = mainIndex(hood.center->pos);
    unsigned mainNdx2 = mainIndex(newLocation);
    unsigned pattNdx1 = patternIndex(hood.center->pos);
    unsigned pattNdx2 = patternIndex(newLocation);

    // automatically accept if in same bin
    if (mainNdx1 == mainNdx2 && pattNdx1 == pattNdx2)
    {
        hood.center->pos = newLocation;
        return;
    }

    // evaluate proposed change
    float deltaLL = computeDeltaLL(mainNdx1, pattNdx1, -1.f * prop.atom1->mass,
        mainNdx2, pattNdx2, prop.atom1->mass);
    if (std::log(rng->uniform()) < deltaLL * mAnnealingTemp)
    {
        hood.center->pos = newLocation;
        safelyChangeMatrixValue(mainNdx1, pattNdx1, -1.f * hood.center->mass);
        mMatrix(mainNdx2, pattNdx2) += hood.center->mass;
        updateAPMatrix(mainNdx2, pattNdx2, hood.center->mass);
    }
}

// exchange some amount of mass between two positions, note it is possible
// for one of the atoms to be deleted if it's mass becomes too small
void GibbsSampler::exchange(GapsRng *rng)
{
    AtomNeighborhood hood = mDomain.randomAtomWithRightNeighbor(rng);
    Atom* exchangeAtom = hood.hasRight() ? hood.right : mDomain.front();

    // get matrix indices for both atom positions
    unsigned mainNdx1 = mainIndex(hood.center->pos);
    unsigned mainNdx2 = mainIndex(exchangeAtom->pos);
    unsigned pattNdx1 = patternIndex(hood.center->pos);
    unsigned pattNdx2 = patternIndex(exchangeAtom->pos);

    // TODO automatically accept if in the same bin
    if (mainNdx1 == mainNdx2 && pattNdx1 == pattNdx2)
    {
        return;
    }

    float m1 = hood.center->mass;
    float m2 = exchangeAtom->mass;

    if (canUseGibbs(pattNdx1, pattNdx2))
    {
        AlphaParameters alpha = alphaParameters(mainNdx1, pattNdx1,
            mainNdx2, pattNdx2);
        OptionalFloat delta = gibbsMass(alpha, m1, m2, rng);
        if (delta.hasValue())
        {
            acceptExchange(hood.center, exchangeAtom, delta, mainNdx1,
                pattNdx1, mainNdx2, pattNdx2);
            return;
        }
    }

    float newMass = rng->truncGammaUpper(m1 + m2, 2.f, 1.f / mLambda);

    float delta = m1 > m2 ? newMass - m1 : m2 - newMass; // change larger mass
    float pOldMass = 2.f * newMass > m1 + m2 ? gaps::max(m1, m2) : gaps::min(m1, m2);

    float pNew = gaps::d_gamma(newMass, 2.f, 1.f / mLambda);
    float pOld = gaps::d_gamma(pOldMass, 2.f, 1.f / mLambda);

    if (pOld == 0.f && pNew != 0.f) // special case
    {
        acceptExchange(hood.center, exchangeAtom, delta, mainNdx1, pattNdx1,
            mainNdx2, pattNdx2);
        return;
    }

    float deltaLL = computeDeltaLL(mainNdx1, pattNdx1, delta, mainNdx2,
        pattNdx2, -1.f * delta);
    float priorLL = (pOld == 0.f) ? 1.f : pOld / pNew;
    if (std::log(prop.rng.uniform() * priorLL) < deltaLL * mAnnealingTemp)
    {
        acceptExchange(hood.center, exchangeAtom, delta, mainNdx1, pattNdx1,
            mainNdx2, pattNdx2);
        return;
    }
}

// helper function for exchange step
void GibbsSampler::acceptExchange(Atom *atom1, Atom *atom2, float delta,
unsigned r1, unsigned c1, unsigned r2, unsigned c2)
{
    // if updating the atom mass makes it too small, erase the entire atom
    float d1 = conditionallyUpdateAtom(atom1, delta) ? delta : -atom1->mass;
    float d2 = conditionallyUpdateAtom(atom2, -delta) ? -delta : -atom2->mass;

    // ensure matrix values don't go negative
    safelyChangeMatrixValue(r1, c1, d1);
    safelyChangeMatrixValue(r2, c2, d2);
}

// helper function for acceptExchange, used to conditionally update the mass
// at a single position, deleting the atom if the new mass is too small
bool GibbsSampler::conditionallyUpdateAtom(Atom *atom, float delta)
{
    if (atom->mass + delta < gaps::epsilon)
    {
        mDomain.erase(atom->pos);
        return false;
    }
    atom->mass += delta;
    return true;
}

OptionalFloat GibbsSampler::gibbsMass(AlphaParameters alpha,
GapsRng *rng)
{        
    alpha.s *= mAnnealingTemp;
    alpha.su *= mAnnealingTemp;

    if (alpha.s > gaps::epsilon)
    {
        float mean = (alpha.su - mLambda) / alpha.s;
        float sd = 1.f / std::sqrt(alpha.s);
        float pLower = gaps::p_norm(0.f, mean, sd);

        if (pLower < 1.f)
        {
            float m = rng->inverseNormSample(pLower, 1.f, mean, sd);
            float gMass = gaps::min(m, mMaxGibbsMass / mLambda);
            if (gMass >= gaps::epsilon)
            {
                return OptionalFloat(gMass);
            }
        }
    }
    return OptionalFloat();
}

std::pair<float, bool> GibbsSampler::gibbsMass(AlphaParameters alpha,
float m1, float m2, GapsRng *rng)
{
    alpha.s *= mAnnealingTemp;
    alpha.su *= mAnnealingTemp;

    if (alpha.s > gaps::epsilon)
    {
        float mean = alpha.su / alpha.s; // lambda cancels out
        float sd = 1.f / std::sqrt(alpha.s);
        float pLower = gaps::p_norm(-m1, mean, sd);
        float pUpper = gaps::p_norm(m2, mean, sd);

        if (!(pLower >  0.95f || pUpper < 0.05f))
        {
            float delta = rng->inverseNormSample(pLower, pUpper, mean, sd);
            float gMass = gaps::min(gaps::max(-m1, delta), m2); // conserve mass
            return OptionalFloat(gMass);
        }
    }
    return OptionalFloat();
}

void GibbsSampler::sync(const GibbsSampler &src)
{
    mOtherMatrix = &(src.mMatrix);
    mAPMatrix = src.mAPMatrix;
}

unsigned GibbsSampler::mainIndex(uint64_t pos) const
{
    return pos / (mBinSize * mNumPatterns);
}

unsigned GibbsSampler::patternIndex(uint64_t pos) const
{
    return (pos / mBinSize) % mNumPatterns;
}

bool GibbsSampler::canUseGibbs(unsigned pattNdx) const
{
    return !gaps::algo::isVectorZero(mOtherMatrix->colPtr(pattNdx),
        mOtherMatrix->nRow());
}

bool GibbsSampler::canUseGibbs(unsigned pattNdx1, unsigned pattNdx2) const
{
    return canUseGibbs(pattNdx1) || canUseGibbs(pattNdx2);
}

AlphaParameters GibbsSampler::alphaParameters(unsigned mainNdx, unsigned patNdx)
{
    return gaps::algo::alphaParameters(mDMatrix.nRow(), mDMatrix.colPtr(main),
        mSMatrix.colPtr(main), mAPMatrix.colPtr(main), mOtherMatrix->colPtr(pat));
}

float AmplitudeGibbsSampler::computeDeltaLL(unsigned mainNdx, unsigned patNdx,
float mass)
{
    return gaps::algo::deltaLL(mDMatrix.nRow(), mDMatrix.colPtr(main),
        mSMatrix.colPtr(main), mAPMatrix.colPtr(main), mOtherMatrix->colPtr(pat),
        mass);
}

// needed to prevent negative values in matrix
void GibbsSampler::safelyChangeMatrixValue(unsigned row, unsigned col, float d)
{
    float newVal = gaps::max(mMatrix(row, col) - d, 0.f);
    updateAPMatrix(row, col, newVal - mMatrix(row, col));
    mMatrix(row, col) = newVal;
}

#ifdef GAPS_DEBUG
bool GibbsSampler::internallyConsistent()
{
    return true;
}
#endif // GAPS_DEBUG

Archive& operator<<(Archive &ar, GibbsSampler &s)
{
    // TODO
    return ar;
}

Archive& operator>>(Archive &ar, GibbsSampler &s)
{
    // TODO
    return ar;
}