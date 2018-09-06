#include "GibbsSampler.h"
#include "math/Algorithms.h"
#include "math/SIMD.h"

static float getDeltaLL(AlphaParameters alpha, float mass)
{
    return mass * (alpha.s_mu - alpha.s * mass / 2.f);
}

unsigned GibbsSampler::dataRows() const
{
    return mDMatrix.nRow();
}

unsigned GibbsSampler::dataCols() const
{
    return mDMatrix.nCol();
}

void GibbsSampler::setSparsity(float alpha, bool singleCell)
{
    float meanD = singleCell ? gaps::algo::nonZeroMean(mDMatrix) :
        gaps::algo::mean(mDMatrix);

    mPropTypeLock.setAlpha(alpha);
    mLambda = alpha * std::sqrt(mNumPatterns / meanD);
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
    GAPS_ASSERT(mat.nRow() == mMatrix.nRow());
    GAPS_ASSERT(mat.nCol() == mMatrix.nCol());

    mMatrix = mat;
}

void GibbsSampler::setSeed(uint64_t seed)
{
    mSeeder.seed(seed);
}

float GibbsSampler::chi2() const
{   
    return 2.f * gaps::algo::loglikelihood(mDMatrix, mSMatrix, mAPMatrix);
}
  
uint64_t GibbsSampler::nAtoms() const
{   
    return mDomain.size();
}

void GibbsSampler::recalculateAPMatrix()
{
    mAPMatrix = gaps::algo::matrixMultiplication(*mOtherMatrix, mMatrix);
}

void GibbsSampler::sync(const GibbsSampler &sampler)
{
    mOtherMatrix = &(sampler.mMatrix);
    gaps::algo::copyTranspose(&mAPMatrix, sampler.mAPMatrix);
}

void GibbsSampler::update(unsigned nSteps, unsigned nCores)
{
    uint64_t seed = mSeeder.next();

    #pragma omp parallel for num_threads(nCores) schedule(static, 1)
    for (unsigned n = 1; n <= nSteps; ++n)
    {
        AtomicProposal prop(mPropTypeLock.createProposal(seed, n));
        mPropLocationLock.fillProposal(&prop, &mDomain, n);
        processProposal(&prop);
    }
    mPropTypeLock.reset();
    mPropLocationLock.reset();
}

// there should be no data races in anything downstream of this function,
// the locks used to create the proposal ensure that no two conflicting
// proposals are sent to this function at the same time, any aditional
// conflicts or data races must be resolved with critical/atomic pragmas
// or some other method of locking
void GibbsSampler::processProposal(AtomicProposal *prop)
{
    uint64_t pos = 0;
    switch (prop->type)
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

// add an atom at a random position, calculate mass either with an
// exponential distribution or with the gibbs mass distribution
void GibbsSampler::birth(AtomicProposal *prop)
{
    unsigned row = getRow(prop->atom1->pos);
    unsigned col = getCol(prop->atom1->pos);

    // calculate proposed mass
    float mass = canUseGibbs(col)
        ? gibbsMass(alphaParameters(row, col), &(prop->rng)).value()
        : prop->rng.exponential(mLambda);

    // accept mass as long as it's non-zero
    if (mass >= gaps::epsilon)
    {
        // no change to atomic domain
        mPropTypeLock.acceptBirth();
        mPropLocationLock.rejectDeath(prop->atom1->pos);
        prop->atom1->mass = mass;

        // change matrix
        changeMatrix(row, col, mass);
        mPropLocationLock.releaseIndex(row);
    }
    else
    {
        // no change to matrix
        mPropLocationLock.releaseIndex(row);

        // change atomic domain
        mPropTypeLock.rejectBirth();
        uint64_t pos = prop->atom1->pos;
        mDomain.erase(pos);
        mPropLocationLock.acceptDeath(pos);
    }
}

// automatically accept death, attempt a rebirth at the same position, using
// the original mass or the gibbs mass distribution
void GibbsSampler::death(AtomicProposal *prop)
{
    unsigned row = getRow(prop->atom1->pos);
    unsigned col = getCol(prop->atom1->pos);

    // calculate alpha parameters assuming atom dies
    AlphaParameters alpha = alphaParametersWithChange(row, col, -prop->atom1->mass);

    // try to calculate rebirth mass using gibbs distribution, otherwise exponetial
    float rebirthMass = prop->atom1->mass;
    if (canUseGibbs(col))
    {
        OptionalFloat gMass = gibbsMass(alpha, &(prop->rng));
        if (gMass.hasValue())
        {
            rebirthMass = gMass.value();
        }
    }

    // accept/reject rebirth
    float deltaLL = getDeltaLL(alpha, rebirthMass) * mAnnealingTemp;
    if (std::log(prop->rng.uniform()) < deltaLL)
    {
        // no change to atomic domain
        mPropTypeLock.rejectDeath();
        mPropLocationLock.rejectDeath(prop->atom1->pos);
        float origMass = prop->atom1->mass;
        prop->atom1->mass = rebirthMass;

        // change matrix
        if (rebirthMass != origMass)
        {
            safelyChangeMatrix(row, col, rebirthMass - origMass);
        }
        mPropLocationLock.releaseIndex(row);
    }
    else
    {
        // change to atomic domain
        mPropTypeLock.acceptDeath();
        uint64_t pos = prop->atom1->pos;
        mDomain.erase(pos);
        mPropLocationLock.acceptDeath(pos);

        // change to matrix
        safelyChangeMatrix(row, col, -prop->atom1->mass);
        mPropLocationLock.releaseIndex(row);
    }
}

// move mass from one position to another
void GibbsSampler::move(AtomicProposal *prop)
{
    unsigned r1 = getRow(prop->atom1->pos);
    unsigned c1 = getCol(prop->atom1->pos);
    unsigned r2 = getRow(prop->pos);
    unsigned c2 = getCol(prop->pos);

    uint64_t pos = prop->atom1->pos;
    if (r1 == r2 && c1 == c2)
    {
        prop->atom1->pos = prop->pos;
    }
    else
    {
        AlphaParameters alpha = alphaParameters(r1, c1, r2, c2);
        if (std::log(prop->rng.uniform()) < getDeltaLL(alpha, -prop->atom1->mass) * mAnnealingTemp)
        {
            prop->atom1->pos = prop->pos;
            safelyChangeMatrix(r1, c1, -prop->atom1->mass);
            changeMatrix(r2, c2, prop->atom1->mass);
        }
    }
    mPropLocationLock.finishMove(pos);
    mPropLocationLock.releaseIndex(r1);
    mPropLocationLock.releaseIndex(r2);
}

// exchange some amount of mass between two positions, note it is possible
// for one of the atoms to be deleted if it's mass becomes too small
void GibbsSampler::exchange(AtomicProposal *prop)
{
    unsigned r1 = getRow(prop->atom1->pos);
    unsigned c1 = getCol(prop->atom1->pos);
    unsigned r2 = getRow(prop->atom2->pos);
    unsigned c2 = getCol(prop->atom2->pos);

    if (r1 == r2 && c1 == c2)
    {
        mPropLocationLock.releaseIndex(r1);
        mPropLocationLock.releaseIndex(r2);
        return;
    }

    // attempt gibbs distribution exchange
    AlphaParameters alpha = alphaParameters(r1, c1, r2, c2);
    if (canUseGibbs(c1, c2))
    {
        OptionalFloat gMass = gibbsMass(alpha, prop->atom1->mass,
            prop->atom2->mass, &(prop->rng));
        if (gMass.hasValue())
        {
            acceptExchange(prop->atom1, prop->atom2, gMass.value(), r1, c1, r2, c2);
            return;
        }
    }

    // resort to metropolis-hastings if gibbs fails
    exchangeUsingMetropolisHastings(prop, alpha, r1, c1, r2, c2);
}

void GibbsSampler::exchangeUsingMetropolisHastings(AtomicProposal *prop,
AlphaParameters alpha, unsigned r1, unsigned c1, unsigned r2, unsigned c2)
{
    // compute amount of mass to be exchanged
    float totalMass = prop->atom1->mass + prop->atom2->mass;
    float newMass = prop->rng.truncGammaUpper(totalMass, 2.f, 1.f / mLambda);

    // compute amount to change atom1 by - always change larger mass to newMass
    float delta = (prop->atom1->mass > prop->atom2->mass)
        ? newMass - prop->atom1->mass
        : prop->atom2->mass - newMass;

    // choose mass for priorLL calculation
    float oldMass = (2.f * newMass > totalMass)
        ? gaps::max(prop->atom1->mass, prop->atom2->mass)
        : gaps::min(prop->atom1->mass, prop->atom2->mass);

    // calculate priorLL
    float pNew = gaps::d_gamma(newMass, 2.f, 1.f / mLambda);
    float pOld = gaps::d_gamma(oldMass, 2.f, 1.f / mLambda);
    float priorLL = (pNew == 0.f) ? 1.f : pOld / pNew;

    // accept/reject
    float deltaLL = getDeltaLL(alpha, delta) * mAnnealingTemp;
    if (priorLL == 0.f || std::log(prop->rng.uniform() * priorLL) < deltaLL)
    {
        acceptExchange(prop->atom1, prop->atom2, delta, r1, c1, r2, c2);
        return;
    }
    mPropLocationLock.releaseIndex(r1);
    mPropLocationLock.releaseIndex(r2);
}

// helper function for exchange step
void GibbsSampler::acceptExchange(Atom *a1, Atom *a2, float delta,
unsigned r1, unsigned c1, unsigned r2, unsigned c2)
{
    if (a1->mass + delta > gaps::epsilon && a2->mass - delta > gaps::epsilon)
    {
        a1->mass += delta;
        a2->mass -= delta;

        changeMatrix(r1, c1, delta);
        changeMatrix(r2, c2, -delta);
    }
    mPropLocationLock.releaseIndex(r1);
    mPropLocationLock.releaseIndex(r2);
}

// helper function for acceptExchange, used to conditionally update the mass
// at a single position, deleting the atom if the new mass is too small
bool GibbsSampler::updateAtomMass(Atom *atom, float delta)
{
    if (atom->mass + delta < gaps::epsilon)
    {
        mPropTypeLock.acceptDeath();
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
    alpha.s_mu *= mAnnealingTemp;

    if (alpha.s > gaps::epsilon)
    {
        float mean = (alpha.s_mu - mLambda) / alpha.s;
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

OptionalFloat GibbsSampler::gibbsMass(AlphaParameters alpha,
float m1, float m2, GapsRng *rng)
{
    alpha.s *= mAnnealingTemp;
    alpha.s_mu *= mAnnealingTemp;

    if (alpha.s > gaps::epsilon)
    {
        float mean = alpha.s_mu / alpha.s; // lambda cancels out
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

// here mass + delta is guaranteed to be positive
void GibbsSampler::changeMatrix(unsigned row, unsigned col, float delta)
{
    mMatrix(row, col) += delta;
    GAPS_ASSERT(mMatrix(row, col) >= 0.f);
    updateAPMatrix(row, col, delta);
}

// delta could be negative, this is needed to prevent negative values in matrix
void GibbsSampler::safelyChangeMatrix(unsigned row, unsigned col, float delta)
{
    float newVal = gaps::max(mMatrix(row, col) + delta, 0.f);
    updateAPMatrix(row, col, newVal - mMatrix(row, col));
    mMatrix(row, col) = newVal;
}

void GibbsSampler::updateAPMatrix(unsigned row, unsigned col, float delta)
{
    const float *other = mOtherMatrix->colPtr(col);
    float *ap = mAPMatrix.colPtr(row);
    unsigned size = mAPMatrix.nRow();

    gaps::simd::PackedFloat pOther, pAP;
    gaps::simd::Index i(0);
    gaps::simd::PackedFloat pDelta(delta);
    for (; i <= size - i.increment(); ++i)
    {
        pOther.load(other + i);
        pAP.load(ap + i);
        pAP += pDelta * pOther;
        pAP.store(ap + i);
    }

    for (unsigned j = i.value(); j < size; ++j)
    {
        ap[j] += delta * other[j];
    }
}

unsigned GibbsSampler::getRow(uint64_t pos) const
{
    return pos / (mBinSize * mNumPatterns); // nCol == nPatterns
}

unsigned GibbsSampler::getCol(uint64_t pos) const
{
    return (pos / mBinSize) % mNumPatterns; // nCol == nPatterns
}

bool GibbsSampler::canUseGibbs(unsigned col) const
{
    return !gaps::algo::isVectorZero(mOtherMatrix->colPtr(col),
        mOtherMatrix->nRow());
}

bool GibbsSampler::canUseGibbs(unsigned c1, unsigned c2) const
{
    return canUseGibbs(c1) || canUseGibbs(c2);
}

AlphaParameters GibbsSampler::alphaParameters(unsigned row, unsigned col)
{
    return gaps::algo::alphaParameters(mDMatrix.nRow(), mDMatrix.colPtr(row),
        mSMatrix.colPtr(row), mAPMatrix.colPtr(row), mOtherMatrix->colPtr(col));
}

AlphaParameters GibbsSampler::alphaParameters(unsigned r1, unsigned c1,
unsigned r2, unsigned c2)
{
    if (r1 == r2)
    {
        return gaps::algo::alphaParameters(mDMatrix.nRow(), mDMatrix.colPtr(r1),
            mSMatrix.colPtr(r1), mAPMatrix.colPtr(r1), mOtherMatrix->colPtr(c1),
            mOtherMatrix->colPtr(c2));
    }
    return alphaParameters(r1, c1) + alphaParameters(r2, c2);
}

AlphaParameters GibbsSampler::alphaParametersWithChange(unsigned row,
unsigned col, float ch)
{
    return gaps::algo::alphaParametersWithChange(mDMatrix.nRow(),
        mDMatrix.colPtr(row), mSMatrix.colPtr(row), mAPMatrix.colPtr(row),
        mOtherMatrix->colPtr(col), ch);
}

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

#ifdef GAPS_DEBUG
bool GibbsSampler::internallyConsistent()
{
    return true;
}
#endif // GAPS_DEBUG
