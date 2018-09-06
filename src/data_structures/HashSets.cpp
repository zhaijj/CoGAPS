#include "HashSets.h"

#ifdef __GAPS_OPENMP__

    #include <omp.h>
    #define gaps_omp_init_lock(x)       omp_init_lock(x)
    #define gaps_omp_destroy_lock(x)    omp_destroy_lock(x)
    #define gaps_omp_set_lock(x)        omp_set_lock(x)
    #define gaps_omp_unset_lock(x)      omp_unset_lock(x)

#else

    #define gaps_omp_init_lock(x)       do {} while(0)
    #define gaps_omp_destroy_lock(x)    do {} while(0)
    #define gaps_omp_set_lock(x)        do {} while(0)
    #define gaps_omp_unset_lock(x)      do {} while(0)

#endif

///////////////////////////// FixedHashSetU32 //////////////////////////////////

FixedHashSetU32::FixedHashSetU32(unsigned size)
    : mSet(std::vector<uint32_t>(size, 0))
{}

void FixedHashSetU32::insert(unsigned n)
{
    mSet[n] = 1;
}

void FixedHashSetU32::erase(unsigned n)
{
    mSet[n] = 0;
}

bool FixedHashSetU32::contains(unsigned n)
{
    return mSet[n] == 1;
}

bool FixedHashSetU32::isEmpty()
{
    unsigned sz = mSet.size();
    for (unsigned i = 0; i < sz; ++i)
    {
        if (mSet[i] == 1)
        {
            return false;
        }
    }
    return true;
}

unsigned FixedHashSetU32::size()
{
    unsigned count = 0;
    unsigned sz = mSet.size();
    for (unsigned i = 0; i < sz; ++i)
    {
        count += mSet[i];
    }
    return count;
}

///////////////////////////// SmallHashSetU64 //////////////////////////////////

SmallHashSetU64::SmallHashSetU64(unsigned capacity)
{
    gaps_omp_init_lock(&mLock);
}

SmallHashSetU64::~SmallHashSetU64()
{
    gaps_omp_destroy_lock(&mLock);
}

bool SmallHashSetU64::isEmpty()
{
    gaps_omp_set_lock(&mLock);
    bool emp = mSet.empty();
    gaps_omp_unset_lock(&mLock);
    return emp;
}

unsigned SmallHashSetU64::size()
{
    gaps_omp_set_lock(&mLock);
    unsigned sz = mSet.size();
    gaps_omp_unset_lock(&mLock);
    return sz;
}

void SmallHashSetU64::insert(uint64_t pos)
{
    gaps_omp_set_lock(&mLock);
    mSet.insert(pos);
    gaps_omp_unset_lock(&mLock);
}

void SmallHashSetU64::erase(uint64_t pos)
{
    gaps_omp_set_lock(&mLock);
    mSet.erase(pos);
    gaps_omp_unset_lock(&mLock);
}

bool SmallHashSetU64::contains(uint64_t pos)
{
    gaps_omp_set_lock(&mLock);
    unsigned ct = mSet.count(pos);
    gaps_omp_unset_lock(&mLock);
    return ct > 0;
}

//////////////////////////// CanaryHashSetU64 //////////////////////////////////

CanaryHashSetU64::CanaryHashSetU64(unsigned capacity)
{
    gaps_omp_init_lock(&mLock);
}

CanaryHashSetU64::~CanaryHashSetU64()
{
    gaps_omp_destroy_lock(&mLock);
}

bool CanaryHashSetU64::isEmpty()
{
    gaps_omp_set_lock(&mLock);
    bool emp = mSet.empty();
    gaps_omp_unset_lock(&mLock);
    return emp;
}

unsigned CanaryHashSetU64::size()
{
    gaps_omp_set_lock(&mLock);
    unsigned sz = mSet.size();
    gaps_omp_unset_lock(&mLock);
    return sz;
}

void CanaryHashSetU64::listen(uint64_t pos)
{
    gaps_omp_set_lock(&mLock);
    mSet.insert(pos);
    gaps_omp_unset_lock(&mLock);
}

void CanaryHashSetU64::signal(uint64_t pos)
{
    gaps_omp_set_lock(&mLock);
    mSet.erase(pos);
    gaps_omp_unset_lock(&mLock);
}

bool CanaryHashSetU64::pop(uint64_t pos)
{
    gaps_omp_set_lock(&mLock);
    unsigned ct = mSet.count(pos);
    if (ct > 0)
    {
        mSet.erase(pos);
    }
    gaps_omp_unset_lock(&mLock);
    return ct == 0;
}
