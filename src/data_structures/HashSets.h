#ifndef __COGAPS_HASH_SETS_H__
#define __COGAPS_HASH_SETS_H__

#include "../utils/GlobalConfig.h"

#include <stdint.h>
#include <vector>

#include <boost/unordered_set.hpp>

#ifdef __GAPS_OPENMP__
#include <omp.h>
#endif

class FixedHashSetU32
{
public:

    FixedHashSetU32(unsigned size);

    void insert(unsigned n);
    void erase(unsigned n);
    bool contains(unsigned n);
    bool isEmpty();
    unsigned size();

private:

    std::vector<uint32_t> mSet;
};

class SmallHashSetU64
{
public:

    SmallHashSetU64(unsigned capacity);
    ~SmallHashSetU64();

    void insert(uint64_t pos);
    void erase(uint64_t pos);
    bool contains(uint64_t pos);
    bool isEmpty();
    unsigned size();

private:

    boost::unordered_set<uint64_t> mSet;

    #ifdef __GAPS_OPENMP__
    omp_lock_t mLock;
    #endif
};

// internally the same as a SmallHashSetU64 but exposes a slightly different
// interface which allows for signalling by means of creating and destroying
// "canaries" - one thread creates a canary at a given position and another
// thread signals it by destroying the canary, if no canary was destroyed
// then no signal was sent
class CanaryHashSetU64
{
public:
    
    CanaryHashSetU64(unsigned capacity);
    ~CanaryHashSetU64();

    void listen(uint64_t pos); // listen at position
    void signal(uint64_t pos); // signal at position
    bool pop(uint64_t pos); // check if signal was given at position
    bool isEmpty();
    unsigned size();

private:

    boost::unordered_set<uint64_t> mSet;

    #ifdef __GAPS_OPENMP__
    omp_lock_t mLock;
    #endif
};

#endif // __COGAPS_HASH_SETS_H__
