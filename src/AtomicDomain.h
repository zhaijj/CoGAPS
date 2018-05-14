#ifndef __COGAPS_ATOMIC_DOMAIN_H__
#define __COGAPS_ATOMIC_DOMAIN_H__

#include "Archive.h"

#include <boost/unordered_set.hpp>
#include <stdint.h>
#include <cstddef>
#include <vector>
#include <map>

class AtomicDomain;

// Want the neighbors to be pointers, but can't use pointers since vector
// is resized, shifted integers allow for 0 to be "null" and 1 to be the
// first index. AtomicDomain is responsible for managing this correctly.
struct Atom
{
//private:    
public:

    friend AtomicDomain;
    uint64_t leftNdx; // shift these up 1, 0 == no neighbor
    uint64_t rightNdx;

public:    

    uint64_t pos;
    float mass;

    Atom(uint64_t p, float m)
        : pos(p), mass(m), leftNdx(0), rightNdx(0)
    {}

    bool operator==(const Atom &other) const
    {
        return pos == other.pos;
    }

    bool operator!=(const Atom &other) const
    {
        return !(pos == other.pos);
    }
};

struct RawAtom
{
    uint64_t pos;
    float mass;

    RawAtom() : pos(0), mass(0.f) {}
    RawAtom(uint64_t p, float m) : pos(p), mass(m) {}
};

// data structure that holds atoms
class AtomicDomain
{
private:

    // size of atomic domain
    uint64_t mDomainSize;

    // domain storage
    std::vector<Atom> mAtoms;
    std::map<uint64_t, uint64_t> mAtomPositions;

    // TODO google_dense_set - first profile and benchmark
    boost::unordered_set<uint64_t> mUsedPositions;

    mutable std::vector<RawAtom> mInsertCache;
    mutable std::vector<uint64_t> mEraseCache;

    mutable unsigned mInsertCacheIndex;
    mutable unsigned mEraseCacheIndex;

    Atom& _left(const Atom &atom);
    Atom& _right(const Atom &atom);

public:

    void setDomainSize(uint64_t size) { mDomainSize = size; }

    // access atoms
    Atom front() const;
    Atom randomAtom() const;
    uint64_t randomFreePosition() const;
    uint64_t size() const;

    const Atom& left(const Atom &atom) const;
    const Atom& right(const Atom &atom) const;
    bool hasLeft(const Atom &atom) const;
    bool hasRight(const Atom &atom) const;

    // modify domain
    Atom insert(uint64_t pos, float mass);
    void erase(uint64_t pos);
    void cacheInsert(uint64_t pos, float mass) const;
    void cacheErase(uint64_t pos) const;
    void updateMass(uint64_t pos, float newMass);
    void flushCache();
    void resetCache(unsigned n);

    // serialization
    friend Archive& operator<<(Archive &ar, AtomicDomain &domain);
    friend Archive& operator>>(Archive &ar, AtomicDomain &domain);
};

#endif