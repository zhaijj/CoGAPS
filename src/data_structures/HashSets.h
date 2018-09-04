#ifndef __COGAPS_HASH_SETS_H__
#define __COGAPS_HASH_SETS_H__

#include <stdint.h>
#include <vector>

#include <boost/unordered_set>


class FixedHashSetU32
{
public:

    FixedHashSetU32(unsigned size);

    void insert(unsigned n);
    void erase(unsigned n);
    bool contains(unsigned n);

private:

    std::vector<uint32_t> mSet;
};

class HashSetU64
{
public:

    HashSetU64();

    void insert(uint64_t pos);
    void erase(uint64_t pos);
    bool contains(uint64_t pos);

private:

    boost::unordered_set<uint64_t> mSet;
}

// __COGAPS_HASH_SETS_H__