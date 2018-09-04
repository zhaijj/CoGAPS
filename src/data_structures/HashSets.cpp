#include "HashSets.h"


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

void HashSetU64::insert(uint64_t pos)
{

}

void HashSetU64::erase(uint64_t pos)
{

}

bool HashSetU64::contains(uint64_t pos)
{

}
