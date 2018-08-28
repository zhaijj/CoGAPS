#include "Vector.h"

Vector::Vector(unsigned size) : mValues(aligned_vector(size, 0.f)) {}

Vector::Vector(const std::vector<float> &v) : mValues(v.size())
{
    for (unsigned i = 0; i < v.size(); ++i)
    {
        mValues[i] = v[i];
    }
}

unsigned Vector::size() const
{
    return mValues.size();
}

float* Vector::ptr()
{
    GAPS_ASSERT(size() > 0);
    return &mValues[0];
}

const float* Vector::ptr() const
{
    GAPS_ASSERT(size() > 0);
    return &mValues[0];
}

float& Vector::operator[](unsigned i)
{
    GAPS_ASSERT(i < size());
    return mValues[i];
}

float Vector::operator[](unsigned i) const
{
    GAPS_ASSERT(i < size());
    return mValues[i];
}

void Vector::operator+=(const Vector &vec)
{
    GAPS_ASSERT(vec.size() == size());
    for (unsigned i = 0; i < size(); ++i)
    {
        mValues[i] += vec[i];
    }
}

void Vector::operator*=(float val)
{
    for (unsigned i = 0; i < mValues.size(); ++i)
    {
        mValues[i] *= val;
    }
}

void Vector::operator/=(float val)
{
    for (unsigned i = 0; i < size(); ++i)
    {
        mValues[i] /= val;
    }
}

Vector Vector::operator-(Vector vec) const
{
    GAPS_ASSERT(vec.size() == size());
    for (unsigned i = 0; i < size(); ++i)
    {
        vec[i] = mValues[i] - vec[i];
    }
    return vec;
}

Vector Vector::operator*(float val) const
{
    Vector vec(*this);
    vec *= val;
    return vec;
}

Vector Vector::operator/(float val) const
{
    Vector vec(*this);
    vec /= val;
    return vec;
}

Archive& operator<<(Archive &ar, Vector &vec)
{
    ar << vec.size();
    for (unsigned i = 0; i < vec.size(); ++i)
    {
        ar << vec[i];
    }
    return ar;
}

Archive& operator>>(Archive &ar, Vector &vec)
{
    // should already be allocated
    unsigned sz = 0;
    ar >> sz;
    GAPS_ASSERT(sz == vec.size());

    // read in data
    for (unsigned i = 0; i < vec.size(); ++i)
    {
        ar >> vec.mValues[i];
    }
    return ar;
}
