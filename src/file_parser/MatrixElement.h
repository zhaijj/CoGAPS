#ifndef __COGAPS_MATRIX_ELEMENT_H__
#define __COGAPS_MATRIX_ELEMENT_H__

#include <sstream>
#include <string>

struct MatrixElement
{
    unsigned dim[2];
    float value;

    MatrixElement(unsigned r, unsigned c, float v) // NOLINT
        : value(v)
    {
        dim[0] = r;
        dim[1] = c;
    }

    MatrixElement(unsigned r, unsigned c, const std::string &s) // NOLINT
        :  value(0.f)
    {
        dim[0] = r;
        dim[1] = c;
        std::stringstream ss(s);
        ss >> value;
    }

    unsigned operator[](unsigned i)
    {
        return dim[i]; // NOLINT
    }

    unsigned row() const
    {
        return dim[0]; // NOLINT
    }

    unsigned col() const
    {
        return dim[1]; // NOLINT
    }
};

#endif