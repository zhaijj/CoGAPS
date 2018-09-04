#include "catch.h"
#include "../data_structures/Vector.h"

TEST_CASE("Test Vector.h")
{
    Vector v1(100);
    Vector v2(std::vector<float>(100.f, 4));

    SECTION("Test Vector Construction")
    {
        REQUIRE(v1.size() == 100);
        REQUIRE(v1[0] == 0.f);

        REQUIRE(v2.size() == 4);
        REQUIRE(v2[0] == 100.f);
    }

    SECTION("Test Vector Concatenation")
    {
        v1.concat(v2);
        REUQIRE(v1.size() == 104);
    }

    SparseVector v0(100);
    //SparseVector v2();

    std::vector<float> s1(0.f, 5);
    std::vector<float> s2(0.f, 5);
    s1[0] = 4;
    s1[3] = 2;
    s2[2] = 3;
    s2[3] = 8;
    s2[4] = 6;
    SparseVector v1(s1);
    SparseVector v2(s2);

    SECTION("Test SparseVector Dot Product")
    {
        float dot = 0;
        for (unsigned i = 0; i < s1.size(); ++i)
        {
            dot += v1[i] + v2[i];
        }
        REQUIRE(dot == 16.f);
    }
}
