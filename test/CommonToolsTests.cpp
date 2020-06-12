//
// Created by artoria on 5/13/20.
//

#include "gtest/gtest.h"


#include "../src/CommonTools.h"

TEST(CommonTools, TransposeRect)
{
    int constexpr len = 8;
    //[1, 2]
    //[3, 4]
    //[5, 6]
    //[7, 8]
    double invec[len] = { 1, 2, 3, 4, 5, 6, 7, 8 };
    //[1, 3, 5, 7]
    //[2, 4, 6, 8]
    double constexpr outvec[len] = { 1, 3, 5, 7, 2, 4, 6, 8 };

    transpose<4, 2>( invec );
    for (int i =0; i < len; ++i)
        ASSERT_TRUE( outvec[i] == invec[i] );
}

TEST(CommonTools, TransposeSquare)
{
    int constexpr len = 9;
    double invec[len] = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    double constexpr outvec[len] = { 1, 4, 7, 2, 5, 8, 3, 6, 9 };

    transpose<3, 3>( invec );
    for (int i =0; i < len; ++i)
        ASSERT_TRUE( outvec[i] == invec[i] );
}