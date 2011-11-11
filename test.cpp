/**
 * @file test.cpp
 * SparseMatrix
 * Murtaza Gulamali (11/11/2011)
 *
 * Test SparseMatrix class using googletest, the Google C++ testing framework.
 * http://code.google.com/p/googletest/
 */

#include <iostream>
#include <climits>
#include "gtest/gtest.h"
#include "SparseMatrix.h"

using namespace std;
using namespace myg;

TEST(SparseMatrixTest,DefaultConstructor) {
    unsigned int dim(5);
    SparseMatrix<unsigned int,double> x(dim);
    for (unsigned int i=0; i<dim; i++) {
        for (unsigned int j=0; j<dim; j++) {
            if (i==j)
                EXPECT_DOUBLE_EQ(1.0,x(i,j));
            else
                EXPECT_DOUBLE_EQ(0.0,x(i,j));
        }
    }
}

GTEST_API_ int main(int argc, char* argv[]) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
