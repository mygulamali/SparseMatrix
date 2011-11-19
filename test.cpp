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
#include <vector>
#include "gtest/gtest.h"
#include "SparseMatrix.h"

using namespace std;
using namespace myg;

//! class to store fixture data for tests
class SparseMatrixTest : public ::testing::Test {
protected:
    SparseMatrixTest() :
    _a(SparseMatrix<unsigned int,double>(4)),
    _b(SparseMatrix<unsigned int,double>(4))
    {}
    
    ~SparseMatrixTest()
    {}
    
    virtual void SetUp() {
        /*
         * [[0.5 0.0 0.0 0.0]
         *  [1.0 2.0 0.0 0.0]
         *  [0.0 0.0 0.5 1.5]
         *  [0.0 0.0 0.0 0.5]]
         */
        _a(1,0,2.0);
        _a(1,1,4.0);
        _a(2,3,3.0);
        _a /= 2.0;
        
        /*
         * [[2.0 4.0 0.0 0.0]
         *  [0.0 2.0 0.0 6.0]
         *  [0.0 0.0 2.0 0.0]
         *  [0.0 0.0 0.0 2.0]]
         */
        _b(0,1,2.0);
        _b(1,3,3.0);
        _b *= 2.0;
    }
    
    //virtual void TearDown() {}
    
    SparseMatrix<unsigned int,double> _a;
    SparseMatrix<unsigned int,double> _b;
};

// TESTS!

TEST_F(SparseMatrixTest,DefaultConstructor) {
    unsigned int dim(3);
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

/*
TEST_F(SparseMatrixTest,SetZeroElement) {
}

TEST_F(SparseMatrixTest,SetNonZeroElement) {
}

TEST_F(SparseMatrixTest,SetZeroDiagonalElement) {
}
*/

TEST_F(SparseMatrixTest,MatrixAddition) {
    SparseMatrix<unsigned int, double> e(4), o(4);
	e(0,0,2.5);
	e(0,1,4.0);
	e(1,0,1.0);
	e(1,1,4.0);
	e(1,3,6.0);
	e(2,2,2.5);
	e(2,3,1.5);
	e(3,3,2.5);
	
    // perform matrix addition
    o = _a+_b;
    
    // perform test
    for (unsigned int i=0; i<4; i++)
        for (unsigned int j=0; j<4; j++)
            EXPECT_DOUBLE_EQ(e(i,j),o(i,j));
}

TEST_F(SparseMatrixTest,MatrixSubtraction) {
    SparseMatrix<unsigned int, double> e(4), o(4);
	e(0,0,-1.5);
	e(0,1,-4.0);
	e(1,0, 1.0);
	e(1,1, 0.0);
	e(1,3,-6.0);
	e(2,2,-1.5);
	e(2,3, 1.5);
	e(3,3,-1.5);
	
    // perform matrix subtraction
    o = _a-_b;
    
    // perform test
    for (unsigned int i=0; i<4; i++)
        for (unsigned int j=0; j<4; j++)
            EXPECT_DOUBLE_EQ(e(i,j),o(i,j));
}

TEST_F(SparseMatrixTest,MatrixMultiplication) {
    SparseMatrix<unsigned int, double> e(4), o(4);
	e(0,0, 1.0);
	e(0,1, 2.0);
	e(1,0, 2.0);
	e(1,1, 8.0);
	e(1,3,12.0);
	e(2,2, 1.0);
	e(2,3, 3.0);
	e(3,3, 1.0);
	
    // perform matrix multiplication
    o = _a*_b;
    
    // perform test
    for (unsigned int i=0; i<4; i++)
        for (unsigned int j=0; j<4; j++)
            EXPECT_DOUBLE_EQ(e(i,j),o(i,j));
}

TEST_F(SparseMatrixTest,VectorSolve) {
    /*
     *  [[1  2  0  0  0  0];
     *   [7  4  3  0  0  0];
     *   [0  0  6  0  0  0];
     *   [9  0  0  1  0  2];
     *   [0  0  8  7  1  0];
     *   [0  0  0  1  0  1]]
     */
    SparseMatrix<unsigned int, double> A(6);
    A(0,1,2.0);
    A(1,0,7.0);
    A(1,1,4.0);
    A(1,2,3.0);
    A(2,2,6.0);
    A(3,0,9.0);
    A(3,5,2.0);
    A(4,2,8.0);
    A(4,3,7.0);
    A(5,3,1.0);        

	// [1 2 3 4 5 6].'
	vector<double> b;
	for (unsigned int i=1; i<=6; i++)
        b.push_back((double) i);
    
	// solve equation Ax = b for x
	vector<double> x = A.solve(b);
	
	// compute residual vector
	vector<double> r = A*x;
	for (unsigned int i=0; i<r.size(); i++)
        r[i] -= b[i];

    // perform test
    for (unsigned int i=0; i<r.size(); i++)
        EXPECT_NEAR(0.0,r[i],1.0e-12);
}

GTEST_API_ int main(int argc, char* argv[]) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
