/**
 * @file SparseMatrixDemo.cpp
 * SparseMatrix
 * Murtaza Gulamali (02/12/2006)
 *
 * Simple demo to set up a sparse matrix object A and vector b, and use them
 * to solve the equation Ax = b for vector x.
 */

#include <iostream>
#include <stdio.h>
#include <vector>
#include "SparseMatrix.h"

using namespace std;
using namespace myg;

int main(int argc, char* argv[]) {
	/* setup matrix A
	 * [[1  2  0  0  0  0];
	 *  [7  4  3  0  0  0];
	 *  [0  0  6  0  0  0];
	 *  [9  0  0  1  0  2];
	 *  [0  0  8  7  1  0];
	 *  [0  0  0  1  0  1]]
	 */
	SparseMatrix<unsigned int,double> A(6); // 6x6 identity matrix
	A(0,1,2.);
	A(1,0,7.);
	A(1,1,4.);
	A(1,2,3.);
	A(2,2,6.);
	A(3,0,9.);
	A(3,5,2.);
	A(4,2,8.);
	A(4,3,7.);
	A(5,3,1.);
    
	// setup vector b
	// [1 2 3 4 5 6]^T
	vector<double> b;
	for (int i=1; i<=6; i++)
		b.push_back((double) i);
    
	// solve equation Ax = b for x
	vector<double> x = A.solve(b);
    
	// output results to stdout
	cout << "Matrix A:" << "\n" << A << endl;
	cout << "Vector b:" << "\n";
	for (unsigned int i=0; i<b.size(); i++)
		cout << b.at(i) << endl;
	cout << endl << "Vector x:" << "\n";
	for (unsigned int i=0; i<x.size(); i++)
		cout << x.at(i) << endl;
    
    // return
    return 0;
}
