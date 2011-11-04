/**
 * @file SparseMatrix.cpp
 * SparseMatrix
 * Murtaza Gulamali (30/11/2006)
 *
 * This software is released under the terms and conditions of The MIT License:
 * http://www.opensource.org/licenses/mit-license.php
 */

#include <cmath>
#include <iterator>
#include <algorithm>

#include "SparseMatrix.h"

using namespace std;
using namespace myg;

// static const values

//! zero
template <typename iT, typename fT>
const fT SparseMatrix<iT,fT>::zero_(0.0);

//! epsilon: a very very small number number greater than zero
template <typename iT, typename fT>
const fT SparseMatrix<iT,fT>::eps_(1.0e-16);

// constructors and destructors

/**
 * default constructor creates identity matrix of specified dimension
 * @param dim dimension of matrix (default is 1)
 */
template <typename iT, typename fT>
SparseMatrix<iT,fT>::SparseMatrix(iT dim) :
nnu_(dim)
{	
	// reserve necessary memory
	ia_.reserve(dim+1);
	ja_.reserve(dim);
	a_.reserve(dim);
	
	// set to identity matrix
	for (iT i=0; i<dim; i++) {
		ia_.push_back(i);
		ja_.push_back(i);
		a_.push_back(static_cast<fT>(1.));
	}
	ia_.push_back(dim);    
}

/**
 * copy constructor
 * @param m SparseMatrix object to copy
 */
template <typename iT, typename fT>
SparseMatrix<iT,fT>::SparseMatrix(const SparseMatrix<iT,fT>& m) :
nnu_(m.Dim()),
ia_(m.IA()),
ja_(m.JA()),
a_(m())
{ }

//! default destructor
template <typename iT, typename fT>
SparseMatrix<iT,fT>::~SparseMatrix() {}

// accessors and mutators for matrix elements

/**
 * return value at specified indices of matrix
 * @param row matrix row
 * @param col matrix column
 * @return value at specified matrix row and column
 */
template <typename iT, typename fT>
const fT& SparseMatrix<iT,fT>::operator()(iT row, iT col) const {
	if (row==col)
        return a_.at(ia_.at(row));
	
	const iT j = findIndex(row,col);
	if (j<ia_.at(row+1))
		return a_.at(j);
	else
		return zero_;		
}

/**
 * set value at specified indices of matrix
 * @param row matrix row
 * @param col matrix column
 * @param val value at specified matrix row and column
 */
template <typename iT, typename fT>
void SparseMatrix<iT,fT>::operator()(iT row, iT col, fT val) {
	if (row==col)
		// set diagonal element
		a_.at(ia_.at(row)) = val;		
	else {
		const iT j = findIndex(row,col);
		if (j<ia_.at(row+1)) {
			// non-zero at (row,col)
			if (fabs(val)>eps_) {
				// change non-zero
				a_.at(j) = val;
			} else {
				// remove non-zero
				for (iT i=(row+1); i<=Dim(); i++) ia_.at(i)--;
				ja_.erase(ja_.begin()+j);
				a_.erase(a_.begin()+j);
			}
		} else {
			// zero at (row,col)
			if (fabs(val)>eps_) {
				// insert non-zero
				for (iT i=(row+1); i<=Dim(); i++) ia_.at(i)++;
				ja_.insert(ja_.begin()+(ia_.at(row+1)-1),col);
				a_.insert(a_.begin()+(ia_.at(row+1)-1),val);
			}
		}
	}
}

// matrix operator methods
/**
 * matrix assignment
 * @param m assign this to m (ie. produce a "deep" copy)
 * @return copy
 */
template <typename iT, typename fT>
SparseMatrix<iT,fT>& SparseMatrix<iT,fT>::operator=(const SparseMatrix<iT,fT>& m) {
	// don't assign to itself
	if (this != &m) {
		// perform copy
		Dim(m.Dim());
		IA(m.IA());
		JA(m.JA());
		(*this)(m());
	}
	
	// return reference
	return (*this);	
}

/**
 * matrix addition
 * @param m matrix to add to this
 * @return result of matrix addition
 */
template <typename iT, typename fT>
SparseMatrix<iT,fT> SparseMatrix<iT,fT>::operator+(const SparseMatrix<iT,fT>& m) const {
	// initialise local variables
	typename std::set<std::pair<iT,iT> > sp = m.sparsityPattern();
	SparseMatrix<iT,fT> ret(*this);
    
	// perform addition
	for (typename std::set<std::pair<iT,iT> >::iterator i=sp.begin(); i!=sp.end(); i++)
		ret((*i).first,(*i).second,
			ret((*i).first,(*i).second) + m((*i).first,(*i).second));
	
	// return result
	return ret;
}

/**
 * matrix subtraction
 * @param m matrix to subtract from this
 * @return result of matrix subtraction
 */
template <typename iT, typename fT>
SparseMatrix<iT,fT> SparseMatrix<iT,fT>::operator-(const SparseMatrix<iT,fT>& m) const {
	// initialise local variables
	typename std::set<std::pair<iT,iT> > sp = m.sparsityPattern();
	SparseMatrix<iT,fT> ret(*this);
	
	// perform subtraction
	for (typename std::set<std::pair<iT,iT> >::iterator i=sp.begin(); i!=sp.end(); i++)
		ret((*i).first,(*i).second,
			ret((*i).first,(*i).second) - m((*i).first,(*i).second));
	
	// return result
	return ret;
}

/**
 * multiply matrix elements by scalar
 * @param alpha scalar multiplier
 * @return result of scalar element multiplication
 */
template <typename iT, typename fT>
SparseMatrix<iT,fT> SparseMatrix<iT,fT>::operator*(fT alpha) const {
	// initalise local variables
	SparseMatrix<iT,fT> ret(*this);
    
	// perform multiplication
	ret *= alpha;
    
	// return result
	return ret;	
}

/**
 * multiply matrix by vector
 * @param v vector multiplier
 * @return result of vector multiplication
 */
template <typename iT, typename fT>
vector<fT> SparseMatrix<iT,fT>::operator*(const vector<fT> v) const {
	// initialise local variables
	fT rowsum;
	vector<fT> ret;
	
	// check dimensions
	if (Dim() == v.size())
		for (iT i=0; i<Dim(); i++) {
			rowsum = zero_;
			for (iT j=ia_.at(i); j<ia_.at(i+1); j++)
				rowsum += v.at(ja_.at(j))*a_.at(j);
			ret.push_back(rowsum);
		}
	
	// return results
	return ret;	
}

/**
 * matrix multiplication
 * @param m matrix multiplier
 * @return result of matrix multiplication
 */
template <typename iT, typename fT>
SparseMatrix<iT,fT> SparseMatrix<iT,fT>::operator*(const SparseMatrix<iT,fT>& m) const {
	// initialise local variables
	vector<fT> col; // column vector
	SparseMatrix<iT,fT> ret(Dim());
	
	// loop over columns
	for (iT j=0; j<Dim(); j++) {
		// extract column to a vector
		col.clear();
		for (iT i=0; i<Dim(); i++)
			col.push_back(static_cast<fT>(m(i,j)));
		
		// perform multiplication
		col = (*this)*col;
		
		// insert results into return matrix
		for (iT i=0; i<Dim(); i++)
			ret(i,j,col.at(i));
	}
	
	// return result
	return ret;	
}

/**
 * divide matrix elements by scalar
 * @param alpha scalar divisor
 * @return result of scalar element division
 */
template <typename iT, typename fT>
SparseMatrix<iT,fT> SparseMatrix<iT,fT>::operator/(fT alpha) const {
	// initalise local variables
	SparseMatrix<iT,fT> ret(*this);
	
	// perform multiplication
	ret /= alpha;
	
	// return result
	return ret;	
}

/**
 * matrix addition
 * @param m matrix to add to this
 * @return result of matrix addition
 */
template <typename iT, typename fT>
SparseMatrix<iT,fT>& SparseMatrix<iT,fT>::operator+=(const SparseMatrix<iT,fT>& m) {
	// initialise local variables
	typename std::set<std::pair<iT,iT> > sp = m.sparsityPattern();
	
	// perform addition
	for (typename std::set<std::pair<iT,iT> >::iterator i=sp.begin(); i!=sp.end(); i++)
		(*this)((*i).first,(*i).second,
				(*this)((*i).first,(*i).second) + m((*i).first,(*i).second));
	
	// return result
	return (*this);
}

/**
 * matrix subtraction
 * @param m matrix to subtract from this
 * @return result of matrix subtraction
 */
template <typename iT, typename fT>
SparseMatrix<iT,fT>& SparseMatrix<iT,fT>::operator-=(const SparseMatrix<iT,fT>& m) {
	// initialise local variables
	typename std::set<std::pair<iT,iT> > sp = m.sparsityPattern();
	
	// perform subtraction
	for (typename std::set<std::pair<iT,iT> >::iterator i=sp.begin(); i!=sp.end(); i++)
		(*this)((*i).first,(*i).second,
				(*this)((*i).first,(*i).second) - m((*i).first,(*i).second));
	
	// return result
	return (*this);
}

/**
 * multiply elements of this matrix by scalar
 * @param alpha scalar multiplier
 * @return result of scalar element multiplication
 */
template <typename iT, typename fT>
SparseMatrix<iT,fT>& SparseMatrix<iT,fT>::operator*=(fT alpha) {
	// loop over elements in the a vector and perform multiplication
	for (iT i=0; i<a_.size(); i++) a_.at(i) *= alpha;
	
	// return result
	return (*this);	
}

/**
 * divide elements of this matrix by scalar
 * @param alpha scalar divisor
 * @return result of scalar element division
 */
template <typename iT, typename fT>
SparseMatrix<iT,fT>& SparseMatrix<iT,fT>::operator/=(fT alpha) {
	// loop over elements in the a vector and perform division
	for (iT i=0; i<a_.size(); i++) a_.at(i) /= alpha;
	
	// return result
	return (*this);
}

// helper methods
/**
 * return solution x of Ax = b, where A is this matrix and b is vector
 * @param b right-hand-side vector
 * @return solution vector x
 */
template <typename iT, typename fT>
vector<fT> SparseMatrix<iT,fT>::solve(const vector<fT> b) const {
	// intialise local variables
	fT val, rowsum;
	vector<iT> ia, ja;
	vector<fT> a, ret = b;
	SparseMatrix<iT,fT> m(*this);
    
	// loop over columns
	for (iT col=0; col<m.Dim(); col++) {
		val = static_cast<fT>(1.)/m(col,col); // get value of diagonal
		if (fabs(val)>eps_) {
			m.rowScale(col,val);  // scale matrix row to make diagonal 1
			ret.at(col) *= val; // scale output vector
			for (iT row=(col+1); row<m.Dim(); row++) { // perform row operations
				val = m(row,col);
				m.rowOp(row,col,val);
				ret.at(row) -= val*ret.at(col);
			}
		} else {
			// OH OH! matrix is singular!  return a vector of NaNs 
			for (iT col=0; col<Dim(); col++)
				ret.at(col) = HUGE_VAL;
			return ret;
		}
	}
	
	// get storage vectors
	ia = m.IA();
	ja = m.JA();
	a  = m();
	
	// perform back substitution
	for (int i=(m.Dim()-1); i>=0; i--) {
		rowsum = zero_;
		for (iT j=(ia.at(i)+1); j<ia.at(i+1); j++) // NB: miss out diagonal element ie. ia(i)
			rowsum += ret.at(ja.at(j))*a.at(j);
		ret.at(i) -= rowsum;
	}
	
	// return results
	return ret;
}

/**
 * print out storage vectors for debugging purposes
 */
template <typename iT, typename fT>
void SparseMatrix<iT,fT>::debug(ostream& os) const {
	iT i;
	fT total = zero_;
	
	os << "ia:";
	for (i=0; i<static_cast<iT>(ia_.size()); i++)
		os << "\t" << ia_.at(i);
	os << endl;
	
	os << "ja:";
	for (i=0; i<static_cast<iT>(ja_.size()); i++)
		os << "\t" << ja_.at(i);
	os << endl;
	
	os << "a:";
	for (i=0; i<static_cast<iT>(a_.size()); i++) {
		if (a_.at(i)>=eps_) total++;
		os << "\t" << a_.at(i);
	}
	os << endl << "Sparsity = " <<
    static_cast<fT>((Dim()*Dim()-total)*100)/static_cast<fT>(Dim()*Dim()) << "\%"<< endl;
}

// protected methods
/**
 * scale matrix row with specified scalar ie. row = alpha*row
 * @param row row to mutiply
 * @param alpha scalar multiplier
 */
template <typename iT, typename fT>
void SparseMatrix<iT,fT>::rowScale(iT row, fT alpha) {
	for (iT j=ia_.at(row); j<ia_.at(row+1); j++)
		a_.at(j) *= alpha;
}

/**
 * subtract first row of matrix with scalar multiple of second row of matrix
 * ie. row_1 = row_1 - alpha*row_2 
 * @param row_1 first row
 * @param row_2 second row
 * @param alpha scalar multiplier
 */
template <typename iT, typename fT>
void SparseMatrix<iT,fT>::rowOp(iT row_1, iT row_2, fT alpha) {
	for (iT j=ia_.at(row_2); j<ia_.at(row_2+1); j++)
		(*this)(row_1,ja_.at(j),
				(*this)(row_1,ja_.at(j)) - alpha*(*this)(row_2,ja_.at(j)));
}

/**
 * return index of a_ vector associated with matrix element at specified row and column
 * @param row row of matrix element
 * @param col column of matrix element
 * @return associated value in a_ vector
 */
template <typename iT, typename fT>
iT SparseMatrix<iT,fT>::findIndex(iT row, iT col) const {
	typename vector<iT>::const_iterator jit = find(ja_.begin()+ia_.at(row),
												   ja_.begin()+ia_.at(row+1),col);
	return static_cast<iT>(distance(ja_.begin(),jit));
}

/**
 * output this matrix as a sparsity pattern
 * @return sparsity pattern (set of pairs, each pair contains row and column of non-zero element)
 */
template <typename iT, typename fT>
std::set<std::pair<iT,iT> > SparseMatrix<iT,fT>::sparsityPattern(void) const {
	// initialise local variables
	std::set<pair<iT,iT> > ret;
	
	// loop over entries
	for (iT i=0; i<Dim(); i++)
		for (iT j=ia_.at(i); j<ia_.at(i+1); j++)
			// add (row,col) pair if non-zero entry
			// NB: need the conditional because could be diagonal zero entry
			if (fabs(a_.at(j))>eps_) ret.insert(make_pair(i,ja_.at(j)));
	
	// return sparsity pattern
	return ret;
}

/**
 * create this matrix from a sparsity pattern
 * @param sp sparsity pattern (set of pairs, each pair contains row and column of non-zero element)
 */
template <typename iT, typename fT>
void SparseMatrix<iT,fT>::sparsityPattern(const std::set<std::pair<iT,iT> >& sp) {
	// initialise local variables
	typename std::set<std::pair<iT,iT> >::iterator it;
    
	// initialise data structure with zeros along diagonal
	ia_.clear();
	ja_.clear();
	a_.clear();
	for (iT i=0; i<Dim(); i++) {
		ia_.push_back(i);
		ja_.push_back(i);
		a_.push_back(zero_);
	}
	ia_.push_back(Dim());
    
	// loop over sparsity pattern inserting 1's
	for (it=sp.begin(); it!=sp.end(); it++)
		(*this)((*it).first,(*it).second,static_cast<fT>(1.0));
}

template class myg::SparseMatrix<int,double>;
template class myg::SparseMatrix<unsigned int,double>;

// global functions
/**
 * overload << to output SparseMatrix object
 */
template <typename iT, typename fT>
ostream& operator<<(ostream& os, const SparseMatrix<iT,fT>& m) {
	for (iT i=0; i<m.Dim(); i++) {
		for (iT j=0; j<m.Dim(); j++)
			os << m(i,j) << "\t";
		os << endl;
	}
	return os;
}

template ostream& operator<<(ostream& os, const SparseMatrix<int,double>& m);
template ostream& operator<<(ostream& os, const SparseMatrix<unsigned int,double>& m);
