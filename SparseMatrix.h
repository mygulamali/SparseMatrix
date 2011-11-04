/**
 * @file SparseMatrix.h
 * SparseMatrix
 * Murtaza Gulamali (30/11/2006)
 *
 * This software is released under the terms and conditions of The MIT License:
 * http://www.opensource.org/licenses/mit-license.php
 */

#ifndef SPARSEMATRIX_H
#define SPARSEMATRIX_H

#include <iostream>
#include <vector>
#include <set>

namespace myg {
    
    //! SparseMatrix class
    template <typename iT, typename fT>
    class SparseMatrix {
        
    public:
        //! zero
        static const fT zero_;
        
        //! epsilon: a very small number greater than zero
        static const fT eps_;
        
        //! default constructor
        SparseMatrix(iT dim=1);
        
        //! copy constructor
        SparseMatrix(const SparseMatrix& m);
        
        //! default destructor
        ~SparseMatrix(void);
        
        // accessors and mutators for class members
        iT                      Dim(void) const;
        void                    Dim(iT dim);
        std::vector<iT>&        IA(void);
        const std::vector<iT>&  IA(void) const;
        void                    IA(const std::vector<iT>& ia);
        std::vector<iT>&        JA(void);
        const std::vector<iT>&  JA(void) const;
        void                    JA(const std::vector<iT>& ja);
        std::vector<fT>&        operator()(void);
        const std::vector<fT>&  operator()(void) const;
        void                    operator()(const std::vector<fT>& a);
        
        // accessors and mutators for matrix elements
        const fT&  operator()(iT row, iT col) const;
        void       operator()(iT row, iT col, fT val);
        
        // matrix operator methods
        SparseMatrix&   operator=(const SparseMatrix& m);
        SparseMatrix    operator+(const SparseMatrix& m) const;
        SparseMatrix    operator-(const SparseMatrix& m) const;
        SparseMatrix    operator*(fT alpha) const;
        std::vector<fT> operator*(const std::vector<fT> v) const;
        SparseMatrix    operator*(const SparseMatrix& m) const;
        SparseMatrix    operator/(fT alpha) const;
        SparseMatrix&   operator+=(const SparseMatrix& m);
        SparseMatrix&   operator-=(const SparseMatrix& m);
        SparseMatrix&   operator*=(fT alpha);
        SparseMatrix&   operator/=(fT alpha);
        
        // helper methods
        std::vector<fT> solve(const std::vector<fT> b) const;
        void            debug(std::ostream& os) const;
        
    protected:
        //! number of rows/columns (ie. dimension) of matrix
        iT nnu_;
        
        //! row index pointer vector, ia
        std::vector<iT> ia_;
        
        //! column index pointer vector, ja
        std::vector<iT> ja_;
        
        //! vector of non-zero elements and diagonal elements of matrix
        std::vector<fT> a_;
        
        // protected helper methods
        void rowScale(iT row, fT alpha);
        void rowOp(iT row_1, iT row_2, fT alpha);
        iT   findIndex(iT row, iT col) const;
        std::set<std::pair<iT,iT> > sparsityPattern(void) const;
        void                        sparsityPattern(const std::set<std::pair<iT,iT> >& sp);
        
    }; // end of SparseMatrix class
    
    //! overload operator to allow output to a stream
    template <typename iT, typename fT>
    std::ostream& operator<<(std::ostream& os,const SparseMatrix<iT,fT>& m);
    
    // inline methods for speed
    
    //! return matrix dimension
    template <typename iT, typename fT>
    inline iT SparseMatrix<iT,fT>::Dim(void) const { return nnu_; }
        
    //! set matrix dimension
    template <typename iT, typename fT>
    inline void SparseMatrix<iT,fT>::Dim(iT dim) { nnu_ = dim; }
    
    //! return row pointer vector
    template <typename iT, typename fT>
    inline std::vector<iT>& SparseMatrix<iT,fT>::IA(void) { return this->ia_; }
    
    //! return copy of row pointer vector
    template <typename iT, typename fT>
    inline const std::vector<iT>& SparseMatrix<iT,fT>::IA(void) const { return ia_; }
    
    //! set row pointer vector to specified vector 
    template <typename iT, typename fT>
    inline void SparseMatrix<iT,fT>::IA(const std::vector<iT>& ia) { ia_ = ia; }
    
    //! return column pointer vector
    template <typename iT, typename fT>
    inline std::vector<iT>& SparseMatrix<iT,fT>::JA(void) { return this->ja_; }
    
    //! return copy of column pointer vector
    template <typename iT, typename fT>
    inline const std::vector<iT>& SparseMatrix<iT,fT>::JA(void) const { return ja_; }
    
    //! set column pointer vector to specified vector
    template <typename iT, typename fT>
    inline void SparseMatrix<iT,fT>::JA(const std::vector<iT>& ja) { ja_ = ja; }
    
    //! return vector of non-zero elements and diagonal elements of matrix
    template <typename iT, typename fT>
    inline std::vector<fT>& SparseMatrix<iT,fT>::operator()(void) { return this->a_; }
    
    //! return copy of vector of non-zero elements and diagonal elements of matrix
    template <typename iT, typename fT>
    inline const std::vector<fT>& SparseMatrix<iT,fT>::operator()(void) const { return a_; }
    
    //! set vector of non-zero elements and diagonal elements of matrix to specified vector
    template <typename iT, typename fT>
    inline void SparseMatrix<iT,fT>::operator()(const std::vector<fT>& a) { a_ = a; }

} // end of namespace myg

#endif
