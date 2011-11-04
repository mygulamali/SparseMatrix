# SparseMatrix

## Synopsis

SparseMatrix is a small C++ class to represent a sparse square matrix in
compressed row format.  It's very useful for FEM codes when coupled with a
decent matrix solver (eg. Algebraic [Multigrid Method](http://en.wikipedia.org/wiki/Multigrid_method "Multigrid Method")).

## Explanation

The SparseMatrix class stores matrix entries in compressed row format.  That is,
each each diagonal element is stored together with any non-zero off-diagonal
elements.  This makes it very memory efficient for storing large sparse matrices
such as those used in the [finite element method](http://en.wikipedia.org/wiki/Finite_element_method "Finite Element Method").

I wrote this class as part of a job interview and thought I'd add open source it
so that others might benefit.  Much of the API is commented in [Doxygen](http://www.doxygen.org "Doxygen")
format, and the class is templated for flexibility.  A very basic solve method
is also included for testing purposes.

Hope it helps!

## License

This software is released under the terms and conditions of [The MIT
License](http://www.opensource.org/licenses/mit-license.php "The MIT License").
Please see the license.txt file for more details.
