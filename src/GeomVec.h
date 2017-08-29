#pragma once

#include "stdafx.h"
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>
namespace ublas = boost::numeric::ublas; // shorten the namespace name
typedef ublas::bounded_vector<double, 4> GeomVec; // define type of geometric vector
typedef ublas::bounded_matrix<double, 4, 4> GeomMat;

//Type of matrix with double
typedef  std::vector<std::vector<double>> Matrix;

//Macros which make matrix type of double size of n*m
#define CreateMatrix(name, n, m) Matrix name(n,std::vector<double>(m, 0))

double length(GeomVec x);                                         // length of the geometric vector
GeomVec x_product(GeomVec v1, GeomVec v2);                        // vector product of two vectors
GeomMat x_product_Vec_Mat(GeomVec v, GeomMat A);                  // returns x_product of Vector $v$ and Tensor $A$
GeomVec get_column(GeomMat  A, int i);                            // returns column $i$ of Tensor $A$
void    put_column(GeomMat& A, int i, GeomVec v);                 // puts Vector $v$ to Tensor $A$ column $i$
GeomMat diada(GeomVec v1, GeomVec v2);                            // returns diada of Vectors $v1$ and $v2$
GeomMat M_rotate(GeomVec& o);                                     // returns Rotation Tensor of corresponding to rotation Vector $o$
GeomVec rotate_Vector_around_vector(GeomVec v, GeomVec o);        // rotates Vector $v$ around Vector $o$ at angle |$o$|
ublas::matrix<double> InvertMat(ublas::matrix<double> A);         // returns inverted Matrix for Matrix $A$
GeomMat InitE();                                                  // initialize unit matrix
