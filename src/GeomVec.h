#pragma once

#include "stdafx.h"
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
namespace ublas = boost::numeric::ublas; // shorten the namespace name
typedef ublas::bounded_vector<double, DIM + 1> GeomVec; // define type of geometric vector

double length(GeomVec x); // length of the geometric vector