#include "stdafx.h"
#include <boost/serialization/array_wrapper.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/numeric/ublas/vector.hpp>

void GetParValue(std::string line, std::string& PAR, std::string& VALUE);   //Gets strings $PAR$ and $VALUE$ separated by sign $=$ from string $line$
std::vector<std::string> split_string(const std::string& str, const std::string& delimiter);