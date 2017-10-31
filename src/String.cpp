#include "stdafx.h"
#include "String.h"


void GetParValue(std::string line, std::string& PAR, std::string& VALUE) {
	size_t i = line.find('=');
	std::string par_tmp(line, 0, i);
	PAR = boost::trim_copy(par_tmp);
	if (i > 0) {
		std::string value_tmp(line, i + 1);
		VALUE = boost::trim_copy(value_tmp);
	}
}
