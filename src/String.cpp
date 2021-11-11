#include "stdafx.h"
#include "String.h"


void GetParValue(std::string line, std::string& PAR, std::string& VALUE) {
	size_t i = line.find('=');
	std::string par_tmp(line, 0, i);
	PAR = boost::trim_copy(par_tmp);
	if (i > 0 && i != std::string::npos) {
		std::string value_tmp(line, i + 1);
		VALUE = boost::trim_copy(value_tmp);
	}
}


std::vector<std::string> split_string(const std::string& str, const std::string& delimiter)
{
	std::vector<std::string> strings;

	std::string::size_type pos = 0;
	std::string::size_type prev = 0;
	while ((pos = str.find(delimiter, prev)) != std::string::npos)
	{
		strings.push_back(str.substr(prev, pos - prev));
		prev = pos + 1;
	}

	// To get the last substring (or only, if delimiter is not found)
	strings.push_back(str.substr(prev));

	return strings;
}