#ifndef PROTEIN_BASIC_PARSING_H_
#define PROTEIN_BASIC_PARSING_H_

#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <stdexcept>

namespace protein
{

namespace basic_parsing
{

template<typename T>
T convert_string(const std::string& str)
{
	std::istringstream input(str);
	input.exceptions(std::istringstream::failbit | std::istringstream::badbit);
	T value;
	input >> value;
	return value;
}

template<typename T>
T safe_convert_string(const std::string& str, const T default_value)
{
	try
	{
		return convert_string<T>(str);
	}
	catch(const std::exception& e)
	{
		return default_value;
	}
}

inline std::string convert_int_to_string(const int value)
{
	std::ostringstream output;
	output << value;
	return output.str();
}

inline std::string convert_double_to_string(const double value, const int precision)
{
	std::ostringstream output;
	output << std::fixed << std::setprecision(precision) << value;
	return output.str();
}

inline std::string substring_of_columned_file_line(const std::string& line, const int start, const int end)
{
	std::string extraction;
	int line_length=static_cast<int>(line.size());
	for(int i=start-1;i<end && i<line_length && i>=0;i++)
	{
		if(line[i]!=32) { extraction.push_back(line[i]); }
	}
	return extraction;
}

inline bool insert_string_to_columned_file_line(const std::string& str, const std::size_t start, const std::size_t end, const bool shift_right, std::string& line)
{
	if(!str.empty() && start>=1 && start<=end && end<=line.size())
	{
		const std::size_t interval_length=(end-start)+1;
		if(str.size()<=interval_length)
		{
			const std::string addition(interval_length-str.size(), ' ');
			line.replace(start-1, interval_length, (shift_right ? addition+str : str+addition));
		}
	}
	return false;
}

}

}

#endif /* PROTEIN_BASIC_PARSING_H_ */
