#ifndef AUXILIARIES_STD_CONTAINERS_IO_H_
#define AUXILIARIES_STD_CONTAINERS_IO_H_

#include <string>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>
#include <set>
#include <map>

namespace auxiliaries
{

class STDContainersIO
{
public:
	static bool check_file_header(std::istream& in, const std::string& header)
	{
		std::string value;
		getline(in, value);
		while(in.good() && (value.empty() || value[0]=='#'))
		{
			value.clear();
			getline(in, value);
		}
		return (header==value);
	}

	static void print_file_comment(std::ostream& out, const std::string& comment)
	{
		out << "# " << comment << "\n";
	}

	template<typename T>
	static void print_vector(std::ostream& out, const std::string& header, const std::vector<T>& v)
	{
		if(!header.empty())
		{
			print_file_header(out, header);
		}

		out << v.size() << "\n";
		for(std::size_t i=0;i<v.size();i++)
		{
			out << v[i] << "\n";
		}
	}

	template<typename T>
	static std::vector<T> read_vector(std::istream& in, const std::string& name, const std::string& header, const bool allow_empty_result)
	{
		std::ostringstream error_output;
		if(name.empty())
		{
			error_output << "Error reading vector: ";
		}
		else
		{
			error_output << "Error reading vector '" << name << "': ";
		}

		if(!in.good())
		{
			error_output << "non-readable input stream";
			throw std::runtime_error(error_output.str());
		}

		if(!header.empty())
		{
			if(!check_file_header(in, header))
			{
				error_output << "missing file header '" << header << "'";
				throw std::runtime_error(error_output.str());
			}
		}

		std::size_t n=0;
		in >> n;

		if(in.fail())
		{
			error_output << "bad format of input stream (missing vector size)";
			throw std::runtime_error(error_output.str());
		}

		if(n>0)
		{
			std::vector<T> v(n);
			for(std::size_t i=0;i<n;i++)
			{
				in >> v[i];

				if(in.fail())
				{
					error_output << "bad format of input stream (failed to read element " << i << ")";
					throw std::runtime_error(error_output.str());
				}
			}
			return v;
		}
		else
		{
			if(!allow_empty_result)
			{
				error_output << "no data in input stream";
				throw std::runtime_error(error_output.str());
			}
			return std::vector<T>();
		}
	}

	template<typename T>
	static void print_set(std::ostream& out, const std::string& header, const std::set<T>& set)
	{
		if(!header.empty())
		{
			print_file_header(out, header);
		}

		out << set.size() << "\n";
		for(typename std::set<T>::const_iterator it=set.begin();it!=set.end();++it)
		{
			out << (*it) << "\n";
		}
	}

	template<typename T>
	static std::set<T> read_set(std::istream& in, const std::string& name, const std::string& header, const bool allow_empty_result)
	{
		std::ostringstream error_output;
		if(name.empty())
		{
			error_output << "Error reading set: ";
		}
		else
		{
			error_output << "Error reading set '" << name << "': ";
		}

		if(!in.good())
		{
			error_output << "non-readable input stream";
			throw std::runtime_error(error_output.str());
		}

		if(!header.empty())
		{
			if(!check_file_header(in, header))
			{
				error_output << "missing file header '" << header << "'";
				throw std::runtime_error(error_output.str());
			}
		}

		std::size_t n=0;
		in >> n;

		if(in.fail())
		{
			error_output << "bad format of input stream (missing set size)";
			throw std::runtime_error(error_output.str());
		}

		std::set<T> set;
		if(n>0)
		{
			typename std::set<T>::iterator prev=set.begin();
			for(std::size_t i=0;i<n;i++)
			{
				T value;
				in >> value;

				if(in.fail())
				{
					error_output << "bad format of input stream (failed to read element " << i << ")";
					throw std::runtime_error(error_output.str());
				}

				prev=set.insert(prev, value);
			}
		}
		else
		{
			if(!allow_empty_result)
			{
				error_output << "no data in input stream";
				throw std::runtime_error(error_output.str());
			}
		}
		return set;
	}

	template<typename A, typename B>
	static void print_map(std::ostream& out, const std::string& header, const std::map<A, B>& map, const bool separate_with_new_line)
	{
		if(!header.empty())
		{
			print_file_header(out, header);
		}

		out << map.size() << "\n";
		for(typename std::map<A, B>::const_iterator it=map.begin();it!=map.end();++it)
		{
			out << it->first << (separate_with_new_line ? "\n" : " ") << it->second << "\n";
		}
	}

	template<typename A, typename B>
	static std::map<A, B> read_map(std::istream& in, const std::string& name, const std::string& header, const bool allow_empty_result)
	{
		std::ostringstream error_output;
		if(name.empty())
		{
			error_output << "Error reading map: ";
		}
		else
		{
			error_output << "Error reading map '" << name << "': ";
		}

		if(!in.good())
		{
			error_output << "non-readable input stream";
			throw std::runtime_error(error_output.str());
		}

		if(!header.empty())
		{
			if(!check_file_header(in, header))
			{
				error_output << "missing file header '" << header << "'";
				throw std::runtime_error(error_output.str());
			}
		}

		std::size_t n=0;
		in >> n;

		if(in.fail())
		{
			error_output << "bad format of input stream (missing map size)";
			throw std::runtime_error(error_output.str());
		}

		std::map<A, B> map;
		if(n>0)
		{
			typename std::map<A, B>::iterator prev=map.begin();
			for(std::size_t i=0;i<n;i++)
			{
				A key;
				in >> key;

				if(in.fail())
				{
					error_output << "bad format of input stream (failed to read key of element " << i << ")";
					throw std::runtime_error(error_output.str());
				}

				B value;
				in >> value;

				if(in.fail())
				{
					error_output << "bad format of input stream (failed to read value of element " << i << ")";
					throw std::runtime_error(error_output.str());
				}

				prev=map.insert(prev, std::make_pair(key, value));
			}
		}
		else
		{
			if(!allow_empty_result)
			{
				error_output << "no data in input stream";
				throw std::runtime_error(error_output.str());
			}
		}
		return map;
	}

private:
	static void print_file_header(std::ostream& out, const std::string& header)
	{
		out << header << "\n";
	}
};

}

#endif /* AUXILIARIES_STD_CONTAINERS_IO_H_ */
