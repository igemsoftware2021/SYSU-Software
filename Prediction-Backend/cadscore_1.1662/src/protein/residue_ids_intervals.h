#ifndef RESIDUE_IDS_INTERVALS_H_
#define RESIDUE_IDS_INTERVALS_H_

#include <sstream>
#include <vector>
#include <limits>

#include "residue_id.h"

namespace protein
{

class ResidueIDsIntervalsReader
{
public:
	static bool read_residue_ids_intervals(const std::string& input, std::vector< std::vector< std::pair<ResidueID, ResidueID> > >& result)
	{
		std::vector< std::vector< std::pair<ResidueID, ResidueID> > > intervals;
		std::string::size_type a=input.find('(', 0);
		std::string::size_type b=input.find(')', 0);
		while(a<input.size() && b<input.size())
		{
			a++;
			if(a>=b)
			{
				return false;
			}
			else
			{
				const std::string block=input.substr(a, b-a);
				std::vector< std::pair<ResidueID, ResidueID> > block_intervals;
				if(!read_residue_ids_intervals(block, block_intervals))
				{
					return false;
				}
				else
				{
					intervals.push_back(block_intervals);
					a=input.find('(', b+1);
					b=input.find(')', b+1);
				}
			}
		}
		if(intervals.empty())
		{
			return false;
		}
		else
		{
			result=intervals;
			return true;
		}
	}

private:
	static bool read_residue_id(const std::string& input, ResidueID& result)
	{
		std::string filtered_input;
		filtered_input.reserve(input.size());
		for(std::string::size_type i=0;i<input.size();i++)
		{
			const char v=input[i];
			if(v!=' ' && v!='\t' && v!='\n' && v!='\r')
			{
				filtered_input.push_back(v);
			}
		}

		if(filtered_input.empty())
		{
			return false;
		}
		else
		{
			std::string chain_id_str;
			std::string residue_number_str;
			if(filtered_input.substr(0,1).find_first_of("0123456789")==0)
			{
				chain_id_str="?";
				residue_number_str=filtered_input;
			}
			else
			{
				chain_id_str=filtered_input.substr(0,1);
				if(filtered_input.size()>1)
				{
					residue_number_str=filtered_input.substr(1);
				}
			}

			int residue_number=std::numeric_limits<int>::min();
			if(!residue_number_str.empty())
			{
				for(std::string::size_type i=0;i<residue_number_str.size();i++)
				{
					if(std::string("0123456789").find(residue_number_str[i], 0)==std::string::npos)
					{
						return false;
					}
				}
				std::istringstream num_stream(residue_number_str);
				num_stream >> residue_number;
				if(num_stream.fail())
				{
					return false;
				}
			}

			result.chain_id=chain_id_str;
			result.residue_number=residue_number;
			return true;
		}
	}

	static bool read_residue_ids_interval(const std::string& input, std::pair<ResidueID, ResidueID>& result)
	{
		std::string::size_type a=input.find('-', 0);
		if(a<input.size())
		{
			ResidueID rid1;
			ResidueID rid2;
			if(read_residue_id(input.substr(0, a), rid1) && read_residue_id(input.substr(a+1), rid2))
			{
				if(rid1.chain_id!=rid2.chain_id || rid1.residue_number==std::numeric_limits<int>::min() || rid2.residue_number==std::numeric_limits<int>::min())
				{
					return false;
				}
				else
				{
					result=std::make_pair(rid1, rid2);
					return true;
				}
			}
			else
			{
				return false;
			}
		}
		else
		{
			ResidueID rid;
			if(read_residue_id(input, rid))
			{
				if(rid.residue_number==std::numeric_limits<int>::min())
				{
					result=std::make_pair(rid, ResidueID(rid.chain_id, std::numeric_limits<int>::max()));
				}
				else
				{
					result=std::make_pair(rid, rid);
				}
				return true;
			}
			else
			{
				return false;
			}
		}
	}

	static bool read_residue_ids_intervals(const std::string& input, std::vector< std::pair<ResidueID, ResidueID> >& result)
	{
		std::vector< std::pair<ResidueID, ResidueID> > intervals;
		std::string::size_type a=0;
		std::string::size_type b=input.find(',', 0);
		while(a<input.size())
		{
			if(a>=b)
			{
				return false;
			}
			else
			{
				const std::string block=input.substr(a, b-a);
				std::pair<ResidueID, ResidueID> block_interval;
				if(read_residue_ids_interval(block, block_interval))
				{
					intervals.push_back(block_interval);
					if(b<input.size())
					{
						a=b+1;
						b=input.find(',', b+1);
					}
					else
					{
						a=b;
					}
				}
				else
				{
					return false;
				}
			}
		}
		if(intervals.empty())
		{
			return false;
		}
		else
		{
			result=intervals;
			return true;
		}
	}
};

}

#endif /* RESIDUE_IDS_INTERVALS_H_ */
