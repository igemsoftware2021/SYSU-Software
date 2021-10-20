#ifndef PROTEIN_VANDERWAALSRADIUSASSIGNER_H_
#define PROTEIN_VANDERWAALSRADIUSASSIGNER_H_

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <stdexcept>

namespace protein
{

class VanDerWaalsRadiusAssigner
{
public:
	VanDerWaalsRadiusAssigner(std::istream& classes_stream, std::istream& members_stream) :
		classes_(read_map<std::string, double>(classes_stream)),
		members_(read_map<std::string, std::string>(members_stream))
	{
		if(classes_.empty())
		{
			throw std::runtime_error("Classes map is empty");
		}
		if(members_.empty())
		{
			throw std::runtime_error("Members map is empty");
		}
	}

	double radius(const std::string& residue, const std::string& atom) const
	{
		const std::string sep="_";
		const std::string any="*";

		std::vector<std::string> residue_patterns;
		residue_patterns.push_back(residue);
		residue_patterns.push_back(any);

		std::vector<std::string> atom_patterns;
		atom_patterns.push_back(atom);
		for(std::size_t i=0;i<atom.size();i++)
		{
			atom_patterns.push_back(atom.substr(0, atom.size()-i)+any);
		}
		atom_patterns.push_back(any);

		for(std::size_t i=0;i<residue_patterns.size();i++)
		{
			for(std::size_t j=0;j<atom_patterns.size();j++)
			{
				const std::string member_pattern=residue_patterns[i]+sep+atom_patterns[j];
				std::map<std::string, std::string>::const_iterator member_iterator=members_.find(member_pattern);
				if(member_iterator!=members_.end())
				{
					const std::string& class_name=member_iterator->second;
					std::map<std::string, double>::const_iterator class_iterator=classes_.find(class_name);
					if(class_iterator!=classes_.end())
					{
						return class_iterator->second;
					}
					else
					{
						throw std::runtime_error(std::string("Missing class ")+class_name+std::string(" for member ")+member_pattern);
					}
				}
			}
		}

		throw std::runtime_error(std::string("Missing member for ")+(residue+sep+atom));

		return 0;
	}

private:
	template<typename A, typename B>
	static std::map<A, B> read_map(std::istream& input)
	{
		std::map<A, B> map;
		while(input.good())
		{
			A key;
			B value;
			input >> key;
			if(!input.fail())
			{
				input >> value;
				if(!input.fail())
				{
					map[key]=value;
				}
			}
		}
		return map;
	}

	std::map<std::string, double> classes_;
	std::map<std::string, std::string> members_;
};

}

#endif /* PROTEIN_VANDERWAALSRADIUSASSIGNER_H_ */
