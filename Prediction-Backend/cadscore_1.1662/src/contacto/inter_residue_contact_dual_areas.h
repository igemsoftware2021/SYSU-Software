#ifndef CONTACTO_INTER_RESIDUE_CONTACT_DUAL_AREAS_H_
#define CONTACTO_INTER_RESIDUE_CONTACT_DUAL_AREAS_H_

#include <map>
#include <string>
#include <iostream>

namespace contacto
{

struct InterResidueContactDualAreas
{
	typedef std::map<std::string, std::pair<double, double> > AreasMap;

	AreasMap areas;

	InterResidueContactDualAreas()
	{
	}

	std::pair<double, double> area(const std::string& area_class) const
	{
		AreasMap::const_iterator it=areas.find(area_class);
		return (it!=areas.end() ? it->second : std::make_pair(0.0, 0.0));
	}

	friend std::ostream& operator<<(std::ostream& output, const InterResidueContactDualAreas& contact)
	{
		output << contact.areas.size() << "\n";
		for(AreasMap::const_iterator it=contact.areas.begin();it!=contact.areas.end();++it)
		{
			output << it->first << " " << it->second.first << " " << it->second.second << "\n";
		}
		return output;
	}

	friend std::istream& operator>>(std::istream& input, InterResidueContactDualAreas& contact)
	{
		contact.areas.clear();
		std::size_t n=0;
		input >> n;
		for(std::size_t i=0;i<n;i++)
		{
			std::string key;
			input >> key;
			double value1=0;
			input >> value1;
			double value2=0;
			input >> value2;
			contact.areas[key]=std::make_pair(value1, value2);
		}
		return input;
	}
};

}

#endif /* CONTACTO_INTER_RESIDUE_CONTACT_DUAL_AREAS_H_ */
