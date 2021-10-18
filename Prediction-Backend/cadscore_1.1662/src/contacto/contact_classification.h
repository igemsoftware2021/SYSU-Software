#ifndef CONTACTO_CONTACT_CLASSIFICATION_H_
#define CONTACTO_CONTACT_CLASSIFICATION_H_

#include <string>
#include <set>
#include <vector>
#include <cstdlib>

namespace contacto
{

class ContactClassification
{
public:
	template<typename AtomType, typename ResidueIDType>
	static std::vector<std::string> classify_atoms_contact(const AtomType& a, const AtomType& b)
	{
		if(a==b)
		{
			std::vector<std::string> atom_classes=classify_atom_by_name(a);
			for(std::size_t i=0;i<atom_classes.size();i++)
			{
				atom_classes[i]+="W";
			}
			return atom_classes;
		}
		else if(ResidueIDType::from_atom(a)!=ResidueIDType::from_atom(b))
		{
			if(!check_if_atoms_contact_is_covalent(a, b))
			{
				return classify_atoms_contact_by_names(a, b);
			}
		}
		return std::vector<std::string>();
	}

	static std::vector<std::string> get_all_classes_list()
	{
		std::vector<std::string> all_classes;
		std::string c1[3]={"A", "M", "S"};
		std::string c2[4]={"A", "M", "S", "W"};
		for(int i=0;i<3;i++)
		{
			for(int j=0;j<4;j++)
			{
				all_classes.push_back(c1[i]+c2[j]);
			}
		}
		return all_classes;
	}

private:
	ContactClassification()
	{
	}

	template<typename AtomType>
	static std::vector<std::string> classify_atom_by_name(const AtomType& atom)
	{
		std::vector<std::string> classes;
		if(atom.location_class==static_cast<int>(AtomType::main_chain))
		{
			classes.push_back("M");
		}
		else if(atom.location_class==static_cast<int>(AtomType::side_chain))
		{
			classes.push_back("S");
		}
		classes.push_back("A");
		return classes;
	}

	template<typename AtomType>
	static std::vector<std::string> classify_atoms_contact_by_names(const AtomType& a, const AtomType& b)
	{
		std::vector<std::string> classes;
		const std::vector<std::string> a_classes=classify_atom_by_name(a);
		const std::vector<std::string> b_classes=classify_atom_by_name(b);
		for(std::size_t i=0;i<a_classes.size();i++)
		{
			for(std::size_t j=0;j<b_classes.size();j++)
			{
				classes.push_back(a_classes[i]+b_classes[j]);
			}
		}
		return classes;
	}

	template<typename AtomType>
	static bool check_if_atoms_contact_is_covalent(const AtomType& a, const AtomType& b)
	{
		if(a.molecule_class==b.molecule_class && abs(a.residue_number-b.residue_number)==1 && a.chain_id==b.chain_id)
		{
			const std::string summary=a.atom_name+"-"+b.atom_name;
			if(a.molecule_class==static_cast<int>(AtomType::amino_acid))
			{
				return (summary=="C-N" || summary=="N-C");
			}
			else if(a.molecule_class==static_cast<int>(AtomType::nucleotide))
			{
				return (summary=="P-O3'" || summary=="P-O3*" || summary=="O3'-P" || summary=="O3*-P");
			}
		}
		return false;
	}
};

}

#endif /* CONTACTO_CONTACT_CLASSIFICATION_H_ */
