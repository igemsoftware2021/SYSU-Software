#ifndef PROTEIN_ATOM_ID_H_
#define PROTEIN_ATOM_ID_H_

#include <string>
#include <iostream>

#include "residue_id.h"

namespace protein
{

struct AtomID
{
	ResidueID residue_id;
	std::string atom_name;

	AtomID()
	{
	}

	AtomID(const ResidueID& residue_id, const std::string& atom_name) : residue_id(residue_id), atom_name(canonical_atom_name(atom_name))
	{
	}

	template<typename AtomType>
	static AtomID from_atom(const AtomType& atom)
	{
		return AtomID(ResidueID::from_atom(atom), atom.atom_name);
	}

	bool operator==(const AtomID& aid) const
	{
		return (residue_id==aid.residue_id && atom_name==aid.atom_name);
	}

	bool operator!=(const AtomID& aid) const
	{
		return (residue_id!=aid.residue_id || atom_name!=aid.atom_name);
	}

	bool operator< (const AtomID& aid) const
	{
		return (residue_id<aid.residue_id || (residue_id==aid.residue_id && atom_name<aid.atom_name));
	}

	friend std::ostream& operator<<(std::ostream &output, const AtomID &aid)
	{
		output << aid.residue_id << " ";
		output << aid.atom_name;
		return output;
	}

	friend std::istream& operator>>(std::istream &input, AtomID &aid)
	{
		input >> aid.residue_id;
		input >> aid.atom_name;
		return input;
	}

	static std::string canonical_atom_name(const std::string& atom_name)
	{
		std::string fixed_atom_name=atom_name;
		if(fixed_atom_name.size()==2 && fixed_atom_name[1]=='\'') { fixed_atom_name[1]='*'; }
		else if(fixed_atom_name.size()==3 && fixed_atom_name[2]=='\'') { fixed_atom_name[2]='*'; }
		else if(fixed_atom_name=="O1P") { fixed_atom_name="OP1"; }
		else if(fixed_atom_name=="O2P") { fixed_atom_name="OP2"; }
		else if(fixed_atom_name=="O3P") { fixed_atom_name="OP3"; }
		return fixed_atom_name;
	}
};

}

#endif /* PROTEIN_ATOM_ID_H_ */
