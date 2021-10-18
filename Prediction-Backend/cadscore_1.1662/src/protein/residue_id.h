#ifndef PROTEIN_RESIDUE_ID_H_
#define PROTEIN_RESIDUE_ID_H_

#include <string>
#include <iostream>

namespace protein
{

struct ResidueID
{
	std::string chain_id;
	int residue_number;

	ResidueID() : residue_number(0)
	{
	}

	ResidueID(const std::string& chain_id, int residue_number) : chain_id(chain_id), residue_number(residue_number)
	{
	}

	template<typename AtomType>
	static ResidueID from_atom(const AtomType& atom)
	{
		return ResidueID(atom.chain_id, atom.residue_number);
	}

	bool operator==(const ResidueID& rid) const
	{
		return (residue_number==rid.residue_number && chain_id==rid.chain_id);
	}

	bool operator!=(const ResidueID& rid) const
	{
		return (residue_number!=rid.residue_number || chain_id!=rid.chain_id);
	}

	bool operator< (const ResidueID& rid) const
	{
		return (chain_id<rid.chain_id || (chain_id==rid.chain_id && residue_number<rid.residue_number));
	}

	friend std::ostream& operator<<(std::ostream &output, const ResidueID &rid)
	{
		output << rid.chain_id << " ";
		output << rid.residue_number;
		return output;
	}

	friend std::istream& operator>>(std::istream &input, ResidueID &rid)
	{
		input >> rid.chain_id;
		input >> rid.residue_number;
		return input;
	}
};

}

#endif /* PROTEIN_RESIDUE_ID_H_ */
