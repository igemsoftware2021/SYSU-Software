#ifndef PROTEIN_ATOMS_H_
#define PROTEIN_ATOMS_H_

#include <string>
#include <vector>
#include <iostream>

#include "atom.h"
#include "pdb_parsing.h"
#include "van_der_waals_radius_assigner.h"

namespace protein
{

class AtomsReading
{
public:
	static std::vector<Atom> read_atoms_from_PDB_file_stream(
			std::istream& pdb_file_stream,
			const VanDerWaalsRadiusAssigner& vdwr_assigner,
			const bool include_heteroatoms,
			const bool include_water,
			const bool include_insertions)
	{
		return collect_atoms_from_PDB_atom_records(
				read_PDB_atom_records_from_PDB_file_stream(pdb_file_stream),
				include_heteroatoms,
				include_water,
				include_insertions,
				vdwr_assigner);
	}

private:
	AtomsReading()
	{
	}

	static Atom atom_from_PDB_atom_record(const PDBAtomRecord& record, const VanDerWaalsRadiusAssigner& vdwr_assigner)
	{
		Atom atom;
		atom.chain_id=record.chain_name;
		atom.atom_number=record.atom_serial_number;
		atom.residue_number=record.residue_sequence_number;
		atom.residue_name=record.residue_name;
		atom.atom_name=record.name;
		atom.x=record.x;
		atom.y=record.y;
		atom.z=record.z;
		atom.r=vdwr_assigner.radius(atom.residue_name, atom.atom_name);
		return atom;
	}

	static std::vector<Atom> collect_atoms_from_PDB_atom_records(
			const std::vector<PDBAtomRecord>& PDB_atom_records,
			const bool include_heteroatoms,
			const bool include_water,
			const bool include_insertions,
			const VanDerWaalsRadiusAssigner& vdwr_assigner)
	{
		std::vector<Atom> atoms;
		for(std::size_t i=0;i<PDB_atom_records.size();i++)
		{
			const PDBAtomRecord& record=PDB_atom_records[i];
			if(!record.name.empty() &&
					(record.name.find("H")!=0 || (!record.element.empty() && record.element!="H")) &&
					record.name.find("1H")!=0 &&
					record.name.find("2H")!=0 &&
					record.name.find("3H")!=0 &&
					record.name.find("4H")!=0)
			{
				if(record.alternate_location_indicator.empty() || record.alternate_location_indicator=="A")
				{
					if(record.label=="ATOM" || (include_heteroatoms && record.label=="HETATM"))
					{
						if(include_water || record.residue_name!="HOH")
						{
							if(record.insertion_code.empty() || include_insertions)
							{
								atoms.push_back(atom_from_PDB_atom_record(record, vdwr_assigner));
							}
						}
					}
				}
			}
		}
		return atoms;
	}
};

}

#endif /* PROTEIN_ATOMS_H_ */
