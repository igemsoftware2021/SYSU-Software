#ifndef ATOMS_CLASSIFICATION_H_
#define ATOMS_CLASSIFICATION_H_

#include <vector>
#include <set>
#include <map>

#include "atom.h"
#include "residue_ids_collection.h"

namespace protein
{

class AtomsClassification
{
public:
	static void classify_atoms(std::vector<Atom>& atoms)
	{
		const std::set<std::string> amino_acid_names=construct_amino_acid_names_set();
		const std::set<std::string> nucleotide_names=construct_nucleotide_names_set();
		const std::set<std::string> amino_acid_main_chain_atom_names=construct_amino_acid_main_chain_atom_names_set();
		const std::set<std::string> nucleotide_main_chain_atom_names=construct_nucleotide_main_chain_atom_names_set();
		const std::map<ResidueID, std::vector<std::size_t> > grouped_atoms_indices=group_atoms_indices_by_residue_ids(atoms);
		for(std::map<ResidueID, std::vector<std::size_t> >::const_iterator it=grouped_atoms_indices.begin();it!=grouped_atoms_indices.end();++it)
		{
			const std::vector<std::size_t>& indices=it->second;
			if(!indices.empty())
			{
				int molecule_class=Atom::unidentified_molecule_class;

				const Atom& first_atom=atoms[indices.front()];
				if(amino_acid_names.count(first_atom.residue_name)>0)
				{
					molecule_class=Atom::amino_acid;
				}
				else if(nucleotide_names.count(first_atom.residue_name)>0)
				{
					molecule_class=Atom::nucleotide;
				}
				else
				{
					std::size_t amino_acid_core_match_count=0;
					std::size_t nucleotide_core_match_count=0;
					for(std::size_t i=0;i<indices.size();i++)
					{
						const Atom& atom=atoms[indices[i]];
						if(amino_acid_main_chain_atom_names.count(atom.atom_name)>0)
						{
							amino_acid_core_match_count++;
						}
						if(nucleotide_main_chain_atom_names.count(atom.atom_name)>0)
						{
							nucleotide_core_match_count++;
						}
					}
					if(amino_acid_core_match_count!=nucleotide_core_match_count)
					{
						if(amino_acid_core_match_count>=AMINO_ACID_MAIN_CHAIN_ATOMS_MIN_COUNT)
						{
							molecule_class=Atom::amino_acid;
						}
						else if(nucleotide_core_match_count>=NUCLEOTIDE_MAIN_CHAIN_ATOMS_MIN_COUNT)
						{
							molecule_class=Atom::nucleotide;
						}
					}
				}

				for(std::size_t i=0;i<indices.size();i++)
				{
					Atom& atom=atoms[indices[i]];
					atom.molecule_class=molecule_class;
					atom.location_class=Atom::unidentified_location_class;
					if(molecule_class==static_cast<int>(Atom::amino_acid))
					{
						atom.location_class=amino_acid_main_chain_atom_names.count(atom.atom_name)>0 ? Atom::main_chain : Atom::side_chain;
					}
					else if(molecule_class==static_cast<int>(Atom::nucleotide))
					{
						atom.location_class=nucleotide_main_chain_atom_names.count(atom.atom_name)>0 ? Atom::main_chain : Atom::side_chain;
					}
				}
			}
		}
	}

private:
	static const std::size_t AMINO_ACID_MAIN_CHAIN_ATOMS_MIN_COUNT=3;
	static const std::size_t NUCLEOTIDE_MAIN_CHAIN_ATOMS_MIN_COUNT=11;

	AtomsClassification()
	{
	}

	static std::set<std::string> construct_amino_acid_names_set()
	{
		std::set<std::string> names;
		names.insert("ALA");
		names.insert("ARG");
		names.insert("ASN");
		names.insert("ASP");
		names.insert("CYS");
		names.insert("GLU");
		names.insert("GLN");
		names.insert("GLY");
		names.insert("HIS");
		names.insert("ILE");
		names.insert("LEU");
		names.insert("LYS");
		names.insert("MET");
		names.insert("PHE");
		names.insert("PRO");
		names.insert("SER");
		names.insert("THR");
		names.insert("TRP");
		names.insert("TYR");
		names.insert("VAL");
		return names;
	}

	static std::set<std::string> construct_nucleotide_names_set()
	{
		std::set<std::string> names;
		names.insert("A");
		names.insert("C");
		names.insert("G");
		names.insert("T");
		names.insert("U");
		return names;
	}

	static std::set<std::string> construct_amino_acid_main_chain_atom_names_set()
	{
		std::set<std::string> names;

		names.insert("C");
		names.insert("CA");
		names.insert("N");

		names.insert("O");
		names.insert("OXT");

		return names;
	}

	static std::set<std::string> construct_nucleotide_main_chain_atom_names_set()
	{
		std::set<std::string> names;

		names.insert("OP3");
		names.insert("O3P");

		names.insert("P");

		names.insert("OP1");
		names.insert("O1P");

		names.insert("OP2");
		names.insert("O2P");

		names.insert("O5'");
		names.insert("O5*");

		names.insert("C5'");
		names.insert("C5*");

		names.insert("C4'");
		names.insert("C4*");

		names.insert("O4'");
		names.insert("O4*");

		names.insert("C3'");
		names.insert("C3*");

		names.insert("O3'");
		names.insert("O3*");

		names.insert("C2'");
		names.insert("C2*");

		names.insert("O2'");
		names.insert("O2*");

		names.insert("C1'");
		names.insert("C1*");
		return names;
	}
};

}

#endif /* ATOMS_CLASSIFICATION_H_ */
