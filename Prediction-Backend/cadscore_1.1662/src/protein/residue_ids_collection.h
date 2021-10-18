#ifndef PROTEIN_RESIDUE_IDS_COLLECTION_H_
#define PROTEIN_RESIDUE_IDS_COLLECTION_H_

#include <vector>
#include <map>
#include <string>

#include "residue_id.h"
#include "residue_summary.h"

namespace protein
{

template<typename AtomType>
std::map<ResidueID, ResidueSummary> collect_residue_ids_from_atoms(const std::vector<AtomType>& atoms)
{
	std::map<ResidueID, ResidueSummary> result;
	for(std::size_t i=0;i<atoms.size();i++)
	{
		ResidueSummary& rs=result[ResidueID::from_atom(atoms[i])];
		rs.name=atoms[i].residue_name;
		rs.atoms_count++;
	}
	return result;
}

template<typename AtomType>
std::map<ResidueID, std::vector<std::size_t> > group_atoms_indices_by_residue_ids(const std::vector<AtomType>& atoms)
{
	std::map<ResidueID, std::vector<std::size_t> > result;
	for(std::size_t i=0;i<atoms.size();i++)
	{
		result[ResidueID::from_atom(atoms[i])].push_back(i);
	}
	return result;
}

}

#endif /* PROTEIN_RESIDUE_IDS_COLLECTION_H_ */
