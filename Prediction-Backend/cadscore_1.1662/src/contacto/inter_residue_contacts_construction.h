#ifndef CONTACTO_INTER_RESIDUE_CONTACTS_CONSTRUCTION_H_
#define CONTACTO_INTER_RESIDUE_CONTACTS_CONSTRUCTION_H_

#include <vector>
#include <map>

#include "inter_atom_contact.h"
#include "contact_classification.h"
#include "contact_id.h"
#include "inter_residue_contact_areas.h"

namespace contacto
{

template<typename Atom, typename ResidueID>
std::map< ContactID<ResidueID>, InterResidueContactAreas > construct_inter_residue_contacts(const std::vector<Atom>& atoms, const std::vector<InterAtomContact>& inter_atom_contacts)
{
	std::map< ContactID<ResidueID>, InterResidueContactAreas > inter_residue_contacts_map;
	for(std::size_t i=0;i<inter_atom_contacts.size();i++)
	{
		const InterAtomContact& inter_atom_contact=inter_atom_contacts[i];
		const Atom& a1=atoms[inter_atom_contact.a];
		const Atom& a2=atoms[inter_atom_contact.b];
		const std::vector<std::string> contact_classes=ContactClassification::classify_atoms_contact<Atom, ResidueID>(a1, a2);
		if(!contact_classes.empty())
		{
			InterResidueContactAreas& inter_residue_contact_areas=inter_residue_contacts_map[ContactID<ResidueID>(ResidueID::from_atom(a1), ResidueID::from_atom(a2))];
			for(std::size_t j=0;j<contact_classes.size();j++)
			{
				inter_residue_contact_areas.areas[contact_classes[j]]+=inter_atom_contact.area;
			}
		}
	}
	return inter_residue_contacts_map;
}

}

#endif /* CONTACTO_INTER_RESIDUE_CONTACTS_CONSTRUCTION_H_ */
