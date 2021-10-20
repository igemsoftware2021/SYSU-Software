#ifndef CONTACTO_INTER_RESIDUE_CONTACTS_COMBINATION_H_
#define CONTACTO_INTER_RESIDUE_CONTACTS_COMBINATION_H_

#include <map>

#include "contact_id.h"
#include "inter_residue_contact_areas.h"
#include "inter_residue_contact_dual_areas.h"
#include "utilities.h"

namespace contacto
{

inline InterResidueContactDualAreas combine_two_contacts(const InterResidueContactAreas& contact_1, const InterResidueContactAreas& contact_2, const bool binarize)
{
	InterResidueContactDualAreas combination;
	combination.areas=combine_two_maps(contact_1.areas, contact_2.areas);
	if(binarize)
	{
		for(InterResidueContactDualAreas::AreasMap::iterator it=combination.areas.begin();it!=combination.areas.end();++it)
		{
			it->second.first=((it->second.first>0.0) ? 1.0 : 0.0);
			it->second.second=((it->second.second>0.0) ? 1.0 : 0.0);
		}
	}
	return combination;
}

template<typename ResidueID>
std::map< ContactID<ResidueID>, InterResidueContactDualAreas > combine_two_inter_residue_contact_maps(const std::map< ContactID<ResidueID>, InterResidueContactAreas >& contact_map_1, const std::map< ContactID<ResidueID>, InterResidueContactAreas >& contact_map_2, const bool binarize)
{
	typedef std::map< ContactID<ResidueID>, std::pair<InterResidueContactAreas, InterResidueContactAreas> > Combination;
	const Combination combination=combine_two_maps(contact_map_1, contact_map_2);
	typedef std::map< ContactID<ResidueID>, InterResidueContactDualAreas > Result;
	Result result;
	typename Result::iterator prev=result.begin();
	for(typename Combination::const_iterator it=combination.begin();it!=combination.end();++it)
	{
		prev=result.insert(prev, std::make_pair(it->first, combine_two_contacts(it->second.first, it->second.second, binarize)));
	}
	return result;
}

}

#endif /* CONTACTO_INTER_RESIDUE_CONTACTS_COMBINATION_H_ */
