#ifndef CONTACTO_RESIDUE_CONTACT_AREA_DIFFERENCE_PROFILE_H_
#define CONTACTO_RESIDUE_CONTACT_AREA_DIFFERENCE_PROFILE_H_

#include <map>
#include <set>

#include "contact_id.h"
#include "inter_residue_contact_dual_areas.h"
#include "residue_contact_area_difference_score.h"

namespace contacto
{

template<typename ResidueID, typename ResidueSummary, typename DifferenceProducer, typename ReferenceProducer>
std::map<ResidueID, ResidueContactAreaDifferenceScore> construct_residue_contact_area_difference_profile(
		const std::map<ContactID<ResidueID>, InterResidueContactDualAreas>& combined_inter_residue_contacts,
		const std::map<ResidueID, ResidueSummary>& residue_ids)
{
	std::map<ResidueID, ResidueContactAreaDifferenceScore> profile;
	{
		typename std::map<ResidueID, ResidueContactAreaDifferenceScore>::iterator prev=profile.begin();
		for(typename std::map<ResidueID, ResidueSummary>::const_iterator it=residue_ids.begin();it!=residue_ids.end();++it)
		{
			prev=profile.insert(prev, std::make_pair(it->first, ResidueContactAreaDifferenceScore()));
		}
	}
	DifferenceProducer difference_producer;
	ReferenceProducer reference_producer;
	for(typename std::map<ContactID<ResidueID>, InterResidueContactDualAreas>::const_iterator it=combined_inter_residue_contacts.begin();it!=combined_inter_residue_contacts.end();++it)
	{
		const ResidueID& residue_id=it->first.a;
		typename std::map<ResidueID, ResidueContactAreaDifferenceScore>::iterator profile_it=profile.find(residue_id);
		if(profile_it!=profile.end())
		{
			ResidueContactAreaDifferenceScore& residue_score=profile_it->second;
			const InterResidueContactDualAreas::AreasMap& areas_map=it->second.areas;
			for(InterResidueContactDualAreas::AreasMap::const_iterator jt=areas_map.begin();jt!=areas_map.end();++jt)
			{
				Ratio& ratio=residue_score.ratios[jt->first];
				ratio.difference+=difference_producer(jt->second.first, jt->second.second);
				ratio.reference+=reference_producer(jt->second.first, jt->second.second);
			}
		}
	}
	return profile;
}

template<typename ResidueID>
ResidueContactAreaDifferenceScore calculate_global_contact_area_difference_score_from_profile(
		const std::map<ResidueID, ResidueContactAreaDifferenceScore>& profile,
		const bool use_min)
{
	ResidueContactAreaDifferenceScore global_score;
	for(typename std::map<ResidueID, ResidueContactAreaDifferenceScore>::const_iterator it=profile.begin();it!=profile.end();++it)
	{
		const ResidueContactAreaDifferenceScore& residue_score=it->second;
		for(ResidueContactAreaDifferenceScore::RatiosMap::const_iterator jt=residue_score.ratios.begin();jt!=residue_score.ratios.end();++jt)
		{
			if(jt->second.reference>0.0)
			{
				Ratio& ratio=global_score.ratios[jt->first];
				ratio.difference+=(use_min ? std::min(jt->second.difference, jt->second.reference) : jt->second.difference);
				ratio.reference+=jt->second.reference;
			}
		}
	}
	return global_score;
}

}

#endif /* CONTACTO_RESIDUE_CONTACT_AREA_DIFFERENCE_PROFILE_H_ */
