#ifndef CONTACTO_RESIDUE_CONTACT_AREA_DIFFERENCE_LOCAL_SCORES_H_
#define CONTACTO_RESIDUE_CONTACT_AREA_DIFFERENCE_LOCAL_SCORES_H_

#include <map>
#include <string>

#include "residue_contact_area_difference_score.h"

namespace contacto
{

template<typename ResidueID>
std::map<ResidueID, double> construct_local_scores_from_profile(
		const std::map<ResidueID, contacto::ResidueContactAreaDifferenceScore>& profile,
		const std::string& category,
		const bool get_differences,
		const bool get_references)
{
	std::map<ResidueID, double> local_scores;
	typename std::map<ResidueID, double>::iterator insertion_it=local_scores.begin();
	for(typename std::map<ResidueID, contacto::ResidueContactAreaDifferenceScore>::const_iterator it=profile.begin();it!=profile.end();++it)
	{
		const ResidueID& residue_id=it->first;
		const contacto::ResidueContactAreaDifferenceScore& residue_score=it->second;

		if(it!=profile.begin())
		{
			typename std::map<ResidueID, contacto::ResidueContactAreaDifferenceScore>::const_iterator previous_it=it;
			previous_it--;
			const ResidueID& prev_residue_id=previous_it->first;
			if((residue_id.residue_number-prev_residue_id.residue_number>1) && residue_id.chain_id==prev_residue_id.chain_id)
			{
				for(int i=prev_residue_id.residue_number+1;i<residue_id.residue_number;i++)
				{
					insertion_it=local_scores.insert(insertion_it, std::make_pair(ResidueID(prev_residue_id.chain_id, i), -1));
				}
			}
			else if(residue_id.chain_id!=prev_residue_id.chain_id)
			{
				insertion_it=local_scores.insert(insertion_it, std::make_pair(ResidueID(residue_id.chain_id, 0), -2));
			}
		}

		const contacto::Ratio ratio=residue_score.ratio(category);
		double local_score_value=-3;
		if(ratio.reference>0)
		{
			if(get_differences && !get_references)
			{
				local_score_value=ratio.difference;
			}
			else if(!get_differences && get_references)
			{
				local_score_value=ratio.reference;
			}
			else
			{
				local_score_value=(ratio.difference/ratio.reference);
			}
		}
		insertion_it=local_scores.insert(insertion_it, std::make_pair(residue_id, local_score_value));
	}
	return local_scores;
}

template<typename ResidueID>
std::map<ResidueID, double> blur_local_scores(const std::map<ResidueID, double>& local_scores, const int window_size)
{
	if(window_size==0)
	{
		return local_scores;
	}

	std::map<ResidueID, double> blured_scores;
	typename std::map<ResidueID, double>::iterator insertion_it=blured_scores.begin();
	for(typename std::map<ResidueID, double>::const_iterator it=local_scores.begin();it!=local_scores.end();++it)
	{
		const ResidueID& residue_id=it->first;
		double sum=it->second;
		int count=1;
		if(sum>=0)
		{
			typename std::map<ResidueID, double>::const_iterator left_it=it;
			double left_last=sum;

			typename std::map<ResidueID, double>::const_iterator right_it=it;
			double right_last=sum;

			for(int i=1;i<=window_size;i++)
			{
				if(left_it!=local_scores.begin())
				{
					left_it--;
					const ResidueID& left_residue_id=left_it->first;
					if(residue_id.residue_number-left_residue_id.residue_number<=window_size
							&& residue_id.chain_id==left_residue_id.chain_id
							&& left_it->second>0)
					{
						left_last=left_it->second;
					}
				}
				sum+=left_last;
				count++;

				if(right_it!=local_scores.end() && (++right_it)!=local_scores.end())
				{
					const ResidueID& right_residue_id=right_it->first;
					if(right_residue_id.residue_number-residue_id.residue_number<=window_size
							&& residue_id.chain_id==right_residue_id.chain_id
							&& right_it->second>0)
					{
						right_last=right_it->second;
					}
				}
				sum+=right_last;
				count++;
			}
		}
		insertion_it=blured_scores.insert(insertion_it, std::make_pair(residue_id, (sum/count)));
	}
	return blured_scores;
}

template<typename ResidueID>
std::map<ResidueID, double> ratios_of_local_scores(const std::map<ResidueID, double>& scores1, const std::map<ResidueID, double>& scores2)
{
	std::map<ResidueID, double> combined_scores;
	if(scores1.size()==scores2.size())
	{
		typename std::map<ResidueID, double>::iterator insertion_it=combined_scores.begin();
		typename std::map<ResidueID, double>::const_iterator it1=scores1.begin();
		typename std::map<ResidueID, double>::const_iterator it2=scores2.begin();
		while(it1!=scores1.end() && it2!=scores2.end() && it1->first==it2->first)
		{
			const ResidueID& rid=it1->first;
			double value=std::min(it1->second, it2->second);
			if(value>0.0)
			{
				value=(it1->second/it2->second);
			}
			insertion_it=combined_scores.insert(insertion_it, std::make_pair(rid, value));
			++it1;
			++it2;
		}
		if(combined_scores.size()!=scores1.size())
		{
			combined_scores.clear();
		}
	}
	return combined_scores;
}

}

#endif /* CONTACTO_RESIDUE_CONTACT_AREA_DIFFERENCE_LOCAL_SCORES_H_ */
