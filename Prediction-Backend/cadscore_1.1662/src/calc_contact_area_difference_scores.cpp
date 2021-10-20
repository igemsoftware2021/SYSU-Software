#include <set>
#include <sstream>

#include "protein/residue_id.h"
#include "protein/residue_summary.h"

#include "contacto/residue_contact_area_difference_local_scores.h"
#include "contacto/residue_contact_area_difference_profile.h"
#include "contacto/contact_classification.h"

#include "auxiliaries/command_line_options.h"
#include "auxiliaries/std_containers_io.h"

void calc_contact_area_difference_local_scores(const auxiliaries::CommandLineOptions& clo)
{
	clo.check_allowed_options("--category: --window: --absolute");

	const std::string category=clo.arg<std::string>("--category");
	const int window_size=clo.arg_in_interval<int>("--window", 0, 1000);
	const bool absolute=clo.isopt("--absolute");

	const std::map<protein::ResidueID, contacto::ResidueContactAreaDifferenceScore> profile=auxiliaries::STDContainersIO::read_map<protein::ResidueID, contacto::ResidueContactAreaDifferenceScore>(std::cin, " CAD profile", "cad_profile", false);

	std::map<protein::ResidueID, double> local_scores=contacto::blur_local_scores(contacto::construct_local_scores_from_profile(profile, category, true, false), window_size);
	if(!absolute)
	{
		local_scores=contacto::ratios_of_local_scores(local_scores, contacto::blur_local_scores(contacto::construct_local_scores_from_profile(profile, category, false, true), window_size));
	}

	if(local_scores.empty())
	{
		throw std::runtime_error("No local scores calculated");
	}
	else
	{
		auxiliaries::STDContainersIO::print_file_comment(std::cout, "This file contains contact area differences for each residue.");
		auxiliaries::STDContainersIO::print_file_comment(std::cout, "");
		auxiliaries::STDContainersIO::print_file_comment(std::cout, "Each line contains residue chain name, residue number and the difference value.");
		if(!absolute)
		{
			auxiliaries::STDContainersIO::print_file_comment(std::cout, "Each non-negative difference value is in the range [0, 1]. 0 means no difference.");
		}
		auxiliaries::STDContainersIO::print_file_comment(std::cout, "Negative difference values mean that contacts are undefined.");
		auxiliaries::STDContainersIO::print_file_comment(std::cout, "");
		auxiliaries::STDContainersIO::print_file_comment(std::cout, category+" contact category was used");
		if(window_size==0)
		{
			auxiliaries::STDContainersIO::print_file_comment(std::cout, "No smoothing window was used.");
		}
		else
		{
			std::ostringstream outstr;
			outstr << "Smoothing window " << window_size << " was used.";
			auxiliaries::STDContainersIO::print_file_comment(std::cout, outstr.str());
		}
		std::cout << "\n";

		auxiliaries::STDContainersIO::print_map(std::cout, "local_scores", local_scores, false);
	}
}

void calc_contact_area_difference_global_scores(const auxiliaries::CommandLineOptions& clo)
{
	clo.check_allowed_options("--use-min --categories: --add-categories:");

	const bool use_min=clo.isopt("--use-min");

	const std::vector<std::string> main_categories=clo.isopt("--categories") ? clo.arg_vector<std::string>("--categories", ',') : contacto::ContactClassification::get_all_classes_list();
	const std::vector<std::string> additional_categories=clo.isopt("--add-categories") ? clo.arg_vector<std::string>("--add-categories", ',') : std::vector<std::string>();

	std::vector<std::string> all_categories=main_categories;
	if(!additional_categories.empty())
	{
		all_categories.insert(all_categories.end(), additional_categories.begin(), additional_categories.end());
	}

	const std::map<protein::ResidueID, contacto::ResidueContactAreaDifferenceScore> residue_contact_area_difference_profile=auxiliaries::STDContainersIO::read_map<protein::ResidueID, contacto::ResidueContactAreaDifferenceScore>(std::cin, " CAD profile", "cad_profile", false);

	const contacto::ResidueContactAreaDifferenceScore global_score=contacto::calculate_global_contact_area_difference_score_from_profile(residue_contact_area_difference_profile, use_min);

	for(std::size_t i=0;i<all_categories.size();i++)
	{
		const std::string& the_class=all_categories[i];
		const contacto::Ratio ratio=global_score.ratio(the_class);
		std::cout << the_class << "_diff " << ratio.difference << "\n";
		std::cout << the_class << "_ref " << ratio.reference << "\n";
		std::cout << the_class << " " << (ratio.reference>0.0 ? (1-(ratio.difference/ratio.reference)) : 0.0) << "\n";
	}
}

void calc_contact_area_difference_size_scores(const auxiliaries::CommandLineOptions& clo)
{
	clo.check_allowed_options("");

	const std::map<protein::ResidueID, contacto::ResidueContactAreaDifferenceScore> profile=auxiliaries::STDContainersIO::read_map<protein::ResidueID, contacto::ResidueContactAreaDifferenceScore>(std::cin, "CAD profile", "cad_profile", false);

	const std::map<protein::ResidueID, protein::ResidueSummary> residue_ids_of_target=auxiliaries::STDContainersIO::read_map<protein::ResidueID, protein::ResidueSummary>(std::cin, "target residue identifiers", "residue_ids", false);

	const std::map<protein::ResidueID, protein::ResidueSummary> residue_ids_of_model=auxiliaries::STDContainersIO::read_map<protein::ResidueID, protein::ResidueSummary>(std::cin, "model residue identifiers", "residue_ids", false);

	std::set<protein::ResidueID> filled_set;
	for(std::map<protein::ResidueID, contacto::ResidueContactAreaDifferenceScore>::const_iterator it=profile.begin();it!=profile.end();++it)
	{
		if(it->second.has_valid_ratio())
		{
			filled_set.insert(it->first);
		}
	}

	int target_residues_count=0;
	int target_atoms_count=0;
	int target_used_residues_count=0;
	int target_used_atoms_count=0;

	for(std::map<protein::ResidueID, protein::ResidueSummary>::const_iterator it=residue_ids_of_target.begin();it!=residue_ids_of_target.end();++it)
	{
		target_residues_count++;
		target_atoms_count+=it->second.atoms_count;
		if(filled_set.count(it->first)>0)
		{
			target_used_residues_count++;
			target_used_atoms_count+=it->second.atoms_count;
		}
	}

	int model_residues_count=0;
	int model_atoms_count=0;
	int model_used_residues_count=0;
	int model_used_atoms_count=0;

	for(std::map<protein::ResidueID, protein::ResidueSummary>::const_iterator it=residue_ids_of_model.begin();it!=residue_ids_of_model.end();++it)
	{
		model_residues_count++;
		model_atoms_count+=it->second.atoms_count;
		if(filled_set.count(it->first)>0)
		{
			model_used_residues_count++;
			model_used_atoms_count+=it->second.atoms_count;
		}
	}

	std::cout << "t_res " << target_residues_count << "\n";
	std::cout << "t_res_used " << target_used_residues_count << "\n";
	std::cout << "t_atoms " << target_atoms_count << "\n";
	std::cout << "t_atoms_used " << target_used_atoms_count << "\n";
	std::cout << "m_res " << model_residues_count << "\n";
	std::cout << "m_res_used " << model_used_residues_count << "\n";
	std::cout << "m_atoms " << model_atoms_count << "\n";
	std::cout << "m_atoms_used " << model_used_atoms_count << "\n";
}
