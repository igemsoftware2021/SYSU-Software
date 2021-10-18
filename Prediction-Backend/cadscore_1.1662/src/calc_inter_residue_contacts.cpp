#include "protein/atom.h"
#include "protein/residue_id.h"
#include "protein/residue_summary.h"
#include "protein/residue_ids_intervals.h"

#include "contacto/inter_residue_contacts_construction.h"
#include "contacto/inter_residue_contacts_combination.h"
#include "contacto/inter_residue_contacts_filtering.h"
#include "contacto/residue_contact_area_difference_profile.h"
#include "contacto/residue_contact_area_difference_basic_scoring_functors.h"

#include "auxiliaries/command_line_options.h"
#include "auxiliaries/std_containers_io.h"

namespace
{

template<typename ContactsMap>
std::vector<std::string> collect_chain_names_from_contacts_map(const ContactsMap& contacts)
{
	std::set<std::string> set_of_names;
	for(typename ContactsMap::const_iterator it=contacts.begin();it!=contacts.end();++it)
	{
		set_of_names.insert(it->first.a.chain_id);
		set_of_names.insert(it->first.b.chain_id);
	}
	std::vector<std::string> vector_of_names;
	vector_of_names.insert(vector_of_names.end(), set_of_names.begin(), set_of_names.end());
	return vector_of_names;
}

}

void calc_inter_residue_contacts(const auxiliaries::CommandLineOptions& clo)
{
	clo.check_allowed_options("--inter-interval: --inter-chain --core --interface-zone --preserve-reflexive");

	const std::vector<protein::Atom> atoms=auxiliaries::STDContainersIO::read_vector<protein::Atom>(std::cin, "atoms", "atoms", false);

	const std::vector<contacto::InterAtomContact> inter_atom_contacts=auxiliaries::STDContainersIO::read_vector<contacto::InterAtomContact>(std::cin, "inter-atom contacts", "contacts", false);

	const std::map< contacto::ContactID<protein::ResidueID>, contacto::InterResidueContactAreas > unfiltered_inter_residue_contacts=contacto::construct_inter_residue_contacts<protein::Atom, protein::ResidueID>(atoms, inter_atom_contacts);

	std::map< contacto::ContactID<protein::ResidueID>, contacto::InterResidueContactAreas > inter_residue_contacts=unfiltered_inter_residue_contacts;
	contacto::filter_custom_contacts<protein::ResidueID, contacto::InterResidueContactAreas, protein::ResidueIDsIntervalsReader>(inter_residue_contacts, clo.isopt("--core"), clo.isopt("--interface-zone"), clo.isopt("--inter-chain"), (clo.isopt("--inter-interval") ? clo.arg<std::string>("--inter-interval") : std::string("")));

	if(clo.isopt("--preserve-reflexive"))
	{
		std::map< contacto::ContactID<protein::ResidueID>, contacto::InterResidueContactAreas > reflexive_contacts;
		for(std::map< contacto::ContactID<protein::ResidueID>, contacto::InterResidueContactAreas >::const_iterator it=inter_residue_contacts.begin();it!=inter_residue_contacts.end();++it)
		{
			{
				std::map< contacto::ContactID<protein::ResidueID>, contacto::InterResidueContactAreas >::const_iterator jt=unfiltered_inter_residue_contacts.find(contacto::ContactID<protein::ResidueID>(it->first.a, it->first.a));
				if(jt!=unfiltered_inter_residue_contacts.end())
				{
					reflexive_contacts.insert(*jt);
				}
			}
			{
				std::map< contacto::ContactID<protein::ResidueID>, contacto::InterResidueContactAreas >::const_iterator jt=unfiltered_inter_residue_contacts.find(contacto::ContactID<protein::ResidueID>(it->first.b, it->first.b));
				if(jt!=unfiltered_inter_residue_contacts.end())
				{
					reflexive_contacts.insert(*jt);
				}
			}
		}
		inter_residue_contacts.insert(reflexive_contacts.begin(), reflexive_contacts.end());
	}

	if(inter_residue_contacts.empty())
	{
		throw std::runtime_error("No inter-residue contacts constructed");
	}
	else
	{
		auxiliaries::STDContainersIO::print_map(std::cout, "residue_contacts", inter_residue_contacts, true);
	}
}

void calc_combined_inter_residue_contacts(const auxiliaries::CommandLineOptions& clo)
{
	typedef std::map< contacto::ContactID<protein::ResidueID>, contacto::InterResidueContactAreas > InterResidueContacts;
	typedef std::map< contacto::ContactID<protein::ResidueID>, contacto::InterResidueContactDualAreas > CombinedInterResidueContacts;

	clo.check_allowed_options("--optimally-rename-chains --binarize");

	const bool binarize=clo.isopt("--binarize");

	InterResidueContacts inter_residue_contacts_1=auxiliaries::STDContainersIO::read_map< contacto::ContactID<protein::ResidueID>, contacto::InterResidueContactAreas >(std::cin, "target inter-residue contacts", "residue_contacts", false);
	InterResidueContacts inter_residue_contacts_2=auxiliaries::STDContainersIO::read_map< contacto::ContactID<protein::ResidueID>, contacto::InterResidueContactAreas >(std::cin, "model inter-residue contacts", "residue_contacts", false);

	CombinedInterResidueContacts resulting_combined_inter_residue_contacts;

	bool renaming_performed=false;
	std::string renaming_comment="";
	if(clo.isopt("--optimally-rename-chains"))
	{
		const std::map<protein::ResidueID, protein::ResidueSummary> residue_ids_1=auxiliaries::STDContainersIO::read_map<protein::ResidueID, protein::ResidueSummary>(std::cin, "target residue identifiers", "residue_ids", false);;
		const std::vector<std::string> chain_names_1=collect_chain_names_from_contacts_map(inter_residue_contacts_1);
		const std::vector<std::string> chain_names_2=collect_chain_names_from_contacts_map(inter_residue_contacts_2);

		bool renaming_allowed=(chain_names_1.size()>1 && chain_names_1.size()==chain_names_2.size());
		for(std::size_t j=0;j<chain_names_1.size() && renaming_allowed;j++)
		{
			renaming_allowed=(chain_names_1[j]==chain_names_2[j]);
		}

		if(renaming_allowed)
		{
			std::pair<double, CombinedInterResidueContacts> best_variation(-1.0, CombinedInterResidueContacts());
			std::vector<std::string> chain_names_permutation=chain_names_1;
			do
			{
				InterResidueContacts inter_residue_contacts_2_with_renamed_chains;
				for(InterResidueContacts::const_iterator it=inter_residue_contacts_2.begin();it!=inter_residue_contacts_2.end();++it)
				{
					contacto::ContactID<protein::ResidueID> cid=it->first;
					bool a_renamed=false;
					bool b_renamed=false;
					for(std::size_t j=0;j<chain_names_2.size() && !(a_renamed && b_renamed);j++)
					{
						if(!a_renamed && cid.a.chain_id==chain_names_2[j])
						{
							cid.a.chain_id=chain_names_permutation[j];
							a_renamed=true;
						}
						if(!b_renamed && cid.b.chain_id==chain_names_2[j])
						{
							cid.b.chain_id=chain_names_permutation[j];
							b_renamed=true;
						}
					}
					inter_residue_contacts_2_with_renamed_chains[cid]=it->second;
				}

				CombinedInterResidueContacts combined_inter_residue_contacts=contacto::combine_two_inter_residue_contact_maps<protein::ResidueID>(inter_residue_contacts_1, inter_residue_contacts_2_with_renamed_chains, binarize);

				const std::map<protein::ResidueID, contacto::ResidueContactAreaDifferenceScore> residue_contact_area_difference_profile=contacto::construct_residue_contact_area_difference_profile<protein::ResidueID, protein::ResidueSummary, contacto::BoundedDifferenceProducer, contacto::SimpleReferenceProducer>(combined_inter_residue_contacts, residue_ids_1);
				const contacto::ResidueContactAreaDifferenceScore global_score=contacto::calculate_global_contact_area_difference_score_from_profile(residue_contact_area_difference_profile, false);
				const contacto::Ratio ratio=global_score.ratio("AA");
				if(ratio.reference>0.0)
				{
					const double score_from_ratio=(1-(ratio.difference/ratio.reference));
					if(score_from_ratio>best_variation.first)
					{
						best_variation.first=score_from_ratio;
						best_variation.second=combined_inter_residue_contacts;
						{
							std::ostringstream renaming_comment_stream;
							renaming_comment_stream << "Renamed chains from ( ";
							for(std::size_t j=0;j<chain_names_2.size();j++)
							{
								renaming_comment_stream << chain_names_2[j] << " ";
							}
							renaming_comment_stream << ") to ( ";
							for(std::size_t j=0;j<chain_names_2.size();j++)
							{
								renaming_comment_stream << chain_names_permutation[j] << " ";
							}
							renaming_comment_stream << ")";
							renaming_comment=renaming_comment_stream.str();
						}
					}
				}

				std::clog << "Tried renaming chains: ( ";
				for(std::size_t j=0;j<chain_names_2.size();j++)
				{
					std::clog << chain_names_2[j] << " ";
				}
				std::clog << ") -> ( ";
				for(std::size_t j=0;j<chain_names_2.size();j++)
				{
					std::clog << chain_names_permutation[j] << " ";
				}
				std::clog << "), got AA-score == " << (ratio.reference>0.0 ? (1-(ratio.difference/ratio.reference)) : 0.0) << "\n";
			}
			while(std::next_permutation(chain_names_permutation.begin(), chain_names_permutation.end()));

			if(best_variation.first>=0.0)
			{
				resulting_combined_inter_residue_contacts=best_variation.second;
			}

			renaming_performed=true;
		}
		else
		{
			std::clog << "Chains renaming was not possible\n";
		}
	}

	if(!renaming_performed)
	{
		resulting_combined_inter_residue_contacts=contacto::combine_two_inter_residue_contact_maps<protein::ResidueID>(inter_residue_contacts_1, inter_residue_contacts_2, binarize);
	}

	if(resulting_combined_inter_residue_contacts.empty())
	{
		throw std::runtime_error("No combined inter-residue contacts constructed");
	}
	else
	{
		auxiliaries::STDContainersIO::print_file_comment(std::cout, "This file contains combined inter-residue contact areas");
		auxiliaries::STDContainersIO::print_file_comment(std::cout, "calculated by accumulating inter-atom contact areas.");
		auxiliaries::STDContainersIO::print_file_comment(std::cout, "");
		auxiliaries::STDContainersIO::print_file_comment(std::cout, "The file is structured as follows:");
		auxiliaries::STDContainersIO::print_file_comment(std::cout, "  file_header (always equals 'combined_residue_contacts')");
		auxiliaries::STDContainersIO::print_file_comment(std::cout, "  N (the number of inter-residue contacts)");
		auxiliaries::STDContainersIO::print_file_comment(std::cout, "  contact_record[1]");
		auxiliaries::STDContainersIO::print_file_comment(std::cout, "  contact_record[2]");
		auxiliaries::STDContainersIO::print_file_comment(std::cout, "  ...");
		auxiliaries::STDContainersIO::print_file_comment(std::cout, "  contact_record[N]");
		auxiliaries::STDContainersIO::print_file_comment(std::cout, "");
		auxiliaries::STDContainersIO::print_file_comment(std::cout, "Each contact record has the following format:");
		auxiliaries::STDContainersIO::print_file_comment(std::cout, "  first_residue_chain_name first_residue_number second_residue_chain_name second_residue_number");
		auxiliaries::STDContainersIO::print_file_comment(std::cout, "  M (the number of contact types with non-zero areas)");
		auxiliaries::STDContainersIO::print_file_comment(std::cout, "  contact_type[1] corresponding_area_in_target corresponding_area_in_model");
		auxiliaries::STDContainersIO::print_file_comment(std::cout, "  contact_type[2] corresponding_area_in_target corresponding_area_in_model");
		auxiliaries::STDContainersIO::print_file_comment(std::cout, "  ...");
		auxiliaries::STDContainersIO::print_file_comment(std::cout, "  contact_type[M] corresponding_area_in_target corresponding_area_in_model");
		auxiliaries::STDContainersIO::print_file_comment(std::cout, "");
		auxiliaries::STDContainersIO::print_file_comment(std::cout, "If chain name is not available, it is denoted as '?'.");
		auxiliaries::STDContainersIO::print_file_comment(std::cout, "If first_residue_chain_name==second_residue_chain_name and first_residue_number==second_residue_number,");
		auxiliaries::STDContainersIO::print_file_comment(std::cout, "then a contact record describes residue solvent-accessible surface areas.");
		auxiliaries::STDContainersIO::print_file_comment(std::cout, "");
		auxiliaries::STDContainersIO::print_file_comment(std::cout, "Contact types are two-letter strings indicating what residue parts are in contact.");
		auxiliaries::STDContainersIO::print_file_comment(std::cout, "Each residue part is coded as a single letter:");
		auxiliaries::STDContainersIO::print_file_comment(std::cout, "  A - all residue atoms");
		auxiliaries::STDContainersIO::print_file_comment(std::cout, "  M - main chain");
		auxiliaries::STDContainersIO::print_file_comment(std::cout, "  S - side chain");
		auxiliaries::STDContainersIO::print_file_comment(std::cout, "There are 9 possible inter-residue contact types: AA, AM, AS, MA, MM, MS, SA, SM, SS.");
		auxiliaries::STDContainersIO::print_file_comment(std::cout, "There are 3 possible solvent-accessible surface types: AW, MW, SW.");
		if(!renaming_comment.empty())
		{
			auxiliaries::STDContainersIO::print_file_comment(std::cout, "");
			auxiliaries::STDContainersIO::print_file_comment(std::cout, renaming_comment);
			auxiliaries::STDContainersIO::print_file_comment(std::cout, "");
		}
		std::cout << "\n";
		auxiliaries::STDContainersIO::print_map(std::cout, "combined_residue_contacts", resulting_combined_inter_residue_contacts, true);
	}
}
