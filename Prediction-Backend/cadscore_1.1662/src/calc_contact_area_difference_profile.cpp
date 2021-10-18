#include "protein/atom.h"
#include "protein/residue_id.h"
#include "protein/residue_summary.h"

#include "contacto/residue_contact_area_difference_profile.h"
#include "contacto/residue_contact_area_difference_basic_scoring_functors.h"

#include "auxiliaries/command_line_options.h"
#include "auxiliaries/std_containers_io.h"

void calc_contact_area_difference_profile(const auxiliaries::CommandLineOptions& clo)
{
	clo.check_allowed_options("--type:");

	const int scoring_mode=clo.isopt("--type") ? clo.arg_in_interval<int>("--type", 0, 3) : 0;

	const std::map< contacto::ContactID<protein::ResidueID>, contacto::InterResidueContactDualAreas > combined_inter_residue_contacts=
			auxiliaries::STDContainersIO::read_map< contacto::ContactID<protein::ResidueID>, contacto::InterResidueContactDualAreas >(std::cin, "combined inter-residue contacts", "combined_residue_contacts", false);

	const std::map<protein::ResidueID, protein::ResidueSummary> residue_ids_1=auxiliaries::STDContainersIO::read_map<protein::ResidueID, protein::ResidueSummary>(std::cin, "target residue identifiers", "residue_ids", false);

	std::map<protein::ResidueID, contacto::ResidueContactAreaDifferenceScore> residue_contact_area_difference_profile;
	if(scoring_mode==0)
	{
		residue_contact_area_difference_profile=
				contacto::construct_residue_contact_area_difference_profile<protein::ResidueID, protein::ResidueSummary, contacto::BoundedDifferenceProducer, contacto::SimpleReferenceProducer>(combined_inter_residue_contacts, residue_ids_1);
	}
	else if(scoring_mode==1)
	{
		residue_contact_area_difference_profile=
				contacto::construct_residue_contact_area_difference_profile<protein::ResidueID, protein::ResidueSummary, contacto::SimpleDifferenceProducer, contacto::SimpleReferenceProducer>(combined_inter_residue_contacts, residue_ids_1);
	}
	else if(scoring_mode==2)
	{
		residue_contact_area_difference_profile=
				contacto::construct_residue_contact_area_difference_profile<protein::ResidueID, protein::ResidueSummary, contacto::SimpleDifferenceProducer, contacto::SummingReferenceProducer>(combined_inter_residue_contacts, residue_ids_1);
	}
	else if(scoring_mode==3)
	{
		residue_contact_area_difference_profile=
				contacto::construct_residue_contact_area_difference_profile<protein::ResidueID, protein::ResidueSummary, contacto::RawDifferenceProducer, contacto::SimpleReferenceProducer>(combined_inter_residue_contacts, residue_ids_1);
	}
	else
	{
		throw std::runtime_error("Invalid profile type");
	}

	if(residue_contact_area_difference_profile.empty())
	{
		throw std::runtime_error("No profile constructed");
	}
	else
	{
		auxiliaries::STDContainersIO::print_map(std::cout, "cad_profile", residue_contact_area_difference_profile, true);
	}
}
