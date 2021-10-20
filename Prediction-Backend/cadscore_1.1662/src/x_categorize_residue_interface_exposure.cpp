#include "protein/residue_id.h"

#include "contacto/contact_id.h"
#include "contacto/inter_residue_contact_areas.h"

#include "auxiliaries/command_line_options.h"
#include "auxiliaries/std_containers_io.h"

void x_categorize_residue_interface_exposure(const auxiliaries::CommandLineOptions& clo)
{
	clo.check_allowed_options("");

	const std::map< contacto::ContactID<protein::ResidueID>, contacto::InterResidueContactAreas > inter_residue_contacts=
			auxiliaries::STDContainersIO::read_map< contacto::ContactID<protein::ResidueID>, contacto::InterResidueContactAreas >(std::cin, "inter-residue contacts", "residue_contacts", false);

	std::map< contacto::ContactID<protein::ResidueID>, contacto::InterResidueContactAreas > extended_inter_residue_contacts=inter_residue_contacts;

	for(std::map< contacto::ContactID<protein::ResidueID>, contacto::InterResidueContactAreas >::const_iterator it=inter_residue_contacts.begin();it!=inter_residue_contacts.end();++it)
	{
		const protein::ResidueID& a_id=it->first.a;
		const protein::ResidueID& b_id=it->first.b;
		if(a_id.chain_id!=b_id.chain_id)
		{
			for(contacto::InterResidueContactAreas::AreasMap::const_iterator jt=it->second.areas.begin();jt!=it->second.areas.end();++jt)
			{
				const std::string& name=jt->first;
				const double area=jt->second;
				const contacto::ContactID<protein::ResidueID> self_contact_id(a_id, a_id);
				extended_inter_residue_contacts[self_contact_id].areas[name+"_exposed"]+=area;
				extended_inter_residue_contacts[self_contact_id].areas[name+"_exposed_to_"+b_id.chain_id]+=area;
			}
		}
	}

	auxiliaries::STDContainersIO::print_map(std::cout, "residue_contacts", extended_inter_residue_contacts, true);
}
