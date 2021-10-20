#include <iostream>
#include <iomanip>
#include <cmath>

#include "protein/residue_id.h"
#include "contacto/contact_id.h"
#include "contacto/inter_residue_contact_areas.h"

#include "auxiliaries/command_line_options.h"
#include "auxiliaries/std_containers_io.h"

void x_print_inter_residue_contacts_graph(const auxiliaries::CommandLineOptions& clo)
{
	typedef std::map< contacto::ContactID<protein::ResidueID>, contacto::InterResidueContactAreas > InterResidueContacts;

	clo.check_allowed_options("--category:");

	const std::string category=clo.arg_or_default_value<std::string>("--category", "AA");

	const InterResidueContacts inter_residue_contacts=auxiliaries::STDContainersIO::read_map< contacto::ContactID<protein::ResidueID>, contacto::InterResidueContactAreas >(std::cin, "inter-residue contacts", "residue_contacts", false);

	for(InterResidueContacts::const_iterator it=inter_residue_contacts.begin();it!=inter_residue_contacts.end();++it)
	{
		if(it->second.areas.count(category)==1)
		{
			std::cout << it->first.a.chain_id << std::setfill('0') << std::setw(5) << it->first.a.residue_number << " ";
			std::cout << it->first.b.chain_id << std::setfill('0') << std::setw(5) << it->first.b.residue_number << " ";
			std::cout << ceil(it->second.areas.find(category)->second*100) << "\n";
		}
	}
}
