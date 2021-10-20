#include <iostream>

#include "protein/atom.h"

#include "contacto/inter_atom_contact.h"

#include "auxiliaries/command_line_options.h"
#include "auxiliaries/std_containers_io.h"

void x_restrict_inter_chain_contacts(const auxiliaries::CommandLineOptions& clo)
{
	clo.check_allowed_options("--no-solvent");

	const bool no_solvent=clo.isopt("--no-solvent");

	const std::vector<protein::Atom> atoms=auxiliaries::STDContainersIO::read_vector<protein::Atom>(std::cin, "atoms", "atoms", false);

	const std::vector<contacto::InterAtomContact> inter_atom_contacts=auxiliaries::STDContainersIO::read_vector<contacto::InterAtomContact>(std::cin, "inter atom contacts", "contacts", false);

	std::vector<contacto::InterAtomContact> restricted_inter_atom_contacts;
	restricted_inter_atom_contacts.reserve(inter_atom_contacts.size()/8);

	for(std::size_t i=0;i<inter_atom_contacts.size();i++)
	{
		const contacto::InterAtomContact& contact=inter_atom_contacts[i];
		const protein::Atom& a=atoms.at(contact.a);
		const protein::Atom& b=atoms.at(contact.b);

		if(a.chain_id!=b.chain_id || (!no_solvent && contact.a==contact.b))
		{
			restricted_inter_atom_contacts.push_back(contact);
		}
	}

	if(atoms.empty() || restricted_inter_atom_contacts.empty())
	{
		throw std::runtime_error("No inter-chain contacts found");
	}
	else
	{
		auxiliaries::STDContainersIO::print_vector(std::cout, "atoms", atoms);
		auxiliaries::STDContainersIO::print_vector(std::cout, "contacts", restricted_inter_atom_contacts);
	}
}
