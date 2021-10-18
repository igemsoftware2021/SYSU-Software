#include <iostream>
#include <cstdlib>

#include "protein/atom.h"
#include "protein/residue_id.h"

#include "contacto/inter_atom_contact.h"
#include "contacto/contact_classification.h"

#include "auxiliaries/command_line_options.h"
#include "auxiliaries/std_containers_io.h"

void summarize_inter_atom_contacts(const auxiliaries::CommandLineOptions& clo)
{
	clo.check_allowed_options("");

	const std::vector<protein::Atom> atoms=auxiliaries::STDContainersIO::read_vector<protein::Atom>(std::cin, "atoms", "atoms", false);

	const std::vector<contacto::InterAtomContact> inter_atom_contacts=auxiliaries::STDContainersIO::read_vector<contacto::InterAtomContact>(std::cin, "inter atom contacts", "contacts", false);

	std::map<std::string, double> values;
	for(std::size_t i=0;i<inter_atom_contacts.size();i++)
	{
		const contacto::InterAtomContact& contact=inter_atom_contacts[i];
		const protein::Atom& a=atoms.at(contact.a);
		const protein::Atom& b=atoms.at(contact.b);

		if(protein::ResidueID::from_atom(a)!=protein::ResidueID::from_atom(b))
		{
			const std::vector<std::string> contact_classes=contacto::ContactClassification::classify_atoms_contact<protein::Atom, protein::ResidueID>(a, b);
			for(std::size_t j=0;j<contact_classes.size();j++)
			{
				const std::string& contact_class=contact_classes[j];
				values[std::string("all_")+contact_class]+=contact.area;
				if(a.chain_id!=b.chain_id)
				{
					values[std::string("inter_chain_")+a.chain_id+"_"+b.chain_id+"_"+contact_class]+=contact.area;
				}
			}
		}
		else if(a==b)
		{
			values["sas"]+=contact.area;
		}
	}

	auxiliaries::STDContainersIO::print_map(std::cout, "summary_values", values, false);
}
