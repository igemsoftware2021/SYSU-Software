#include <iostream>
#include <map>
#include <cmath>

#include "protein/atom.h"
#include "protein/atom_id.h"

#include "contacto/contact_id.h"
#include "contacto/inter_atom_contact.h"
#include "contacto/utilities.h"
#include "contacto/ratio.h"

#include "auxiliaries/command_line_options.h"
#include "auxiliaries/std_containers_io.h"

namespace
{

std::map<contacto::ContactID<protein::AtomID>, double> construct_inter_atom_contacts_map(const std::vector<protein::Atom>& atoms, const std::vector<contacto::InterAtomContact>& inter_atom_contacts, bool inter_chain)
{
	std::map<contacto::ContactID<protein::AtomID>, double> contacts_map;
	for(std::size_t i=0;i<inter_atom_contacts.size();i++)
	{
		const contacto::InterAtomContact& contact=inter_atom_contacts[i];
		const protein::AtomID a_id=protein::AtomID::from_atom(atoms[contact.a]);
		const protein::AtomID b_id=protein::AtomID::from_atom(atoms[contact.b]);
		if(a_id.residue_id!=b_id.residue_id)
		{
			if(!inter_chain || a_id.residue_id.chain_id!=b_id.residue_id.chain_id)
			{
				contacts_map[contacto::ContactID<protein::AtomID>(a_id, b_id)]=contact.area;
			}
		}
	}
	return contacts_map;
}

void print_atoms_local_scores_as_pdb(const std::vector<protein::Atom>& atoms, const std::map<protein::AtomID, contacto::Ratio>& local_ratios, bool normalized)
{
	for(std::size_t i=0;i<atoms.size();i++)
	{
		const protein::Atom& a=atoms[i];
		double value=-1.0;
		const std::map<protein::AtomID, contacto::Ratio>::const_iterator it=local_ratios.find(protein::AtomID::from_atom(a));
		if(it!=local_ratios.end())
		{
			const contacto::Ratio& local_ratio=it->second;
			if(local_ratio.reference>0.0)
			{
				if(normalized)
				{
					value=local_ratio.difference/local_ratio.reference*100.0;
				}
				else
				{
					value=local_ratio.difference;
				}
			}
			else
			{
				value=-2.0;
			}
		}
		std::cout << a.string_for_PDB_file(value) << "\n";
	}
	std::cout << "END\n";
}

}

void calc_inter_atom_contact_area_difference_score(const auxiliaries::CommandLineOptions& clo)
{
	clo.check_allowed_options("--local --local-as-pdb-of-target --local-as-pdb-of-model --local-normalized-for-pdb --global --inter-chain");

	const bool print_local=clo.isopt("--local");
	const bool print_local_as_pdb_of_target=clo.isopt("--local-as-pdb-of-target");
	const bool print_local_as_pdb_of_model=clo.isopt("--local-as-pdb-of-model");
	const bool local_normalized_for_pdb=clo.isopt("--local-normalized-for-pdb");
	const bool print_any_local=(print_local || print_local_as_pdb_of_target || print_local_as_pdb_of_model);
	const bool print_global=clo.isopt("--global") || !print_any_local;
	const bool inter_chain=clo.isopt("--inter-chain");

	const std::vector<protein::Atom> atoms_1=auxiliaries::STDContainersIO::read_vector<protein::Atom>(std::cin, "target atoms", "atoms", false);

	const std::vector<contacto::InterAtomContact> inter_atom_contacts_1=auxiliaries::STDContainersIO::read_vector<contacto::InterAtomContact>(std::cin, "target inter-atom contacts", "contacts", false);

	const std::vector<protein::Atom> atoms_2=auxiliaries::STDContainersIO::read_vector<protein::Atom>(std::cin, "model atoms", "atoms", false);

	const std::vector<contacto::InterAtomContact> inter_atom_contacts_2=auxiliaries::STDContainersIO::read_vector<contacto::InterAtomContact>(std::cin, "model inter-atom contacts", "contacts", false);

	typedef std::map< contacto::ContactID<protein::AtomID>, std::pair<double, double> > CombinedContactsMap;

	const CombinedContactsMap combined_contacts=contacto::combine_two_maps(
			construct_inter_atom_contacts_map(atoms_1, inter_atom_contacts_1, inter_chain),
			construct_inter_atom_contacts_map(atoms_2, inter_atom_contacts_2, inter_chain));

	if(print_global)
	{
		contacto::Ratio global_ratio;
		for(CombinedContactsMap::const_iterator it=combined_contacts.begin();it!=combined_contacts.end();++it)
		{
			const double t=it->second.first;
			const double m=it->second.second;
			global_ratio.difference+=std::min(fabs(t-m), t);
			global_ratio.reference+=t;
		}
		std::cout << "inter_atom_cad_global_score " << (global_ratio.reference>0.0 ? (1-(global_ratio.difference/global_ratio.reference)) : 0.0) << "\n";
	}

	if(print_any_local)
	{
		std::map<protein::AtomID, contacto::Ratio> local_ratios;
		for(CombinedContactsMap::const_iterator it=combined_contacts.begin();it!=combined_contacts.end();++it)
		{
			const double t=it->second.first;
			const double m=it->second.second;
			contacto::Ratio& local_ratio=local_ratios[it->first.a];
			local_ratio.difference+=std::min(fabs(t-m), t);
			local_ratio.reference+=t;
		}

		if(print_local)
		{
			auxiliaries::STDContainersIO::print_map(std::cout, "inter_atom_cad_local_scores", local_ratios, false);
		}

		if(print_local_as_pdb_of_target)
		{
			print_atoms_local_scores_as_pdb(atoms_1, local_ratios, local_normalized_for_pdb);
		}

		if(print_local_as_pdb_of_model)
		{
			print_atoms_local_scores_as_pdb(atoms_2, local_ratios, local_normalized_for_pdb);
		}
	}
}
