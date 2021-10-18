#include <iostream>

#include "protein/atom.h"
#include "protein/residue_ids_collection.h"

#include "auxiliaries/command_line_options.h"
#include "auxiliaries/std_containers_io.h"

void filter_atoms_by_target(const auxiliaries::CommandLineOptions& clo)
{
	clo.check_allowed_options("--detailed --print-rejected --allow-unmatched-residue-names");

	const bool detailed=clo.isopt("--detailed");
	const bool print_rejected=clo.isopt("--print-rejected");
	const bool allow_unmatched_residue_names=clo.isopt("--allow-unmatched-residue-names");

	const std::vector<protein::Atom> atoms_of_model=auxiliaries::STDContainersIO::read_vector<protein::Atom>(std::cin, "model atoms", "atoms", false);

	const std::vector<protein::Atom> atoms_of_target=auxiliaries::STDContainersIO::read_vector<protein::Atom>(std::cin, "target atoms", "atoms", false);

	std::vector<protein::Atom> result;
	result.reserve(atoms_of_model.size());

	if(detailed)
	{
		const std::map<protein::ResidueID, std::vector<std::size_t> > residue_ids_atoms_of_target=group_atoms_indices_by_residue_ids(atoms_of_target);
		for(std::size_t i=0;i<atoms_of_model.size();i++)
		{
			protein::Atom atom_of_model=atoms_of_model[i];
			std::map<protein::ResidueID, std::vector<std::size_t> >::const_iterator it=residue_ids_atoms_of_target.find(protein::ResidueID::from_atom(atom_of_model));
			if(it!=residue_ids_atoms_of_target.end())
			{
				const std::vector<std::size_t>& target_residue_atoms_ids=it->second;
				bool passed=false;
				for(std::size_t j=0;j<target_residue_atoms_ids.size() && !passed;j++)
				{
					const protein::Atom& atom_of_target=atoms_of_target[target_residue_atoms_ids[j]];
					if(atom_of_model.residue_name==atom_of_target.residue_name)
					{
						if(atom_of_model.atom_name==atom_of_target.atom_name)
						{
							passed=true;
						}
						else if(atom_of_model.atom_name.size()==3 &&
								atom_of_model.atom_name.size()==atom_of_target.atom_name.size() &&
								atom_of_model.atom_name.substr(0,2)==atom_of_target.atom_name.substr(0,2) &&
								((atom_of_model.atom_name[2]=='\'' && atom_of_target.atom_name[2]=='*') || (atom_of_model.atom_name[2]=='*' && atom_of_target.atom_name[2]=='\'')))
						{
							atom_of_model.atom_name=atom_of_target.atom_name;
							passed=true;
						}
						else if(
								(atom_of_model.atom_name=="OP1" && atom_of_target.atom_name=="O1P") ||
								(atom_of_model.atom_name=="OP2" && atom_of_target.atom_name=="O2P") ||
								(atom_of_model.atom_name=="OP3" && atom_of_target.atom_name=="O3P") ||
								(atom_of_model.atom_name=="O1P" && atom_of_target.atom_name=="OP1") ||
								(atom_of_model.atom_name=="O2P" && atom_of_target.atom_name=="OP2") ||
								(atom_of_model.atom_name=="O3P" && atom_of_target.atom_name=="OP3"))
						{
							atom_of_model.atom_name=atom_of_target.atom_name;
							passed=true;
						}
					}
				}
				if(passed)
				{
					result.push_back(atom_of_model);
				}
				else if(print_rejected)
				{
					std::clog << "Rejected: " << atom_of_model.string_for_human_reading() << "\n";
				}
			}
		}
	}
	else
	{
		const std::map<protein::ResidueID, protein::ResidueSummary> residue_ids_of_target=protein::collect_residue_ids_from_atoms(atoms_of_target);
		for(std::size_t i=0;i<atoms_of_model.size();i++)
		{
			const protein::Atom& atom=atoms_of_model[i];
			std::map<protein::ResidueID, protein::ResidueSummary>::const_iterator it=residue_ids_of_target.find(protein::ResidueID::from_atom(atom));
			if(it!=residue_ids_of_target.end())
			{
				if(atom.residue_name==it->second.name || allow_unmatched_residue_names)
				{
					result.push_back(atom);
				}
				else
				{
					std::ostringstream output;
					output << "Model atom chain name and residue number matched the target, but model atom residue name did not: " << atom.string_for_human_reading();
					throw std::runtime_error(output.str());
				}
			}
		}
	}

	if(result.empty())
	{
		throw std::runtime_error("Chain naming and/or residue numbering of the model did not match the target, therefore all the model atoms were rejected");
	}
	else
	{
		auxiliaries::STDContainersIO::print_vector(std::cout, "atoms", result);
	}
}

void filter_atoms_by_name(const auxiliaries::CommandLineOptions& clo)
{
	clo.check_allowed_options("--name:");

	const std::string name=clo.arg<std::string>("--name");

	const std::vector<protein::Atom> atoms=auxiliaries::STDContainersIO::read_vector<protein::Atom>(std::cin, "atoms", "atoms", false);

	std::vector<protein::Atom> result;
	result.reserve(atoms.size());

	for(std::size_t i=0;i<atoms.size();i++)
	{
		const protein::Atom& atom=atoms[i];
		if(atom.atom_name==name)
		{
			result.push_back(atom);
		}
	}

	if(result.empty())
	{
		throw std::runtime_error("No atoms found for the given name");
	}
	else
	{
		auxiliaries::STDContainersIO::print_vector(std::cout, "atoms", result);
	}
}
