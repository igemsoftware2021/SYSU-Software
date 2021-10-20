#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <memory>

#include "protein/atoms_reading.h"
#include "protein/atoms_classification.h"
#include "protein/residue_ids_collection.h"

#include "auxiliaries/command_line_options.h"
#include "auxiliaries/std_containers_io.h"

#include "resources/vdwr.h"

protein::VanDerWaalsRadiusAssigner construct_radius_assigner(const std::string& radius_classes_file_name, const std::string& radius_members_file_name)
{

#ifdef FOR_OLDER_COMPILERS
	typedef std::auto_ptr<std::istream> AutoPtrToInputStream;
#else
	typedef std::unique_ptr<std::istream> AutoPtrToInputStream;
#endif

	AutoPtrToInputStream radius_classes_stream;
	if(radius_classes_file_name.empty())
	{
		radius_classes_stream.reset(new std::istringstream(std::string(reinterpret_cast<const char*>(resources::vdwr_classes), resources::vdwr_classes_len)));
	}
	else
	{
		radius_classes_stream.reset(new std::ifstream(radius_classes_file_name.c_str()));
	}

	AutoPtrToInputStream radius_members_stream;
	if(radius_members_file_name.empty())
	{
		radius_members_stream.reset(new std::istringstream(std::string(reinterpret_cast<const char*>(resources::vdwr_members), resources::vdwr_members_len)));
	}
	else
	{
		radius_members_stream.reset(new std::ifstream(radius_members_file_name.c_str()));
	}

	const protein::VanDerWaalsRadiusAssigner radius_assigner(*radius_classes_stream, *radius_members_stream);

	return radius_assigner;
}

void collect_atoms(const auxiliaries::CommandLineOptions& clo)
{
	clo.check_allowed_options("--radius-classes: --radius-members: --HETATM --HOH --rename-chain: --auto-rename-chains --include-insertions");

	std::string radius_classes_file_name="";
	std::string radius_members_file_name="";
	if(clo.isopt("--radius-classes") || clo.isopt("--radius-members"))
	{
		radius_classes_file_name=clo.arg<std::string>("--radius-classes");
		radius_members_file_name=clo.arg<std::string>("--radius-members");
	}

	const bool include_heteroatoms=clo.isopt("--HETATM");
	const bool include_water=clo.isopt("--HOH");
	const std::string simple_chain_renaming=clo.isopt("--rename-chain") ? clo.arg<std::string>("--rename-chain") : std::string("");
	const bool auto_rename_chains=clo.isopt("--auto-rename-chains");
	const bool include_insertions=clo.isopt("--include-insertions");

	const protein::VanDerWaalsRadiusAssigner radius_assigner=construct_radius_assigner(radius_classes_file_name, radius_members_file_name);

	std::vector<protein::Atom> atoms=protein::AtomsReading::read_atoms_from_PDB_file_stream(std::cin, radius_assigner, include_heteroatoms, include_water, include_insertions);

	protein::AtomsClassification::classify_atoms(atoms);

	if(!simple_chain_renaming.empty())
	{
		for(std::size_t i=0;i<atoms.size();i++)
		{
			atoms[i].chain_id=simple_chain_renaming;
		}
	}

	if(auto_rename_chains)
	{
		std::map<std::string, std::string> chain_names_map;
		for(std::size_t i=0;i<atoms.size();i++)
		{
			chain_names_map[atoms[i].chain_id]="";
		}
		{
			const std::string letters="ABCDEFGHIJKLMNOPQRSTUVWXYZ";
			std::size_t letter_id=0;
			for(std::map<std::string, std::string>::iterator it=chain_names_map.begin();it!=chain_names_map.end();++it)
			{
				if(chain_names_map.size()<letters.size())
				{
					it->second=letters.substr(letter_id, 1);
				}
				else
				{
					std::ostringstream number_output;
					number_output << "n" << letter_id;
					it->second=number_output.str();
				}
				letter_id++;
			}
		}
		for(std::size_t i=0;i<atoms.size();i++)
		{
			std::string& atom_chain_id=atoms[i].chain_id;
			atom_chain_id=chain_names_map[atom_chain_id];
			if(atom_chain_id.empty())
			{
				atom_chain_id="x";
			}
		}
	}

	if(atoms.empty())
	{
		throw std::runtime_error("No atoms were collected from the provided PDB file stream");
	}
	else
	{
		auxiliaries::STDContainersIO::print_vector(std::cout, "atoms", atoms);
	}
}

void collect_residue_ids(const auxiliaries::CommandLineOptions& clo)
{
	clo.check_allowed_options("");

	const std::vector<protein::Atom> atoms=auxiliaries::STDContainersIO::read_vector<protein::Atom>(std::cin, "atoms", "atoms", false);

	const std::map<protein::ResidueID, protein::ResidueSummary> residue_ids=protein::collect_residue_ids_from_atoms(atoms);

	if(residue_ids.empty())
	{
		throw std::runtime_error("No residue identifiers were collected from the provided atoms stream");
	}
	else
	{
		auxiliaries::STDContainersIO::print_map(std::cout, "residue_ids", residue_ids, false);
	}
}

void merge_atoms(const auxiliaries::CommandLineOptions& clo)
{
	clo.check_allowed_options("");
	std::vector<protein::Atom> atoms=auxiliaries::STDContainersIO::read_vector<protein::Atom>(std::cin, "atoms", "atoms", true);
	while(auxiliaries::STDContainersIO::check_file_header(std::cin, "atoms"))
	{
		std::vector<protein::Atom> more_atoms=auxiliaries::STDContainersIO::read_vector<protein::Atom>(std::cin, "more atoms", "", true);
		atoms.insert(atoms.end(), more_atoms.begin(), more_atoms.end());
	}
	auxiliaries::STDContainersIO::print_vector(std::cout, "atoms", atoms);
}
