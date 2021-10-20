#include <iostream>
#include <sstream>
#include <map>

#include "protein/pdb_parsing.h"

#include "auxiliaries/command_line_options.h"

void x_rename_chains(const auxiliaries::CommandLineOptions& clo)
{
	clo.check_allowed_options("--before: --after:");

	const std::vector<std::string> before_list=clo.arg_vector<std::string>("--before", ',');
	const std::vector<std::string> after_list=clo.arg_vector<std::string>("--after", ',');

	if(before_list.empty() || before_list.size()!=after_list.size())
	{
		std::cerr << "Invalid parameters.\n";
		return;
	}

	std::vector<protein::PDBAtomRecord> records=protein::read_PDB_atom_records_from_PDB_file_stream(std::cin);

	std::map<std::string, std::string> before_to_after;
	for(std::size_t i=0;i<before_list.size();i++)
	{
		before_to_after[before_list[i]]=after_list[i];
		std::clog << before_list[i] << " -> " << after_list[i] << "\n";
	}

	for(std::size_t i=0;i<records.size();i++)
	{
		if(before_to_after.count(records[i].chain_name)>0)
		{
			records[i].chain_name=before_to_after.find(records[i].chain_name)->second;
		}
		std::cout << records[i].generate_PDB_file_line() << "\n";
	}
}
