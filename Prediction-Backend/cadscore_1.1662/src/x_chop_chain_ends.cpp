#include <iostream>
#include <fstream>
#include <map>

#include "protein/pdb_parsing.h"

#include "auxiliaries/command_line_options.h"

void x_chop_chain_ends(const auxiliaries::CommandLineOptions& clo)
{
	typedef std::map<int, std::vector<std::size_t> > ChainResiduesMap;
	typedef std::map<std::string, ChainResiduesMap> ChainsMap;

	clo.check_allowed_options("--output-prefix:");

	const std::string output_prefix=clo.arg<std::string>("--output-prefix");

	const std::vector<protein::PDBAtomRecord> records=protein::read_PDB_atom_records_from_PDB_file_stream(std::cin);

	ChainsMap chains_map;

	for(std::size_t i=0;i<records.size();i++)
	{
		if(records[i].label=="ATOM")
		{
			chains_map[records[i].chain_name][records[i].residue_sequence_number].push_back(i);
		}
	}

	if(chains_map.empty())
	{
		std::cerr << "Exiting because there are no chains.\n";
		return;
	}

	for(ChainsMap::const_iterator it=chains_map.begin();it!=chains_map.end();++it)
	{
		if(it->second.size()<40)
		{
			std::cerr << "Exiting because there is a too short chain.\n";
			return;
		}
	}

	for(std::size_t chunk_percents=0;chunk_percents<=24;chunk_percents+=3)
	{
		std::ostringstream filename_stream;
		filename_stream << output_prefix << "_" << chunk_percents;
		std::string filename=filename_stream.str();
		std::ofstream output(filename.c_str(), std::ios::out);
		for(ChainsMap::const_iterator it=chains_map.begin();it!=chains_map.end();++it)
		{
			const ChainResiduesMap& chain_residues_map=it->second;
			const std::size_t chunk_size=chunk_percents*chain_residues_map.size()/100;
			if(chunk_size<chain_residues_map.size())
			{
				std::size_t num=0;
				for(ChainResiduesMap::const_iterator jt=chain_residues_map.begin();jt!=chain_residues_map.end();++jt)
				{
					if((num>=chunk_size) && ((num+chunk_size)<chain_residues_map.size()))
					{
						const std::vector<std::size_t>& records_ids=jt->second;
						for(std::size_t i=0;i<records_ids.size();i++)
						{
							output << records[records_ids[i]].generate_PDB_file_line() << "\n";
						}
					}
					num++;
				}
			}
		}
		output << "END\n";
	}
}
