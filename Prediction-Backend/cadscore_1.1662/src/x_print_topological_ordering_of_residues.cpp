#include <iostream>
#include <map>
#include <iomanip>

#include "protein/pdb_parsing.h"

#include "apollota/basic_operations_on_points.h"

#include "auxiliaries/command_line_options.h"

namespace
{

struct FullResidueID
{
	std::string chain_name;
	int residue_sequence_number;
	std::string insertion_code;
	std::string residue_name;

	FullResidueID(const protein::PDBAtomRecord& record) :
		chain_name(record.chain_name),
		residue_sequence_number(record.residue_sequence_number),
		insertion_code(record.insertion_code),
		residue_name(record.residue_name)
	{
	}

	bool operator==(const FullResidueID& rid) const
	{
		return (residue_sequence_number==rid.residue_sequence_number && chain_name==rid.chain_name && insertion_code==rid.insertion_code);
	}

	bool operator< (const FullResidueID& rid) const
	{
		return (chain_name<rid.chain_name
				|| (chain_name==rid.chain_name && residue_sequence_number<rid.residue_sequence_number)
				|| (chain_name==rid.chain_name && residue_sequence_number==rid.residue_sequence_number && insertion_code<rid.insertion_code));
	}

	friend std::ostream& operator<<(std::ostream &output, const FullResidueID &rid)
	{
		output << std::setw(2) << (rid.chain_name.empty() ? std::string("?") : rid.chain_name) << " ";
		output << std::setw(6) << rid.residue_sequence_number << " ";
		output << std::setw(2) << (rid.insertion_code.empty() ? std::string(".") : rid.insertion_code) << " ";
		output << std::setw(4) << (rid.residue_name.empty() ? std::string("_") : rid.residue_name);
		return output;
	}
};

struct DoubleID
{
	int left;
	int right;

	DoubleID() : left(-1), right(-1)
	{
	}
};

}

void x_print_topological_ordering_of_residues(const auxiliaries::CommandLineOptions& clo)
{
	clo.check_allowed_options("--bond-distance: --print-ordered-pdb-file");

	const double bond_distance=clo.arg_or_default_value<double>("--bond-distance", 1.8);
	const bool print_ordered_pdb_file=clo.isopt("--print-ordered-pdb-file");

	const std::vector<protein::PDBAtomRecord> records=protein::read_PDB_atom_records_from_PDB_file_stream(std::cin);

	for(int t=0;t<2;t++)
	{
		std::map<FullResidueID, DoubleID> nodes_map;
		for(std::size_t i=0;i<records.size();i++)
		{
			const protein::PDBAtomRecord& record=records[i];
			if((record.alternate_location_indicator.empty() || record.alternate_location_indicator=="A") && record.residue_name!="HOH")
			{
				if(t==0)
				{
					if(record.name=="N")
					{
						nodes_map[FullResidueID(record)].left=i;
					}
					else if(record.name=="C")
					{
						nodes_map[FullResidueID(record)].right=i;
					}
				}
				else if(t==1)
				{
					if(record.name=="P")
					{
						nodes_map[FullResidueID(record)].left=i;
					}
					else if(record.name=="O3'" || record.name=="O3*")
					{
						nodes_map[FullResidueID(record)].right=i;
					}
				}
			}
		}

		std::vector< std::pair<FullResidueID, DoubleID> > nodes_vector;
		for(std::map<FullResidueID, DoubleID>::const_iterator it=nodes_map.begin();it!=nodes_map.end();++it)
		{
			if(it->second.left>=0 || it->second.right>=0)
			{
				nodes_vector.push_back(*it);
			}
		}

		std::vector<DoubleID> bonds(nodes_vector.size());
		for(std::size_t i=0;i<nodes_vector.size();i++)
		{
			bool found=false;
			for(std::size_t l=0;l<nodes_vector.size() && !found;l++)
			{
				const std::size_t j=(i+l<nodes_vector.size() ? (i+l) : (i+l-nodes_vector.size()));
				if(j!=i)
				{
					const int i_right=nodes_vector[i].second.right;
					const int j_left=nodes_vector[j].second.left;
					if(i_right>=0 && j_left>=0 && apollota::distance_from_point_to_point(records[i_right], records[j_left])<bond_distance)
					{
						found=true;
						bonds[i].right=j;
						bonds[j].left=i;
					}
				}
			}
		}

		std::vector< std::vector<FullResidueID> > chains;
		std::vector<bool> visited(nodes_vector.size(), false);
		for(std::size_t i=0;i<bonds.size();i++)
		{
			if(bonds[i].left<0 && !visited[i])
			{
				chains.push_back(std::vector<FullResidueID>());
				chains.back().push_back(nodes_vector[i].first);
				visited[i]=true;
				int next=bonds[i].right;
				while(next>=0 && !visited[next])
				{
					chains.back().push_back(nodes_vector[next].first);
					visited[next]=true;
					next=bonds[next].right;
				}
			}
		}

		if(print_ordered_pdb_file)
		{
			if(!chains.empty())
			{
				const std::string possible_chain_names="ABCDEFGHIJKLMNOPQRSTUVWXYZ";
				if(chains.size()>possible_chain_names.size())
				{
					std::cerr << "Too many chains to name.\n";
				}
				else
				{
					std::vector<protein::PDBAtomRecord> named_records=records;
					std::map< FullResidueID, std::vector<protein::PDBAtomRecord> > named_records_map;
					for(std::size_t i=0;i<records.size();i++)
					{
						named_records_map[FullResidueID(records[i])].push_back(records[i]);
					}
					for(std::size_t i=0;i<chains.size();i++)
					{
						const std::string chain_name=possible_chain_names.substr(i, 1);
						for(std::size_t j=0;j<chains[i].size();j++)
						{
							std::map< FullResidueID, std::vector<protein::PDBAtomRecord> >::iterator it=named_records_map.find(chains[i][j]);
							if(it!=named_records_map.end())
							{
								std::vector<protein::PDBAtomRecord>& residue_records=it->second;
								for(std::size_t u=0;u<residue_records.size();u++)
								{
									residue_records[u].chain_name=chain_name;
								}
							}
						}
					}
					std::clog << (t==0 ? "protein chains" : "nucleic acid chains") << " count " << chains.size() << "\n";
					for(std::map< FullResidueID, std::vector<protein::PDBAtomRecord> >::const_iterator it=named_records_map.begin();it!=named_records_map.end();++it)
					{
						const std::vector<protein::PDBAtomRecord>& residue_records=it->second;
						for(std::size_t u=0;u<residue_records.size();u++)
						{
							std::cout << residue_records[u].generate_PDB_file_line() << "\n";
						}
					}
				}
			}
		}
		else
		{
			std::cout << (t==0 ? "protein chains" : "nucleic acid chains") << " count " << chains.size() << "\n\n";
			for(std::size_t i=0;i<chains.size();i++)
			{
				std::cout << "chain size " << chains[i].size() << "\n";
				for(std::size_t j=0;j<chains[i].size();j++)
				{
					std::cout << chains[i][j] << "\n";
				}
				std::cout << "\n";
			}
		}
	}
}
