#include <iostream>

#include "protein/basic_parsing.h"

#include "auxiliaries/command_line_options.h"

namespace
{

struct DSSPResidueRecord
{
	int residue_sequence_number;
	std::string insertion_code;
	std::string chain_name;
	std::string short_residue_name;
	std::string secondary_structure_summary;

	DSSPResidueRecord(const std::string& DSSP_file_line) :
		residue_sequence_number(protein::basic_parsing::convert_string<int>(protein::basic_parsing::substring_of_columned_file_line(DSSP_file_line, 6, 10))),
		insertion_code(protein::basic_parsing::substring_of_columned_file_line(DSSP_file_line, 11, 11)),
		chain_name(protein::basic_parsing::substring_of_columned_file_line(DSSP_file_line, 12, 12)),
		short_residue_name(protein::basic_parsing::substring_of_columned_file_line(DSSP_file_line, 14, 14)),
		secondary_structure_summary(protein::basic_parsing::substring_of_columned_file_line(DSSP_file_line, 17, 17))
	{
		if(chain_name.empty())
		{
			chain_name="?";
		}
		if(short_residue_name.empty())
		{
			secondary_structure_summary=" ";
		}
		if(secondary_structure_summary.empty())
		{
			secondary_structure_summary=" ";
		}
	}
};

std::vector<DSSPResidueRecord> read_DSSP_atom_records_from_DSSP_file_stream(std::istream& dssp_file_stream)
{
	std::vector<DSSPResidueRecord> records;
	bool records_started=false;
	while(dssp_file_stream.good())
	{
		std::string line;
		std::getline(dssp_file_stream, line);
		if(!records_started)
		{
			if(line.find("  #  RESIDUE AA STRUCTURE")==0)
			{
				records_started=true;
			}
		}
		else
		{
			const std::string short_residue_name=protein::basic_parsing::substring_of_columned_file_line(line, 14, 14);
			if(!short_residue_name.empty() && short_residue_name!="!")
			{
				try
				{
					const DSSPResidueRecord record(line);
					records.push_back(record);
				}
				catch(const std::exception& e)
				{
					std::cerr << "Invalid DSSP record in line: " << line << "\n";
				}
			}
		}
	}
	return records;
}

inline bool is_in_helix(const DSSPResidueRecord& dr)
{
	return (dr.secondary_structure_summary=="H" || dr.secondary_structure_summary=="G" || dr.secondary_structure_summary=="I");
}

inline bool is_in_strand(const DSSPResidueRecord& dr)
{
	return (dr.secondary_structure_summary=="B" || dr.secondary_structure_summary=="E");
}

}

void x_summarize_dssp_file(const auxiliaries::CommandLineOptions& clo)
{
	clo.check_allowed_options("--print-sequence --print-sequence-colored");

	std::vector<DSSPResidueRecord> dssp_records=read_DSSP_atom_records_from_DSSP_file_stream(std::cin);

	if(clo.isopt("--print-sequence"))
	{
		for(std::size_t i=0;i<dssp_records.size();i++)
		{
			std::cout << dssp_records[i].short_residue_name;
		}
		std::cout << "\n";
		for(std::size_t i=0;i<dssp_records.size();i++)
		{
			std::cout << dssp_records[i].secondary_structure_summary;
		}
		std::cout << "\n";
	}
	else if(clo.isopt("--print-sequence-colored"))
	{
		for(std::size_t i=0;i<dssp_records.size();i++)
		{
			const DSSPResidueRecord& dr=dssp_records[i];
			std::cout << "\033[30;";
			if(is_in_helix(dr))
			{
				std::cout << "41m";
			}
			else if(is_in_strand(dr))
			{
				std::cout << "42m";
			}
			else
			{
				std::cout << "47m";
			}
			std::cout << dr.short_residue_name;
			std::cout << "\033[0m";
		}
		std::cout << "\n";
	}
	else
	{
		std::size_t count_of_residues_in_helix=0;
		std::size_t count_of_residues_in_strand=0;
		for(std::size_t i=0;i<dssp_records.size();i++)
		{
			const DSSPResidueRecord& dr=dssp_records[i];
			if(is_in_helix(dr))
			{
				count_of_residues_in_helix++;
			}
			else if(is_in_strand(dr))
			{
				count_of_residues_in_strand++;
			}
		}
		std::cout << "all_residues_count " << dssp_records.size() << "\n";
		std::cout << "helix_residues_count " << count_of_residues_in_helix << "\n";
		std::cout << "strand_residues_count " << count_of_residues_in_strand << "\n";
	}
}
