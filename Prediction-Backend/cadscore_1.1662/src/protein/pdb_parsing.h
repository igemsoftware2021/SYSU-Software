#ifndef PROTEIN_PDBPARSING_H_
#define PROTEIN_PDBPARSING_H_

#include <vector>

#include "basic_parsing.h"

namespace protein
{

struct PDBAtomRecord
{
	std::string label;
	int atom_serial_number;
	std::string name;
	std::string alternate_location_indicator;
	std::string residue_name;
	std::string chain_name;
	int residue_sequence_number;
	std::string insertion_code;
	double x;
	double y;
	double z;
	double temperature_factor;
	std::string element;

	PDBAtomRecord(const std::string& PDB_file_line) :
		label(basic_parsing::substring_of_columned_file_line(PDB_file_line, 1, 6)),
		atom_serial_number(basic_parsing::convert_string<int>(basic_parsing::substring_of_columned_file_line(PDB_file_line, 7, 11))),
		name(basic_parsing::substring_of_columned_file_line(PDB_file_line, 13, 16)),
		alternate_location_indicator(basic_parsing::substring_of_columned_file_line(PDB_file_line, 17, 17)),
		residue_name(basic_parsing::substring_of_columned_file_line(PDB_file_line, 18, 20)),
		chain_name(basic_parsing::substring_of_columned_file_line(PDB_file_line, 22, 22)),
		residue_sequence_number(basic_parsing::convert_string<int>(basic_parsing::substring_of_columned_file_line(PDB_file_line, 23, 26))),
		insertion_code(basic_parsing::substring_of_columned_file_line(PDB_file_line, 27, 27)),
		x(basic_parsing::convert_string<double>(basic_parsing::substring_of_columned_file_line(PDB_file_line, 31, 38))),
		y(basic_parsing::convert_string<double>(basic_parsing::substring_of_columned_file_line(PDB_file_line, 39, 46))),
		z(basic_parsing::convert_string<double>(basic_parsing::substring_of_columned_file_line(PDB_file_line, 47, 54))),
		temperature_factor(basic_parsing::safe_convert_string<double>(basic_parsing::substring_of_columned_file_line(PDB_file_line, 61, 66), 0)),
		element(basic_parsing::substring_of_columned_file_line(PDB_file_line, 77, 78))
	{
		if (label.empty() || name.empty() || residue_name.empty())
		{
			throw std::runtime_error("Atom record has not enough string data");
		}
		if(chain_name.empty())
		{
			chain_name="?";
		}
	}

	std::string generate_PDB_file_line() const
	{
		std::string line(80, ' ');
		basic_parsing::insert_string_to_columned_file_line(label, 1, 6, false, line);
		basic_parsing::insert_string_to_columned_file_line(basic_parsing::convert_int_to_string(atom_serial_number), 7, 11, true, line);
		basic_parsing::insert_string_to_columned_file_line(name, (name.size()>3 ? 13 : 14), 16, false, line);
		basic_parsing::insert_string_to_columned_file_line(alternate_location_indicator, 17, 17, false, line);
		basic_parsing::insert_string_to_columned_file_line(residue_name, 18, 20, false, line);
		basic_parsing::insert_string_to_columned_file_line((chain_name=="?" ? std::string(" ") : chain_name), 22, 22, false, line);
		basic_parsing::insert_string_to_columned_file_line(basic_parsing::convert_int_to_string(residue_sequence_number), 23, 26, true, line);
		basic_parsing::insert_string_to_columned_file_line(insertion_code, 27, 27, false, line);
		basic_parsing::insert_string_to_columned_file_line(basic_parsing::convert_double_to_string(x, 3), 31, 38, true, line);
		basic_parsing::insert_string_to_columned_file_line(basic_parsing::convert_double_to_string(y, 3), 39, 46, true, line);
		basic_parsing::insert_string_to_columned_file_line(basic_parsing::convert_double_to_string(z, 3), 47, 54, true, line);
		basic_parsing::insert_string_to_columned_file_line(basic_parsing::convert_double_to_string(temperature_factor, 2), 61, 66, true, line);
		//basic_parsing::insert_string_to_columned_file_line(element, 77, 78, true, line);
		return line;
	}

	std::string generate_TER_PDB_file_line() const
	{
		std::string line(80, ' ');
		basic_parsing::insert_string_to_columned_file_line("TER", 1, 6, false, line);
		basic_parsing::insert_string_to_columned_file_line(basic_parsing::convert_int_to_string(atom_serial_number+1), 7, 11, true, line);
		basic_parsing::insert_string_to_columned_file_line(residue_name, 18, 20, false, line);
		basic_parsing::insert_string_to_columned_file_line((chain_name=="?" ? std::string(" ") : chain_name), 22, 22, false, line);
		basic_parsing::insert_string_to_columned_file_line(basic_parsing::convert_int_to_string(residue_sequence_number), 23, 26, true, line);
		return line;
	}
};


inline std::vector<PDBAtomRecord> read_PDB_atom_records_from_PDB_file_stream(std::istream& pdb_file_stream)
{
	std::vector<PDBAtomRecord> records;
	while(pdb_file_stream.good())
	{
		std::string line;
		std::getline(pdb_file_stream, line);
		const std::string label=basic_parsing::substring_of_columned_file_line(line, 1, 6);
		if(label=="ATOM" || label=="HETATM")
		{
			try
			{
				const PDBAtomRecord record(line);
				records.push_back(record);
			}
			catch(const std::exception& e)
			{
				std::cerr << "Invalid atom record in line: " << line << "\n";
			}
		}
		else if(label=="END" || label=="ENDMDL")
		{
			return records;
		}
	}
	return records;
}

}

#endif /* PROTEIN_PDBPARSING_H_ */
