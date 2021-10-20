#ifndef PROTEIN_ATOM_H_
#define PROTEIN_ATOM_H_

#include <string>
#include <iostream>
#include <sstream>

#include "basic_parsing.h"

namespace protein
{

struct Atom
{
	enum MoleculeClass {unidentified_molecule_class, amino_acid, nucleotide};
	enum LocationClass {unidentified_location_class, main_chain, side_chain};

	std::string chain_id;
	int atom_number;
	int residue_number;
	std::string residue_name;
	std::string atom_name;
	double x;
	double y;
	double z;
	double r;
	int molecule_class;
	int location_class;

	Atom() :
		chain_id("?"),
		atom_number(0),
		residue_number(0),
		residue_name("?"),
		atom_name("?"),
		x(0),
		y(0),
		z(0),
		r(0),
		molecule_class(unidentified_molecule_class),
		location_class(unidentified_location_class) {}

	bool operator==(const Atom& a) const
	{
		return (this==&a ||
				(atom_number==a.atom_number &&
				residue_number==a.residue_number &&
				chain_id==a.chain_id &&
				residue_name==a.residue_name &&
				atom_name==a.atom_name &&
				x==a.x &&
				y==a.y &&
				z==a.z &&
				r==a.r &&
				molecule_class==a.molecule_class &&
				location_class==a.location_class));
	}

	friend std::ostream& operator<<(std::ostream &output, const Atom &atom)
	{
		output << atom.chain_id << " ";
		output << atom.atom_number << " ";
		output << atom.residue_number << " ";
		output << atom.residue_name << " ";
		output << atom.atom_name << " ";
		output << atom.x << " ";
		output << atom.y << " ";
		output << atom.z << " ";
		output << atom.r << " ";
		output << atom.molecule_class << " ";
		output << atom.location_class;
		return output;
	}

	friend std::istream& operator>>(std::istream &input, Atom &atom)
	{
		input >> atom.chain_id;
		input >> atom.atom_number;
		input >> atom.residue_number;
		input >> atom.residue_name;
		input >> atom.atom_name;
		input >> atom.x;
		input >> atom.y;
		input >> atom.z;
		input >> atom.r;
		input >> atom.molecule_class;
		input >> atom.location_class;
		return input;
	}

	std::string string_for_human_reading() const
	{
		std::ostringstream output;
		output << "atom number = " << atom_number << ", ";
		output << "chain = " << chain_id << ", ";
		output << "residue number = " << residue_number << ", ";
		output << "residue name = " << residue_name << ", ";
		output << "atom name = " << atom_name << ", ";
		output << "x = " << x << ", ";
		output << "y = " << y << ", ";
		output << "z = " << z << ", ";
		output << "radius = " << r;
		return output.str();
	}

	std::string string_for_PDB_file(double temperature_factor) const
	{
		std::string line(80, ' ');
		basic_parsing::insert_string_to_columned_file_line("ATOM", 1, 6, false, line);
		basic_parsing::insert_string_to_columned_file_line(basic_parsing::convert_int_to_string(atom_number), 7, 11, true, line);
		basic_parsing::insert_string_to_columned_file_line(atom_name, (atom_name.size()>3 ? 13 : 14), 16, false, line);
		basic_parsing::insert_string_to_columned_file_line(residue_name, 18, 20, false, line);
		basic_parsing::insert_string_to_columned_file_line((chain_id=="?" ? std::string(" ") : chain_id), 22, 22, false, line);
		basic_parsing::insert_string_to_columned_file_line(basic_parsing::convert_int_to_string(residue_number), 23, 26, true, line);
		basic_parsing::insert_string_to_columned_file_line(basic_parsing::convert_double_to_string(x, 3), 31, 38, true, line);
		basic_parsing::insert_string_to_columned_file_line(basic_parsing::convert_double_to_string(y, 3), 39, 46, true, line);
		basic_parsing::insert_string_to_columned_file_line(basic_parsing::convert_double_to_string(z, 3), 47, 54, true, line);
		basic_parsing::insert_string_to_columned_file_line(basic_parsing::convert_double_to_string(temperature_factor, 2), 61, 66, true, line);
		return line;
	}
};

}

#endif /* PROTEIN_ATOM_H_ */
