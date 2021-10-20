#ifndef RESIDUE_SUMMARY_H_
#define RESIDUE_SUMMARY_H_

#include <string>
#include <iostream>

namespace protein
{

struct ResidueSummary
{
	std::string name;
	int atoms_count;

	ResidueSummary() : atoms_count(0)
	{
	}

	ResidueSummary(const std::string& name, int atoms_count) : name(name), atoms_count(atoms_count)
	{
	}

	friend std::ostream& operator<<(std::ostream &output, const ResidueSummary &rs)
	{
		output << rs.name << " ";
		output << rs.atoms_count;
		return output;
	}

	friend std::istream& operator>>(std::istream &input, ResidueSummary &rs)
	{
		input >> rs.name;
		input >> rs.atoms_count;
		return input;
	}
};

}

#endif /* RESIDUE_SUMMARY_H_ */
