#include <iostream>

#include "protein/atom.h"
#include "protein/residue_ids_collection.h"
#include "protein/sequence_tools.h"

#include "contacto/inter_atom_contact.h"

#include "auxiliaries/command_line_options.h"
#include "auxiliaries/std_containers_io.h"

namespace
{

std::string collect_string_without_gaps(const std::string& input)
{
	std::string output;
	for(std::size_t i=0;i<input.size();i++)
	{
		if(input[i]!='-')
		{
			output.push_back(input[i]);
		}
	}
	return output;
}

template<typename T>
std::vector< std::pair<typename T::key_type, typename T::mapped_type> > filter_residue_ids_map_by_chain_id(const T& residue_ids_map, const std::string& chain_id)
{
	std::vector< std::pair<typename T::key_type, typename T::mapped_type> > filtered_residue_ids_map;
	for(typename T::const_iterator it=residue_ids_map.begin();it!=residue_ids_map.end();++it)
	{
		if(it->first.chain_id==chain_id)
		{
			filtered_residue_ids_map.push_back(*it);
		}
	}
	return filtered_residue_ids_map;
}

template<typename T>
std::string collect_sequence_string_from_residue_ids_map(const T& residue_ids_map)
{
	std::ostringstream output;
	for(typename T::const_iterator it=residue_ids_map.begin();it!=residue_ids_map.end();++it)
	{
		output << protein::sequence_tools::LettersCoding::convert_residue_codes_big_to_small(it->second.name);
	}
	return output.str();
}

}

void x_renumber_residues_in_inter_atom_contacts(const auxiliaries::CommandLineOptions& clo)
{
	clo.check_allowed_options("--target-chain-names: --model-chain-names: --match: --mismatch: --gap-start: --gap-extension:");

	const std::string chain_names_1=clo.arg<std::string>("--target-chain-names");
	const std::string chain_names_2=clo.arg<std::string>("--model-chain-names");
	const int match_score=clo.arg_or_default_value<int>("--match", 10);
	const int mismatch_score=clo.arg_or_default_value<int>("--mismatch", -10);
	const int gap_start_score=clo.arg_or_default_value<int>("--gap-start", -11);
	const int gap_extension_score=clo.arg_or_default_value<int>("--gap-extension", -1);

	std::vector<protein::Atom> atoms_1=auxiliaries::STDContainersIO::read_vector<protein::Atom>(std::cin, "target atoms", "atoms", false);
	const std::vector<contacto::InterAtomContact> inter_atom_contacts_1=auxiliaries::STDContainersIO::read_vector<contacto::InterAtomContact>(std::cin, "target inter-atom contacts", "contacts", false);
	std::vector<protein::Atom> atoms_2=auxiliaries::STDContainersIO::read_vector<protein::Atom>(std::cin, "model atoms", "atoms", false);
	const std::vector<contacto::InterAtomContact> inter_atom_contacts_2=auxiliaries::STDContainersIO::read_vector<contacto::InterAtomContact>(std::cin, "model inter-atom contacts", "contacts", false);

	const std::map<protein::ResidueID, protein::ResidueSummary> residue_ids_1=protein::collect_residue_ids_from_atoms(atoms_1);
	const std::map<protein::ResidueID, std::vector<std::size_t> > residue_ids_indices_1=protein::group_atoms_indices_by_residue_ids(atoms_1);
	const std::map<protein::ResidueID, protein::ResidueSummary> residue_ids_2=protein::collect_residue_ids_from_atoms(atoms_2);
	const std::map<protein::ResidueID, std::vector<std::size_t> > residue_ids_indices_2=protein::group_atoms_indices_by_residue_ids(atoms_2);

	for(std::size_t l=0;l<chain_names_1.size() && l<chain_names_2.size();l++)
	{
		const std::string chain_name_1=chain_names_1.substr(l, 1);
		const std::string chain_name_2=chain_names_2.substr(l, 1);

		std::clog << "Chain in target: " << chain_name_1 << "\n";
		std::clog << "Chain in model: " << chain_name_2 << "\n";

		const std::vector< std::pair<protein::ResidueID, protein::ResidueSummary> > filtered_residue_ids_1=filter_residue_ids_map_by_chain_id(residue_ids_1, chain_name_1);
		const std::vector< std::pair<protein::ResidueID, protein::ResidueSummary> > filtered_residue_ids_2=filter_residue_ids_map_by_chain_id(residue_ids_2, chain_name_2);

		if(filtered_residue_ids_1.empty() || filtered_residue_ids_2.empty())
		{
			throw std::runtime_error("No atoms for chain name");
		}

		const std::string sequence_string_from_atoms_1=collect_sequence_string_from_residue_ids_map(filtered_residue_ids_1);
		const std::string sequence_string_from_atoms_2=collect_sequence_string_from_residue_ids_map(filtered_residue_ids_2);

		if(sequence_string_from_atoms_1.size()!=filtered_residue_ids_1.size() || sequence_string_from_atoms_2.size()!=filtered_residue_ids_2.size())
		{
			throw std::runtime_error("Unconvertible sequence");
		}

		std::clog << "Sequences from atoms:\n";
		std::clog << sequence_string_from_atoms_1 << "\n";
		std::clog << sequence_string_from_atoms_2 << "\n";

		std::string input_alignment_string_1;
		std::cin >> input_alignment_string_1;
		std::string input_alignment_string_2;
		std::cin >> input_alignment_string_2;

		const std::string sequence_string_from_input_alignment_1=collect_string_without_gaps(input_alignment_string_1);
		const std::string sequence_string_from_input_alignment_2=collect_string_without_gaps(input_alignment_string_2);

		if(input_alignment_string_1.empty() || input_alignment_string_1.size()!=input_alignment_string_2.size() || sequence_string_from_input_alignment_1.empty() || sequence_string_from_input_alignment_2.empty())
		{
			throw std::runtime_error("Invalid input alignment");
		}

		std::clog << "Input alignment:\n";
		std::clog << input_alignment_string_1 << "\n";
		std::clog << input_alignment_string_2 << "\n";

		const protein::sequence_tools::PairwiseSequenceAlignment::SimpleScorer scorer(match_score, mismatch_score, gap_start_score, gap_extension_score);
		std::vector< std::pair<int, int> > alignments[3];
		alignments[0]=protein::sequence_tools::PairwiseSequenceAlignment::construct_sequence_alignment(sequence_string_from_atoms_1, sequence_string_from_input_alignment_1, scorer);
		alignments[1]=protein::sequence_tools::PairwiseSequenceAlignment::read_sequence_alignment(input_alignment_string_1, input_alignment_string_2);
		alignments[2]=protein::sequence_tools::PairwiseSequenceAlignment::construct_sequence_alignment(sequence_string_from_input_alignment_2, sequence_string_from_atoms_2, scorer);

		for(int j=0;j<2;j++)
		{
			for(std::size_t i=0;i<std::max(alignments[j].size(), alignments[j+1].size());i++)
			{
				if(i<alignments[j].size() && alignments[j][i].second<0 && (i>=alignments[j+1].size() || alignments[j+1][i].first>=0))
				{
					alignments[j+1].insert(alignments[j+1].begin()+i, std::pair<int, int>(-1, -1));
				}
				else if(i<alignments[j+1].size() && alignments[j+1][i].first<0 && (i>=alignments[j].size() || alignments[j][i].second>=0))
				{
					for(int e=0;e<=j;e++)
					{
						alignments[e].insert(alignments[e].begin()+i, std::pair<int, int>(-1, -1));
					}
				}
			}
		}

		if(alignments[0].size()!=alignments[1].size() || alignments[1].size()!=alignments[2].size())
		{
			throw std::runtime_error("Failed to align alignments");
		}

		std::clog << "Result alignment:\n";
		protein::sequence_tools::PairwiseSequenceAlignment::print_sequence_alignment(sequence_string_from_atoms_1, sequence_string_from_input_alignment_1, alignments[0], std::clog);
		protein::sequence_tools::PairwiseSequenceAlignment::print_sequence_alignment(sequence_string_from_input_alignment_2, sequence_string_from_atoms_2, alignments[2], std::clog);

		const std::vector< std::pair<protein::ResidueID, std::vector<std::size_t> > > filtered_residue_ids_indices_1=filter_residue_ids_map_by_chain_id(residue_ids_indices_1, chain_name_1);
		const std::vector< std::pair<protein::ResidueID, std::vector<std::size_t> > > filtered_residue_ids_indices_2=filter_residue_ids_map_by_chain_id(residue_ids_indices_2, chain_name_2);

		for(std::size_t i=0;i<alignments[0].size();i++)
		{
			if(alignments[0][i].first>=0)
			{
				const std::vector<std::size_t>& atoms_indices=filtered_residue_ids_indices_1.at(alignments[0][i].first).second;
				for(std::size_t j=0;j<atoms_indices.size();j++)
				{
					protein::Atom& atom=atoms_1[atoms_indices[j]];
					atom.chain_id=chain_name_1;
					atom.residue_number=i+1;
				}
			}
			if(alignments[2][i].second>=0)
			{
				const std::vector<std::size_t>& atoms_indices=filtered_residue_ids_indices_2.at(alignments[2][i].second).second;
				for(std::size_t j=0;j<atoms_indices.size();j++)
				{
					protein::Atom& atom=atoms_2[atoms_indices[j]];
					atom.chain_id=chain_name_2;
					atom.residue_number=i+1;
				}
			}
		}

		std::clog << "\n";
	}

	auxiliaries::STDContainersIO::print_vector(std::cout, "atoms", atoms_1);
	auxiliaries::STDContainersIO::print_vector(std::cout, "contacts", inter_atom_contacts_1);
	auxiliaries::STDContainersIO::print_vector(std::cout, "atoms", atoms_2);
	auxiliaries::STDContainersIO::print_vector(std::cout, "contacts", inter_atom_contacts_2);
}
