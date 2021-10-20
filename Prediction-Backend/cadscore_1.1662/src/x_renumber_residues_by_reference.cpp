#include <iostream>

#include "protein/atom.h"
#include "protein/residue_ids_collection.h"
#include "protein/sequence_tools.h"

#include "auxiliaries/command_line_options.h"
#include "auxiliaries/std_containers_io.h"

namespace
{

struct Residue
{
	protein::ResidueID id;
	protein::ResidueSummary summary;
	Residue(const protein::ResidueID& id, const protein::ResidueSummary& summary) : id(id), summary(summary)
	{
	}
};

typedef std::vector<Residue> ResidueVector;

struct Overlay
{
	int score;
	std::map<protein::ResidueID, protein::ResidueID> model_to_target;
	Overlay() : score(0)
	{
	}

	static Overlay combine_overlays(const Overlay& a, const Overlay& b)
	{
		Overlay ab=a;
		ab.score+=b.score;
		ab.model_to_target.insert(b.model_to_target.begin(), b.model_to_target.end());
		return ab;
	}
};

ResidueVector collect_residues(const std::vector<protein::Atom>& atoms)
{
	const std::map<protein::ResidueID, protein::ResidueSummary> residues_map=protein::collect_residue_ids_from_atoms(atoms);
	ResidueVector residues;
	for(std::map<protein::ResidueID, protein::ResidueSummary>::const_iterator it=residues_map.begin();it!=residues_map.end();++it)
	{
		residues.push_back(Residue(it->first, it->second));
	}
	return residues;
}

std::vector<ResidueVector> subdivide_residues_by_chain(const ResidueVector& residues)
{
	std::map<std::string, ResidueVector> chains_map;
	for(ResidueVector::const_iterator it=residues.begin();it!=residues.end();++it)
	{
		chains_map[it->id.chain_id].push_back(*it);
	}
	std::vector<ResidueVector> divided_residues;
	for(std::map<std::string, ResidueVector>::const_iterator it=chains_map.begin();it!=chains_map.end();++it)
	{
		divided_residues.push_back(it->second);
	}
	return divided_residues;
}

std::string collect_sequence_from_residues(const ResidueVector& residues)
{
	std::ostringstream output;
	for(ResidueVector::const_iterator it=residues.begin();it!=residues.end();++it)
	{
		output << protein::sequence_tools::LettersCoding::convert_residue_codes_big_to_small(it->summary.name);
	}
	return output.str();
}

Overlay produce_overlay(const protein::sequence_tools::PairwiseSequenceAlignment::SimpleScorer& scorer, const ResidueVector& target_residues, const ResidueVector& model_residues, std::vector<ResidueVector>& leftovers)
{
	const std::string target_sequence=collect_sequence_from_residues(target_residues);
	const std::string model_sequence=collect_sequence_from_residues(model_residues);
	Overlay overlay;
	const std::vector< std::pair<int, int> > alignment=protein::sequence_tools::PairwiseSequenceAlignment::construct_sequence_alignment(target_sequence, model_sequence, scorer, true, &overlay.score);
	leftovers.clear();
	for(std::size_t i=0;i<alignment.size();i++)
	{
		const std::pair<int, int>& pairing=alignment[i];
		if(pairing.first>=0 && pairing.second>=0)
		{
			overlay.model_to_target[model_residues[pairing.second].id]=target_residues[pairing.first].id;
		}
		if(pairing.first<0 && pairing.second>=0)
		{
			if(leftovers.empty() || (i>0 && alignment[i-1].first>=0))
			{
				leftovers.push_back(ResidueVector());
			}
			leftovers.back().push_back(model_residues[pairing.second]);
		}
	}
	return overlay;
}

void recursive_calculate_best_combined_overlay(
		const protein::sequence_tools::PairwiseSequenceAlignment::SimpleScorer& scorer,
		const std::vector<ResidueVector>& target_divided_residues,
		const std::size_t i,
		const std::vector<ResidueVector>& model_divided_residues,
		const Overlay& previous_overlay,
		Overlay& currently_best_overlay)
{
	if(i>=target_divided_residues.size() || model_divided_residues.empty())
	{
		if(previous_overlay.score>currently_best_overlay.score)
		{
			currently_best_overlay=previous_overlay;
		}
	}
	else
	{
		bool ideal_produced_overlay=false;
		for(std::size_t j=0;j<model_divided_residues.size() && !ideal_produced_overlay;j++)
		{
			std::vector<ResidueVector> leftovers;
			const Overlay new_overlay=produce_overlay(scorer, target_divided_residues[i], model_divided_residues[j], leftovers);
			ideal_produced_overlay=(target_divided_residues[i].size()==model_divided_residues[j].size() && new_overlay.score==(scorer.match(1, 1)*static_cast<int>(target_divided_residues[i].size())));
			std::vector<ResidueVector> new_model_divided_residues;
			for(std::size_t e=0;e<model_divided_residues.size();e++)
			{
				if(e!=j)
				{
					new_model_divided_residues.push_back(model_divided_residues[e]);
				}
				else if(!leftovers.empty())
				{
					new_model_divided_residues.insert(new_model_divided_residues.end(), leftovers.begin(), leftovers.end());
				}
			}
			recursive_calculate_best_combined_overlay(scorer, target_divided_residues, i+1, new_model_divided_residues, Overlay::combine_overlays(previous_overlay, new_overlay), currently_best_overlay);
		}
		recursive_calculate_best_combined_overlay(scorer, target_divided_residues, i+1, model_divided_residues, previous_overlay, currently_best_overlay);
	}
}

}

void x_renumber_residues_by_reference(const auxiliaries::CommandLineOptions& clo)
{
	clo.check_allowed_options("--match: --mismatch: --gap-start: --gap-extension: --replace-residue-names --output-in-pdb-format --print-summary-log --print-detailed-log");

	const int match_score=clo.arg_or_default_value<int>("--match", 10);
	const int mismatch_score=clo.arg_or_default_value<int>("--mismatch", -10);
	const int gap_start_score=clo.arg_or_default_value<int>("--gap-start", -11);
	const int gap_extension_score=clo.arg_or_default_value<int>("--gap-extension", -1);
	const bool replace_residue_names=clo.isopt("--replace-residue-names");
	const bool output_in_pdb_format=clo.isopt("--output-in-pdb-format");
	const bool print_summary_log=clo.isopt("--print-summary-log");
	const bool print_detailed_log=clo.isopt("--print-detailed-log");

	const protein::sequence_tools::PairwiseSequenceAlignment::SimpleScorer scorer(match_score, mismatch_score, gap_start_score, gap_extension_score);

	const std::vector<protein::Atom> target_atoms=auxiliaries::STDContainersIO::read_vector<protein::Atom>(std::cin, "target atoms", "atoms", false);
	const std::vector<protein::Atom> model_atoms=auxiliaries::STDContainersIO::read_vector<protein::Atom>(std::cin, "model atoms", "atoms", false);

	const ResidueVector target_residues=collect_residues(target_atoms);
	const ResidueVector model_residues=collect_residues(model_atoms);

	const std::vector<ResidueVector> target_divided_residues=subdivide_residues_by_chain(target_residues);
	const std::vector<ResidueVector> model_divided_residues=subdivide_residues_by_chain(model_residues);

	Overlay best_combined_overlay;
	recursive_calculate_best_combined_overlay(scorer, target_divided_residues, 0, model_divided_residues, Overlay(), best_combined_overlay);

	std::vector<protein::Atom> renumbered_model_atoms;
	for(std::vector<protein::Atom>::const_iterator it=model_atoms.begin();it!=model_atoms.end();++it)
	{
		protein::Atom atom=(*it);
		const protein::ResidueID residue_id=protein::ResidueID::from_atom(atom);
		if(best_combined_overlay.model_to_target.count(residue_id)>0)
		{
			const protein::ResidueID new_residue_id=best_combined_overlay.model_to_target.find(residue_id)->second;
			atom.chain_id=new_residue_id.chain_id;
			atom.residue_number=new_residue_id.residue_number;
			renumbered_model_atoms.push_back(atom);
		}
	}

	if(replace_residue_names)
	{
		const std::map<protein::ResidueID, protein::ResidueSummary> target_residue_summaries=protein::collect_residue_ids_from_atoms(target_atoms);
		for(std::vector<protein::Atom>::iterator it=renumbered_model_atoms.begin();it!=renumbered_model_atoms.end();++it)
		{
			protein::Atom& atom=(*it);
			const protein::ResidueID residue_id=protein::ResidueID::from_atom(atom);
			if(target_residue_summaries.count(residue_id)>0)
			{
				atom.residue_name=target_residue_summaries.find(residue_id)->second.name;
			}
		}
	}

	if(!renumbered_model_atoms.empty())
	{
		if(output_in_pdb_format)
		{
			for(std::vector<protein::Atom>::const_iterator it=renumbered_model_atoms.begin();it!=renumbered_model_atoms.end();++it)
			{
				std::cout << it->string_for_PDB_file(0.0) << "\n";
			}
			std::cout << "END\n";
		}
		else
		{
			auxiliaries::STDContainersIO::print_vector(std::cout, "atoms", renumbered_model_atoms);
		}
	}

	if(print_summary_log)
	{
		std::clog << "target_atoms " << target_atoms.size() << "\n";
		std::clog << "model_atoms " << model_atoms.size() << "\n";
		std::clog << "target_residues " << target_residues.size() << "\n";
		std::clog << "model_residues " << model_residues.size() << "\n";
		std::clog << "target_chains " << target_divided_residues.size() << "\n";
		std::clog << "model_chains " << model_divided_residues.size() << "\n";
		std::clog << "best_alignments_score " << best_combined_overlay.score << "\n";
		std::clog << "mapped_model_residues " << best_combined_overlay.model_to_target.size() << "\n";

		std::clog << "target_sequence         ";
		for(std::size_t i=0;i<target_divided_residues.size();i++)
		{
			if(!target_divided_residues[i].empty())
			{
				std::clog << target_divided_residues[i][0].id.chain_id << ":" << collect_sequence_from_residues(target_divided_residues[i]) << " ";
			}
		}
		std::clog << "\n";
		std::clog << "model_sequence          ";
		for(std::size_t i=0;i<model_divided_residues.size();i++)
		{
			if(!model_divided_residues[i].empty())
			{
				std::clog << model_divided_residues[i][0].id.chain_id << ":" << collect_sequence_from_residues(model_divided_residues[i]) << " ";
			}
		}
		std::clog << "\n";
		std::clog << "model_assigned_chains   ";
		for(std::size_t i=0;i<model_divided_residues.size();i++)
		{
			if(!model_divided_residues[i].empty())
			{
				std::clog << "  ";
				for(std::size_t j=0;j<model_divided_residues[i].size();j++)
				{
					if(best_combined_overlay.model_to_target.count(model_divided_residues[i][j].id)>0)
					{
						std::clog << best_combined_overlay.model_to_target.find(model_divided_residues[i][j].id)->second.chain_id;
					}
					else
					{
						std::clog << "-";
					}
				}
				std::clog << " ";
			}
		}
		std::clog << "\n";
	}

	if(print_detailed_log)
	{
		for(std::map<protein::ResidueID, protein::ResidueID>::const_iterator it=best_combined_overlay.model_to_target.begin();it!=best_combined_overlay.model_to_target.end();++it)
		{
			std::clog << it->first.chain_id << " " << it->first.residue_number << "   " << it->second.chain_id << " " << it->second.residue_number << "\n";
		}
	}
}
