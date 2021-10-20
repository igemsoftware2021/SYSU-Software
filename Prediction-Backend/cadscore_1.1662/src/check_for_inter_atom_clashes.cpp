#include <iostream>

#include "protein/atom.h"
#include "protein/residue_id.h"

#include "apollota/search_for_spherical_collisions.h"

#include "auxiliaries/command_line_options.h"
#include "auxiliaries/std_containers_io.h"

namespace
{

inline bool check_if_atomic_collision_is_tolerable(const protein::Atom& a, const protein::Atom& b)
{
	if(protein::ResidueID::from_atom(a)==protein::ResidueID::from_atom(b))
	{
		return true;
	}
	else if(a.molecule_class==b.molecule_class)
	{
		if(a.residue_name=="CYS" && b.residue_name=="CYS" && a.atom_name=="SG" && b.atom_name=="SG")
		{
			return true;
		}
		if(a.chain_id==b.chain_id)
		{
			const std::string summary=a.atom_name+"-"+b.atom_name;
			if(a.molecule_class==static_cast<int>(protein::Atom::amino_acid))
			{
				return (summary=="C-N" || summary=="N-C");
			}
			else if(a.molecule_class==static_cast<int>(protein::Atom::nucleotide))
			{
				return (summary=="P-O3'" || summary=="P-O3*" || summary=="O3'-P" || summary=="O3*-P");
			}
		}
	}
	return false;
}

}

void check_for_inter_atom_clashes(const auxiliaries::CommandLineOptions& clo)
{
	typedef apollota::BoundingSpheresHierarchy Hierarchy;

	clo.check_allowed_options("--cutoff: --inter-chain --only-main-chain");

	const double cutoff=clo.isopt("--cutoff") ? clo.arg_with_min_value<double>("--cutoff", 0) : 0;
	const bool inter_chain=clo.isopt("--inter-chain");
	const bool only_main_chain=clo.isopt("--only-main-chain");

	std::vector<protein::Atom> atoms=auxiliaries::STDContainersIO::read_vector<protein::Atom>(std::cin, "atoms", "atoms", false);

	if(only_main_chain)
	{
		std::vector<protein::Atom> main_chain_atoms;
		for(std::size_t i=0;i<atoms.size();i++)
		{
			const protein::Atom& atom=atoms[i];
			if(atom.location_class==protein::Atom::main_chain)
			{
				main_chain_atoms.push_back(atom);
			}
		}
		if(main_chain_atoms.empty())
		{
			throw std::runtime_error("No main chain atoms provided");
		}
		atoms=main_chain_atoms;
	}

	const Hierarchy hierarchy(atoms, 3.5, 1);

	std::vector< std::pair< double, std::pair<std::size_t, std::size_t> > > collisions;
	std::set<std::size_t> colliding_atoms;
	for(std::size_t i=0;i+1<atoms.size();i++)
	{
		const protein::Atom& atom=atoms[i];
		const std::vector<std::size_t> neighbour_ids=apollota::SearchForSphericalCollisions::find_all_collisions(hierarchy, apollota::SimpleSphere(atom));
		for(std::size_t j=0;j<neighbour_ids.size();j++)
		{
			const std::size_t& neighbour_id=neighbour_ids[j];
			if(i<neighbour_id)
			{
				const protein::Atom& neighbour=atoms[neighbour_id];
				if(!check_if_atomic_collision_is_tolerable(atom, neighbour)
						&& (!inter_chain || atom.chain_id!=neighbour.chain_id))
				{
					const double overlap=0.0-apollota::minimal_distance_from_sphere_to_sphere(atom, neighbour);
					if(overlap>cutoff)
					{
						collisions.push_back(std::make_pair(overlap, std::make_pair(i, neighbour_id)));
						colliding_atoms.insert(i);
						colliding_atoms.insert(neighbour_id);
					}
				}
			}
		}
	}

	if(!collisions.empty())
	{
		std::cout << "collisions " << collisions.size() << "\n";
		std::cout << "colliding_atoms " << colliding_atoms.size() << "\n\n";

		std::sort(collisions.rbegin(), collisions.rend());

		for(std::size_t i=0;i<collisions.size();i++)
		{
			std::cout << "overlap " << collisions[i].first << "\n";
			const std::pair<std::size_t, std::size_t>& collision=collisions[i].second;
			std::cout << atoms[collision.first].string_for_human_reading() << "\n";
			std::cout << atoms[collision.second].string_for_human_reading() << "\n\n";
		}
	}
}
