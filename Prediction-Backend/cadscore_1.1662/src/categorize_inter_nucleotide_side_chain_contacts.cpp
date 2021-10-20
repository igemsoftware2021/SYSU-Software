#include "protein/atom.h"
#include "protein/residue_ids_collection.h"

#include "apollota/basic_operations_on_spheres.h"

#include "contacto/contact_id.h"
#include "contacto/inter_residue_contact_areas.h"

#include "auxiliaries/command_line_options.h"
#include "auxiliaries/std_containers_io.h"

namespace
{

struct NucleotidePlane
{
	apollota::SimplePoint point;
	apollota::SimplePoint normal;

	static std::map<protein::ResidueID, NucleotidePlane> calc_nucleotides_planes(const std::vector<protein::Atom>& atoms)
	{
		std::map<protein::ResidueID, NucleotidePlane> planes;
		const std::map<protein::ResidueID, std::vector<std::size_t> > residue_ids_atoms=protein::group_atoms_indices_by_residue_ids(atoms);
		for(std::map<protein::ResidueID, std::vector<std::size_t> >::const_iterator it=residue_ids_atoms.begin();it!=residue_ids_atoms.end();++it)
		{
			const protein::ResidueID& residue_id=it->first;
			const std::vector<std::size_t> atoms_ids=it->second;
			if(!atoms_ids.empty())
			{
				if(atoms.at(atoms_ids.front()).molecule_class==protein::Atom::nucleotide)
				{
					std::vector<std::size_t> abc;
					for(std::size_t i=0;i<atoms_ids.size() && abc.size()<3;i++)
					{
						const std::size_t a_id=atoms_ids[i];
						const protein::Atom& a=atoms.at(a_id);
						if(a.atom_name=="C2" || a.atom_name=="C4" || a.atom_name=="C6")
						{
							abc.push_back(a_id);
						}
					}
					if(abc.size()==3)
					{
						const protein::Atom& a=atoms[abc[0]];
						const protein::Atom& b=atoms[abc[1]];
						const protein::Atom& c=atoms[abc[2]];
						const apollota::SimplePoint pa(a.x, a.y, a.z);
						const apollota::SimplePoint pb(b.x, b.y, b.z);
						const apollota::SimplePoint pc(c.x, c.y, c.z);
						NucleotidePlane& plane=planes[residue_id];
						plane.point=((pa+pb)+pc)*(1.0/3.0);
						plane.normal=((pb-pa)&(pc-pa)).unit();
					}
				}
			}
		}
		return planes;
	}

	friend std::ostream& operator<<(std::ostream &output, const NucleotidePlane& plane)
	{
		output << plane.point.x << " " << plane.point.y << " " << plane.point.z << " ";
		output << plane.normal.x << " " << plane.normal.y << " " << plane.normal.z;
		return output;
	}

	friend std::istream& operator>>(std::istream& input, NucleotidePlane& plane)
	{
		input >> plane.point.x >> plane.point.y >> plane.point.z;
		input >> plane.normal.x >> plane.normal.y >> plane.normal.z;
		return input;
	}
};

}

void categorize_inter_nucleotide_side_chain_contacts(const auxiliaries::CommandLineOptions& clo)
{
	clo.check_allowed_options("--diagnostic-output --use-atom-centers --output-normals-cos --output-rings-shift");

	const bool diagnostic_output=clo.isopt("--diagnostic-output");
	const bool use_atom_centers=clo.isopt("--use-atom-centers");
	const bool output_normals_cos=clo.isopt("--output-normals-cos");
	const bool output_rings_shift=clo.isopt("--output-rings-shift");

	const std::vector<protein::Atom> atoms=auxiliaries::STDContainersIO::read_vector<protein::Atom>(std::cin, "atoms", "atoms", false);

	std::map< contacto::ContactID<protein::ResidueID>, contacto::InterResidueContactAreas > inter_residue_contacts=
			auxiliaries::STDContainersIO::read_map< contacto::ContactID<protein::ResidueID>, contacto::InterResidueContactAreas >(std::cin, "inter-residue contacts", "residue_contacts", false);

	const std::map< protein::ResidueID, NucleotidePlane > nucleotides_planes=NucleotidePlane::calc_nucleotides_planes(atoms);

	const std::map< protein::ResidueID, std::vector<std::size_t> > residue_ids_atoms=protein::group_atoms_indices_by_residue_ids(atoms);

	for(std::map< contacto::ContactID<protein::ResidueID>, contacto::InterResidueContactAreas >::iterator it=inter_residue_contacts.begin();it!=inter_residue_contacts.end();++it)
	{
		const double area=it->second.area("SS");
		if(area>0.0)
		{
			const protein::ResidueID& rid_a=it->first.a;
			const protein::ResidueID& rid_b=it->first.b;
			if(nucleotides_planes.count(rid_a)==1 && residue_ids_atoms.count(rid_b)==1)
			{
				const NucleotidePlane& plane_a=nucleotides_planes.find(rid_a)->second;
				const std::vector<std::size_t> atoms_ids_b=residue_ids_atoms.find(rid_b)->second;
				std::vector<std::size_t> sc_atoms_ids_b;
				for(std::size_t i=0;i<atoms_ids_b.size();i++)
				{
					std::size_t atom_id=atoms_ids_b[i];
					if(atoms[atom_id].molecule_class==protein::Atom::nucleotide && atoms[atom_id].location_class==protein::Atom::side_chain)
					{
						sc_atoms_ids_b.push_back(atom_id);
					}
				}
				if(!sc_atoms_ids_b.empty())
				{
					const protein::Atom& first_atom=atoms[sc_atoms_ids_b[0]];
					const int first_halfspace=apollota::halfspace_of_sphere(plane_a.point, plane_a.normal, apollota::SimpleSphere(first_atom, (use_atom_centers ? 0.0 : first_atom.r)));
					bool one_halfspace=(first_halfspace!=0);
					for(std::size_t i=1;i<sc_atoms_ids_b.size() && one_halfspace;i++)
					{
						const protein::Atom& another_atom=atoms[sc_atoms_ids_b[i]];
						const int atom_halfspace=apollota::halfspace_of_sphere(plane_a.point, plane_a.normal, apollota::SimpleSphere(another_atom, (use_atom_centers ? 0.0 : another_atom.r)));
						if(atom_halfspace!=first_halfspace)
						{
							one_halfspace=false;
						}
					}
					if(one_halfspace)
					{
						if(first_halfspace<0)
						{
							it->second.areas["na_stacking_down"]=area;
						}
						else
						{
							it->second.areas["na_stacking_up"]=area;
						}
						it->second.areas["na_stacking"]=area;
					}
					else
					{
						it->second.areas["na_siding"]=area;
					}
				}
			}
		}
	}

	if(diagnostic_output)
	{
		std::vector<std::string> contact_types_of_interest;
		contact_types_of_interest.push_back("na_stacking");
		contact_types_of_interest.push_back("na_siding");
		for(std::map< contacto::ContactID<protein::ResidueID>, contacto::InterResidueContactAreas >::iterator it=inter_residue_contacts.begin();it!=inter_residue_contacts.end();++it)
		{
			for(std::size_t i=0;i<contact_types_of_interest.size();i++)
			{
				const std::string& contact_type=contact_types_of_interest[i];
				const double area=it->second.area(contact_type);
				if(area>0.0)
				{
					const protein::ResidueID& rid_a=it->first.a;
					const protein::ResidueID& rid_b=it->first.b;
					if(residue_ids_atoms.count(rid_a)>0 && residue_ids_atoms.count(rid_b)>0)
					{
						const std::vector<std::size_t>& atoms_ids_a=residue_ids_atoms.find(rid_a)->second;
						const std::vector<std::size_t>& atoms_ids_b=residue_ids_atoms.find(rid_b)->second;
						if(!atoms_ids_a.empty() && !atoms_ids_b.empty())
						{
							const std::size_t atom_id_a=atoms_ids_a.front();
							const std::size_t atom_id_b=atoms_ids_b.front();
							if(atom_id_a<atoms.size() && atom_id_b<atoms.size())
							{
								std::cout << rid_a.chain_id << " " << rid_a.residue_number << " " << atoms[atom_id_a].residue_name << " ";
								std::cout << rid_b.chain_id << " " << rid_b.residue_number << " " << atoms[atom_id_b].residue_name << " ";
								std::cout << contact_type << " " << area << " " << (rid_a<rid_b ? "left_to_right" : "right_to_left");
								if(output_normals_cos)
								{
									std::cout << " " << apollota::dot_product(nucleotides_planes.find(rid_a)->second.normal, nucleotides_planes.find(rid_b)->second.normal);
								}
								if(output_rings_shift)
								{
									const NucleotidePlane& plane_a=nucleotides_planes.find(rid_a)->second;
									const NucleotidePlane& plane_b=nucleotides_planes.find(rid_b)->second;
									const apollota::SimplePoint ab=(plane_b.point-plane_a.point);
									const double ab_proj=apollota::dot_product(ab, plane_a.normal);
									const double shift_dist=sqrt(apollota::squared_point_module(ab)-(ab_proj*ab_proj));
									const double shift_cos=ab_proj/ab.module();
									std::cout << " " << shift_dist << " " << shift_cos;
								}
								std::cout << "\n";
							}
						}
					}
				}
			}
		}
	}
	else
	{
		auxiliaries::STDContainersIO::print_map(std::cout, "residue_contacts", inter_residue_contacts, true);
	}
}
