#include <iostream>

#include "protein/atom.h"

#include "apollota/triangulation.h"
#include "apollota/utilities_for_triangulation.h"
#include "apollota/inter_sphere_contact_face_on_hyperboloid.h"
#include "apollota/inter_sphere_contact_surface_on_sphere.h"

#include "contacto/inter_atom_contact.h"

#include "auxiliaries/command_line_options.h"
#include "auxiliaries/std_containers_io.h"

void calc_inter_atom_faces(const auxiliaries::CommandLineOptions& clo)
{
	typedef apollota::InterSphereContactFaceOnHyperboloid CellFace;

	clo.check_allowed_options("--probe: --step: --projections: --add-sas");

	const double probe_radius=clo.isopt("--probe") ? clo.arg_with_min_value<double>("--probe", 0) : 1.4;
	const double step_length=clo.isopt("--step") ? clo.arg_with_min_value<double>("--step", 0.1) : 0.7;
	const int projections_count=clo.isopt("--projections") ? clo.arg_with_min_value<int>("--projections", 5) : 5;
	const bool add_sas=clo.isopt("--add-sas");

	const std::vector<protein::Atom> atoms=auxiliaries::STDContainersIO::read_vector<protein::Atom>(std::cin, "atoms", "atoms", false);

	if(atoms.size()<4)
	{
		throw std::runtime_error("Less than 4 atoms provided");
	}

	const apollota::Triangulation::Result triangulation_result=apollota::Triangulation::construct_result(apollota::UtilitiesForTriangulation::collect_simple_spheres(atoms), 3.5, true, false);

	const apollota::UtilitiesForTriangulation::PairsNeighborsMap pairs_neighbours_map=apollota::UtilitiesForTriangulation::collect_pairs_neighbors_map_from_quadruples_map(triangulation_result.quadruples_map);

	std::set<contacto::InterAtomContact> inter_atom_contacts;

	for(apollota::UtilitiesForTriangulation::PairsNeighborsMap::const_iterator it=pairs_neighbours_map.begin();it!=pairs_neighbours_map.end();++it)
	{
		const std::size_t a_id=it->first.get(0);
		const std::size_t b_id=it->first.get(1);

		const protein::Atom& a=atoms[a_id];
		const protein::Atom& b=atoms[b_id];

		std::vector<const protein::Atom*> cs;
		cs.reserve(it->second.size());
		for(apollota::UtilitiesForTriangulation::PairsNeighborsMap::mapped_type::const_iterator jt=it->second.begin();jt!=it->second.end();++jt)
		{
			cs.push_back(&(atoms[*jt]));
		}

		const double face_area=CellFace::construct(a, b, cs, probe_radius, step_length, projections_count).area();
		if(apollota::greater(face_area, 0.0))
		{
			inter_atom_contacts.insert(contacto::InterAtomContact(a_id, b_id, face_area));
			inter_atom_contacts.insert(contacto::InterAtomContact(b_id, a_id, face_area));
		}
	}

	if(add_sas)
	{
		const std::vector< std::vector<std::size_t> > graph=apollota::UtilitiesForTriangulation::collect_neighbors_graph_from_neighbors_map(apollota::UtilitiesForTriangulation::collect_neighbors_map_from_quadruples_map(triangulation_result.quadruples_map), atoms.size());
		const std::vector<contacto::InterAtomContact> inter_atom_spherical_contacts=apollota::InterSphereContactSurfaceOnSphere::construct_inter_sphere_contacts_from_surface_areas<contacto::InterAtomContact>(apollota::InterSphereContactSurfaceOnSphere::calculate_surface_areas(atoms, graph, 3, probe_radius));
		for(std::size_t i=0;i<inter_atom_spherical_contacts.size();i++)
		{
			if(inter_atom_spherical_contacts[i].a==inter_atom_spherical_contacts[i].b)
			{
				inter_atom_contacts.insert(inter_atom_spherical_contacts[i]);
			}
		}
	}

	if(atoms.empty() || inter_atom_contacts.empty())
	{
		throw std::runtime_error("No inter-atom faces constructed");
	}
	else
	{
		auxiliaries::STDContainersIO::print_vector(std::cout, "atoms", atoms);
		auxiliaries::STDContainersIO::print_set(std::cout, "contacts", inter_atom_contacts);
	}
}
