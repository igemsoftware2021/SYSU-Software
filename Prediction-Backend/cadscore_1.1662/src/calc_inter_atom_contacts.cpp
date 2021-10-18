#include <iostream>

#include "protein/atom.h"

#include "apollota/triangulation.h"
#include "apollota/utilities_for_triangulation.h"
#include "apollota/inter_sphere_contact_surface_on_sphere.h"

#include "contacto/inter_atom_contact.h"

#include "auxiliaries/command_line_options.h"
#include "auxiliaries/std_containers_io.h"

void calc_inter_atom_contacts(const auxiliaries::CommandLineOptions& clo)
{
	clo.check_allowed_options("--depth: --probe:");

	const std::size_t subdivision_depth=clo.isopt("--depth") ? clo.arg_in_interval<std::size_t>("--depth", 1, 4) : 3;
	const double probe_radius=clo.isopt("--probe") ? clo.arg_with_min_value<double>("--probe", 0) : 1.4;

	const std::vector<protein::Atom> atoms=auxiliaries::STDContainersIO::read_vector<protein::Atom>(std::cin, "atoms", "atoms", false);

	if(atoms.size()<4)
	{
		throw std::runtime_error("Less than 4 atoms provided");
	}

	const std::vector< std::vector<std::size_t> > graph=apollota::UtilitiesForTriangulation::collect_neighbors_graph_from_neighbors_map(apollota::UtilitiesForTriangulation::collect_neighbors_map_from_quadruples_map(apollota::Triangulation::construct_result(apollota::UtilitiesForTriangulation::collect_simple_spheres(atoms), 3.5, true, false).quadruples_map), atoms.size());

	for(std::size_t i=0;i<graph.size();i++)
	{
		if(graph[i].empty())
		{
			std::clog << "Sphere was not included into the Voronoi diagram: " << atoms[i].string_for_human_reading() << "\n";
		}
	}

	const std::vector<contacto::InterAtomContact> inter_atom_contacts=apollota::InterSphereContactSurfaceOnSphere::construct_inter_sphere_contacts_from_surface_areas<contacto::InterAtomContact>(
			apollota::InterSphereContactSurfaceOnSphere::calculate_surface_areas(atoms, graph, subdivision_depth, probe_radius));

	if(atoms.empty() || inter_atom_contacts.empty())
	{
		throw std::runtime_error("No inter-atom contacts constructed");
	}
	else
	{
		auxiliaries::STDContainersIO::print_vector(std::cout, "atoms", atoms);
		auxiliaries::STDContainersIO::print_vector(std::cout, "contacts", inter_atom_contacts);
	}
}
