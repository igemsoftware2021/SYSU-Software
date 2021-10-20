#include <iostream>

#include "protein/atom.h"

#include "apollota/triangulation.h"
#include "apollota/triangulation_output.h"
#include "apollota/utilities_for_triangulation.h"

#include "auxiliaries/command_line_options.h"
#include "auxiliaries/std_containers_io.h"

namespace
{

std::vector<apollota::SimpleSphere> read_raw_spheres_from_stream(std::istream& input)
{
	std::vector<apollota::SimpleSphere> result;
	while(input.good())
	{
		std::string line;
		std::getline(input, line);
		if(!line.empty())
		{
			std::istringstream line_input(line);
			apollota::SimpleSphere sphere;
			line_input >> sphere.x >> sphere.y >> sphere.z >> sphere.r;
			if(!line_input.fail())
			{
				result.push_back(sphere);
			}
		}
	}
	return result;
}

template<typename SphereType>
void calc_quadruples(
		const std::vector<SphereType>& input_atoms,
		const double bsi_init_radius,
		const bool use_one_radius,
		const bool augment,
		const bool skip_output,
		const bool print_log,
		const bool check)
{
	if(input_atoms.size()<4)
	{
		throw std::runtime_error("Less than 4 atoms provided");
	}

	const std::vector<SphereType>* atoms_ptr=(&input_atoms);

	std::vector<SphereType> modified_atoms;
	if(use_one_radius)
	{
		modified_atoms=input_atoms;
		for(std::size_t i=0;i<modified_atoms.size();i++)
		{
			modified_atoms[i].r=1;
		}
		atoms_ptr=(&modified_atoms);
	}

	const std::vector<SphereType>& atoms=(*atoms_ptr);

	const apollota::Triangulation::Result apollonius_triangulation_result=apollota::Triangulation::construct_result(apollota::UtilitiesForTriangulation::collect_simple_spheres(atoms), bsi_init_radius, false, augment);

	if(!skip_output)
	{
		apollota::TriangulationOutput::print_vertices_vector(apollota::Triangulation::collect_vertices_vector_from_quadruples_map(apollonius_triangulation_result.quadruples_map), std::cout);
	}

	if(print_log)
	{
		std::clog << "atoms " <<  atoms.size() << "\n";
		apollota::TriangulationOutput::print_status(apollonius_triangulation_result, std::clog);
	}

	if(check)
	{
		std::clog << "check " << (apollota::Triangulation::check_quadruples_map(atoms, apollonius_triangulation_result.quadruples_map)) << "\n";
	}
}

}

void calc_quadruples(const auxiliaries::CommandLineOptions& clo)
{
	clo.check_allowed_options("--bsi-init-radius: --raw-input --use-one-radius --augment --skip-output --print-log --check");

	const double bsi_init_radius=clo.isopt("--bsi-init-radius") ? clo.arg_with_min_value<double>("--bsi-radius", 1) : 3.5;
	const bool raw_input=clo.isopt("--raw-input");
	const bool use_one_radius=clo.isopt("--use-one-radius");
	const bool augment=clo.isopt("--augment");
	const bool skip_output=clo.isopt("--skip-output");
	const bool print_log=clo.isopt("--print-log");
	const bool check=clo.isopt("--check");

	if(raw_input)
	{
		calc_quadruples(read_raw_spheres_from_stream(std::cin), bsi_init_radius, use_one_radius, augment, skip_output, print_log, check);
	}
	else
	{
		calc_quadruples(auxiliaries::STDContainersIO::read_vector<protein::Atom>(std::cin, "atoms", "atoms", false), bsi_init_radius, use_one_radius, augment, skip_output, print_log, check);
	}
}
