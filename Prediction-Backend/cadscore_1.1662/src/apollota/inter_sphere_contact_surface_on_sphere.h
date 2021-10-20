#ifndef APOLLOTA_INTER_SPHERE_CONTACT_SURFACE_ON_SPHERE_H_
#define APOLLOTA_INTER_SPHERE_CONTACT_SURFACE_ON_SPHERE_H_

#include <vector>
#include <map>

#include "basic_operations_on_spheres.h"
#include "subdivided_icosahedron.h"
#include "hyperboloid_between_two_spheres.h"

namespace apollota
{

class InterSphereContactSurfaceOnSphere
{
public:
	typedef std::map<std::size_t, double> SurfaceArea;

	template<typename SphereType>
	static std::vector<SurfaceArea> calculate_surface_areas(
			const std::vector<SphereType>& spheres,
			const std::vector< std::vector<std::size_t> >& graph,
			const std::size_t subdivision_depth,
			const double probe_radius)
	{
		return construct_surfaces<SurfaceAreaOutputFunctor>(spheres, graph, subdivision_depth, probe_radius);
	}

	template<typename SphereType>
	static double check_surface_area(const SphereType& sphere, const double probe_radius, const SurfaceArea& surface_area)
	{
		const double PI=3.14159265;
		const double theoretical_surface_area=4*PI*(sphere.r+probe_radius)*(sphere.r+probe_radius);
		double sum=0.0;
		for(SurfaceArea::const_iterator it=surface_area.begin();it!=surface_area.end();++it)
		{
			sum+=it->second;
		}
		return (theoretical_surface_area-sum)/theoretical_surface_area;
	}

	template<typename InterSphereContact>
	static std::vector<InterSphereContact> construct_inter_sphere_contacts_from_surface_areas(const std::vector<SurfaceArea>& surface_areas)
	{
		std::vector<InterSphereContact> inter_sphere_contacts;
		std::size_t contacts_count=0;
		for(std::size_t i=0;i<surface_areas.size();i++)
		{
			contacts_count+=surface_areas[i].size();
		}
		inter_sphere_contacts.reserve(contacts_count);
		for(std::size_t i=0;i<surface_areas.size();i++)
		{
			const SurfaceArea& surface_area=surface_areas[i];
			for(SurfaceArea::const_iterator it=surface_area.begin();it!=surface_area.end();it++)
			{
				inter_sphere_contacts.push_back(InterSphereContact(i, it->first, it->second));
			}
		}
		return inter_sphere_contacts;
	}

private:
	struct SurfaceAreaOutputFunctor
	{
		typedef SurfaceArea ResultType;
		ResultType result;

		void operator()(const std::size_t id, const SimplePoint& a, const SimplePoint& b, const SimplePoint& c)
		{
			result[id]+=triangle_area(a, b, c);
		}
	};

	template<typename OutputFunctor, typename SphereType>
	static std::vector<typename OutputFunctor::ResultType> construct_surfaces(
			const std::vector<SphereType>& spheres,
			const std::vector< std::vector<std::size_t> >& graph,
			const std::size_t subdivision_depth,
			const double probe_radius)
	{
		std::vector<typename OutputFunctor::ResultType> surfaces;
		SubdividedIcosahedron sih(subdivision_depth);
		surfaces.reserve(graph.size());
		for(std::size_t i=0;i<spheres.size();i++)
		{
			sih.fit_into_sphere(spheres[i], spheres[i].r+probe_radius);
			OutputFunctor output_functor;
			construct_surface(
					sih,
					spheres,
					collect_influences(sih, spheres, i, graph[i]),
					output_functor);
			surfaces.push_back(output_functor.result);
		}
		return surfaces;
	}

	template<typename SphereType>
	static std::vector<std::size_t> collect_influences(
			const SubdividedIcosahedron& sih,
			const std::vector<SphereType>& spheres,
			const std::size_t self_id,
			const std::vector<size_t>& neighbours)
	{
		std::vector<std::size_t> influences(sih.vertices().size());
		for(std::size_t i=0;i<influences.size();i++)
		{
			double min_distance=minimal_distance_from_point_to_sphere(sih.vertices().at(i), spheres.at(self_id));
			influences[i]=self_id;
			for(std::size_t j=0;j<neighbours.size();j++)
			{
				double distance=minimal_distance_from_point_to_sphere(sih.vertices().at(i), spheres.at(neighbours[j]));
				if(distance<min_distance)
				{
					min_distance=distance;
					influences[i]=neighbours[j];
				}
			}
		}
		return influences;
	}

	template<typename SphereType, typename OutputFunctor>
	static void construct_surface(
			const SubdividedIcosahedron& sih,
			const std::vector<SphereType>& spheres,
			const std::vector<std::size_t>& influences,
			OutputFunctor& output_functor)
	{
		for(std::size_t e=0;e<sih.triples().size();e++)
		{
			const Triple& triple=sih.triples()[e];
			const std::size_t a=triple.get(0);
			const std::size_t b=triple.get(1);
			const std::size_t c=triple.get(2);
			if(influences[a]==influences[b] && influences[a]==influences[c])
			{
				output_functor(influences[a], sih.vertices()[a], sih.vertices()[b], sih.vertices()[c]);
			}
			else if(influences[a]!=influences[b] && influences[a]!=influences[c] && influences[b]!=influences[c])
			{
				const SimplePoint& pa=sih.vertices()[a];
				const SimplePoint& pb=sih.vertices()[b];
				const SimplePoint& pc=sih.vertices()[c];

				const SimplePoint a_b_border=pa+((pb-pa).unit()*HyperboloidBetweenTwoSpheres::intersect_vector_with_hyperboloid(pa, pb, spheres[influences[a]], spheres[influences[b]]));
				const SimplePoint a_c_border=pa+((pc-pa).unit()*HyperboloidBetweenTwoSpheres::intersect_vector_with_hyperboloid(pa, pc, spheres[influences[a]], spheres[influences[c]]));
				const SimplePoint b_c_border=pb+((pc-pb).unit()*HyperboloidBetweenTwoSpheres::intersect_vector_with_hyperboloid(pb, pc, spheres[influences[b]], spheres[influences[c]]));

				const SimplePoint middle=(a_b_border+a_c_border+b_c_border)*(1.0/3.0);

				output_functor(influences[a], pa, a_b_border, a_c_border);
				output_functor(influences[a], middle, a_b_border, a_c_border);

				output_functor(influences[b], pb, a_b_border, b_c_border);
				output_functor(influences[b], middle, a_b_border, b_c_border);

				output_functor(influences[c], pc, a_c_border, b_c_border);
				output_functor(influences[c], middle, a_c_border, b_c_border);
			}
			else
			{
				std::size_t s=a;
				std::size_t d1=b;
				std::size_t d2=c;
				if(influences[b]!=influences[a] && influences[b]!=influences[c])
				{
					s=b;
					d1=a;
					d2=c;
				}
				else if(influences[c]!=influences[a] && influences[c]!=influences[b])
				{
					s=c;
					d1=a;
					d2=b;
				}

				const SimplePoint& ps=sih.vertices()[s];
				const SimplePoint& pd1=sih.vertices()[d1];
				const SimplePoint& pd2=sih.vertices()[d2];

				const SimplePoint s_d1_border=ps+((pd1-ps).unit()*HyperboloidBetweenTwoSpheres::intersect_vector_with_hyperboloid(ps, pd1, spheres[influences[s]], spheres[influences[d1]]));
				const SimplePoint s_d2_border=ps+((pd2-ps).unit()*HyperboloidBetweenTwoSpheres::intersect_vector_with_hyperboloid(ps, pd2, spheres[influences[s]], spheres[influences[d2]]));

				output_functor(influences[s], ps, s_d1_border, s_d2_border);
				output_functor(influences[d1], pd1, s_d1_border, s_d2_border);
				output_functor(influences[d2], pd2, pd1, s_d2_border);
			}
		}
	}
};

}

#endif /* APOLLOTA_INTER_SPHERE_CONTACT_SURFACE_ON_SPHERE_H_ */
