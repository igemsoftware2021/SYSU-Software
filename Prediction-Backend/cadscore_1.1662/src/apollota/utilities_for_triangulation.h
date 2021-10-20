#ifndef APOLLOTA_UTILITIES_FOR_TRIANGULATION_H_
#define APOLLOTA_UTILITIES_FOR_TRIANGULATION_H_

#include "triangulation.h"

namespace apollota
{

class UtilitiesForTriangulation
{
public:
	typedef std::tr1::unordered_map<std::size_t, std::tr1::unordered_set<std::size_t> > NeighborsMap;
	typedef std::vector< std::vector<std::size_t> > NeighborsGraph;
	typedef std::tr1::unordered_map<Pair, std::tr1::unordered_set<std::size_t>, Pair::HashFunctor> PairsNeighborsMap;

	template<typename T>
	static std::vector<SimpleSphere> collect_simple_spheres(const T& spheres)
	{
		std::vector<SimpleSphere> result;
		result.reserve(spheres.size());
		for(typename T::const_iterator it=spheres.begin();it!=spheres.end();++it)
		{
			result.push_back(SimpleSphere(*it));
		}
		return result;
	}

	static NeighborsMap collect_neighbors_map_from_quadruples_map(const Triangulation::QuadruplesMap& quadruples_map)
	{
		NeighborsMap neighbors_map;
		for(Triangulation::QuadruplesMap::const_iterator it=quadruples_map.begin();it!=quadruples_map.end();++it)
		{
			const Quadruple& quadruple=it->first;
			for(int a=0;a<4;a++)
			{
				for(int b=a+1;b<4;b++)
				{
					neighbors_map[quadruple.get(a)].insert(quadruple.get(b));
					neighbors_map[quadruple.get(b)].insert(quadruple.get(a));
				}
			}
		}
		return neighbors_map;
	}

	static NeighborsGraph collect_neighbors_graph_from_neighbors_map(const NeighborsMap& neighbors_map, const std::size_t number_of_vertices)
	{
		NeighborsGraph neighbors_graph(number_of_vertices);
		for(NeighborsMap::const_iterator it=neighbors_map.begin();it!=neighbors_map.end();++it)
		{
			if((it->first)<neighbors_graph.size())
			{
				neighbors_graph[it->first].insert(neighbors_graph[it->first].end(), it->second.begin(), it->second.end());
			}
		}
		return neighbors_graph;
	}

	static PairsNeighborsMap collect_pairs_neighbors_map_from_quadruples_map(const Triangulation::QuadruplesMap& quadruples_map)
	{
		PairsNeighborsMap pairs_neighbors_map;
		for(Triangulation::QuadruplesMap::const_iterator it=quadruples_map.begin();it!=quadruples_map.end();++it)
		{
			const Quadruple& quadruple=it->first;
			for(int a=0;a<4;a++)
			{
				for(int b=a+1;b<4;b++)
				{
					const Pair pair(quadruple.get(a), quadruple.get(b));
					for(int c=0;c<4;c++)
					{
						if(c!=a && c!=b)
						{
							pairs_neighbors_map[pair].insert(quadruple.get(c));
						}
					}
				}
			}
		}
		return pairs_neighbors_map;
	}
};

}

#endif /* APOLLOTA_UTILITIES_FOR_TRIANGULATION_H_ */
