#ifndef CONTACTO_UTILITIES_H_
#define CONTACTO_UTILITIES_H_

#include <map>
#include <set>
#include <algorithm>
#include <stdexcept>

namespace contacto
{

template<typename Key, typename Value1, typename Value2>
static std::map< Key, std::pair<Value1, Value2> > combine_two_maps(const std::map<Key, Value1>& map_1, const std::map<Key, Value2>& map_2)
{
	typedef std::map< Key, std::pair<Value1, Value2> > CombinedMap;
	CombinedMap combined_map;
	typename std::map<Key, Value1>::const_iterator it_1=map_1.begin();
	typename std::map<Key, Value2>::const_iterator it_2=map_2.begin();
	typename CombinedMap::iterator prev=combined_map.begin();
	while(it_1!=map_1.end() || it_2!=map_2.end())
	{
		if(it_1!=map_1.end() && it_2!=map_2.end() && it_1->first==it_2->first)
		{
			prev=combined_map.insert(prev, std::make_pair(it_1->first, std::make_pair(it_1->second, it_2->second)));
			++it_1;
			++it_2;
		}
		else if(it_1!=map_1.end() && (it_2==map_2.end() || it_1->first<it_2->first))
		{
			prev=combined_map.insert(prev, std::make_pair(it_1->first, std::make_pair(it_1->second, Value2())));
			++it_1;
		}
		else if(it_2!=map_2.end() && (it_1==map_1.end() || it_2->first<it_1->first))
		{
			prev=combined_map.insert(prev, std::make_pair(it_2->first, std::make_pair(Value1(), it_2->second)));
			++it_2;
		}
		else
		{
			throw std::logic_error("Two maps iteration error");
		}
	}
	return combined_map;
}

template<typename T>
static std::set<T> intersect_two_sets(const std::set<T>& set_1, const std::set<T>& set_2)
{
	std::set<T> result;
	typename std::set<T>::const_iterator it_1=set_1.begin();
	typename std::set<T>::const_iterator it_2=set_2.begin();
	typename std::set<T>::iterator prev=result.begin();
	while(it_1!=set_1.end() && it_2!=set_2.end())
	{
		if( (*it_1) < (*it_2) )
		{
			++it_1;
		}
		else if( (*it_2) < (*it_1) )
		{
			++it_2;
		}
		else
		{
			prev=result.insert(prev, (*it_1));
			++it_1;
			++it_2;
		}
	}
	return result;
}

}

#endif /* CONTACTO_UTILITIES_H_ */
