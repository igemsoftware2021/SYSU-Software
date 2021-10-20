#ifndef CONTACTO_INTER_RESIDUE_CONTACTS_FILTERING_H_
#define CONTACTO_INTER_RESIDUE_CONTACTS_FILTERING_H_

#include <vector>
#include <set>
#include <map>
#include <string>
#include <stdexcept>

#include "contact_id.h"

namespace contacto
{

template<typename ResidueIDType, typename ValueType>
std::map< contacto::ContactID<ResidueIDType>, ValueType > filter_core_contacts(const std::map< contacto::ContactID<ResidueIDType>, ValueType >& all_contacts)
{
	typedef std::map< contacto::ContactID<ResidueIDType>, ValueType > ContactsMap;
	ContactsMap core_contacts;
	typename ContactsMap::iterator prev=core_contacts.begin();
	for(typename ContactsMap::const_iterator it=all_contacts.begin();it!=all_contacts.end();++it)
	{
		typename ContactsMap::const_iterator aa_it=all_contacts.find(contacto::ContactID<ResidueIDType>(it->first.a, it->first.a));
		typename ContactsMap::const_iterator bb_it=all_contacts.find(contacto::ContactID<ResidueIDType>(it->first.b, it->first.b));
		if((aa_it==all_contacts.end() || aa_it->second.areas.count("AW")==0) && (bb_it==all_contacts.end() || bb_it->second.areas.count("AW")==0))
		{
			prev=core_contacts.insert(prev, *it);
		}
	}
	return core_contacts;
}

template<typename ResidueIDType, typename ValueType>
std::map< contacto::ContactID<ResidueIDType>, ValueType > filter_interface_zone_contacts(const std::map< contacto::ContactID<ResidueIDType>, ValueType >& all_contacts)
{
	typedef std::map< contacto::ContactID<ResidueIDType>, ValueType > ContactsMap;
	std::set<ResidueIDType> allowed_residue_ids;
	for(typename ContactsMap::const_iterator it=all_contacts.begin();it!=all_contacts.end();++it)
	{
		if(it->first.a.chain_id!=it->first.b.chain_id)
		{
			allowed_residue_ids.insert(it->first.a);
			allowed_residue_ids.insert(it->first.b);
		}
	}
	ContactsMap interface_zone_contacts;
	typename ContactsMap::iterator prev=interface_zone_contacts.begin();
	for(typename ContactsMap::const_iterator it=all_contacts.begin();it!=all_contacts.end();++it)
	{
		if(allowed_residue_ids.count(it->first.a)>0 && allowed_residue_ids.count(it->first.b))
		{
			prev=interface_zone_contacts.insert(prev, *it);
		}
	}
	return interface_zone_contacts;
}

template<typename ResidueIDType, typename ValueType>
std::map< contacto::ContactID<ResidueIDType>, ValueType > filter_inter_chain_contacts(const std::map< contacto::ContactID<ResidueIDType>, ValueType >& all_contacts)
{
	typedef std::map< contacto::ContactID<ResidueIDType>, ValueType > ContactsMap;
	ContactsMap interface_contacts;
	typename ContactsMap::iterator prev=interface_contacts.begin();
	for(typename ContactsMap::const_iterator it=all_contacts.begin();it!=all_contacts.end();++it)
	{
		if(it->first.a.chain_id!=it->first.b.chain_id)
		{
			prev=interface_contacts.insert(prev, *it);
		}
	}
	return interface_contacts;
}

template<typename ResidueIDType, typename ValueType, typename InterfaceReaderType>
std::map< contacto::ContactID<ResidueIDType>, ValueType > filter_inter_interval_contacts(const std::map< contacto::ContactID<ResidueIDType>, ValueType >& all_contacts, const std::string& intervals_string)
{
	typedef std::map< contacto::ContactID<ResidueIDType>, ValueType > ContactsMap;
	std::vector< std::vector< std::pair<ResidueIDType, ResidueIDType> > > intervals;
	if(!InterfaceReaderType::read_residue_ids_intervals(intervals_string, intervals) || intervals.size()<2)
	{
		throw std::runtime_error(std::string("Invalid intervals string: ")+intervals_string);
	}
	ContactsMap interface_contacts;
	typename ContactsMap::iterator prev=interface_contacts.begin();
	for(typename ContactsMap::const_iterator it=all_contacts.begin();it!=all_contacts.end();++it)
	{
		const ResidueIDType& a=it->first.a;
		const ResidueIDType& b=it->first.b;
		int a_group=-1;
		int b_group=-1;
		for(std::size_t i=0;i<intervals.size() && (a_group<0 || b_group<0);i++)
		{
			for(std::size_t j=0;j<intervals[i].size() && (a_group<0 || b_group<0);j++)
			{
				const ResidueIDType& r1=intervals[i][j].first;
				const ResidueIDType& r2=intervals[i][j].second;
				if(a_group<0 && a.chain_id==r1.chain_id && a.residue_number>=r1.residue_number && a.residue_number<=r2.residue_number)
				{
					a_group=i;
				}
				else if(b_group<0 && b.chain_id==r1.chain_id && b.residue_number>=r1.residue_number && b.residue_number<=r2.residue_number)
				{
					b_group=i;
				}
			}
		}
		if(a_group>=0 && b_group>=0 && a_group!=b_group)
		{
			prev=interface_contacts.insert(prev, *it);
		}
	}
	return interface_contacts;
}

template<typename ResidueIDType, typename ValueType, typename InterfaceReaderType>
void filter_custom_contacts(std::map< contacto::ContactID<ResidueIDType>, ValueType >& all_contacts, const bool core, const bool interface_zone, const bool inter_chain, const std::string& intervals_string)
{
	if(core)
	{
		all_contacts=filter_core_contacts<ResidueIDType, ValueType>(all_contacts);
	}

	if(interface_zone)
	{
		all_contacts=filter_interface_zone_contacts<ResidueIDType, ValueType>(all_contacts);
	}

	if(inter_chain)
	{
		all_contacts=filter_inter_chain_contacts<ResidueIDType, ValueType>(all_contacts);
	}

	if(!intervals_string.empty())
	{
		all_contacts=filter_inter_interval_contacts<ResidueIDType, ValueType, InterfaceReaderType>(all_contacts, intervals_string);
	}
}

}

#endif /* CONTACTO_INTER_RESIDUE_CONTACTS_FILTERING_H_ */
