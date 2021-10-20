#include <iostream>
#include <sstream>
#include <iomanip>
#include <memory>

#include "protein/atom.h"
#include "protein/residue_id.h"
#include "protein/residue_ids_intervals.h"

#include "contacto/contact_id.h"
#include "contacto/inter_residue_contact_dual_areas.h"
#include "contacto/contact_classification.h"

#include "apollota/triangulation.h"
#include "apollota/utilities_for_triangulation.h"
#include "apollota/inter_sphere_contact_face_on_hyperboloid.h"
#include "apollota/opengl_printer.h"

#include "auxiliaries/command_line_options.h"
#include "auxiliaries/std_containers_io.h"
#include "auxiliaries/color.h"

namespace
{

class ColorManagementForMapping
{
public:
	static auxiliaries::Color default_color()
	{
		return auxiliaries::Color::from_code(0xFFFFFF);
	}

	template<typename T>
	static auxiliaries::Color color_from_map(const std::map<T, auxiliaries::Color>& map_of_colors, const T& name)
	{
		typename std::map<T, auxiliaries::Color>::const_iterator it=map_of_colors.find(name);
		return (it==map_of_colors.end() ? default_color() : it->second);
	}
};

template<typename T>
class NameColorizer
{
public:
	NameColorizer()
	{
	}

	const std::map<T, auxiliaries::Color>& map_of_colors() const
	{
		return map_of_colors_;
	}

	auxiliaries::Color color(const T& name) const
	{
		return ColorManagementForMapping::color_from_map(map_of_colors_, name);
	}

	void set_map_of_colors(const std::map<T, auxiliaries::Color> map_of_colors)
	{
		map_of_colors_=map_of_colors;
	}

	void add_name_color(const T& name, const auxiliaries::Color& color)
	{
		map_of_colors_[name]=color;
	}

private:
	std::map<T, auxiliaries::Color> map_of_colors_;
};

class ColorManagementForPymol
{
public:
	static std::string color_to_string_id(const auxiliaries::Color& color)
	{
		std::ostringstream output;
		output << "custom_color_" << static_cast<int>(color.r) << "_" << static_cast<int>(color.g) << "_" << static_cast<int>(color.b);
		return output.str();
	}

	static std::string color_to_string_value(const auxiliaries::Color& color)
	{
		std::ostringstream output;
		output << "[ " << color.r_double() << ", " << color.g_double() << ", " << color.b_double() << " ]";
		return output.str();
	}

	static void list_color(const auxiliaries::Color& color)
	{
		std::cout << "cmd.set_color('" << color_to_string_id(color) << "', " << color_to_string_value(color) << ")\n";
	}

	template<typename ColorsMapType>
	static void list_colors_from_map(const ColorsMapType& map_of_colors)
	{
		for(typename ColorsMapType::const_iterator it=map_of_colors.begin();it!=map_of_colors.end();++it)
		{
			list_color(it->second);
		}
		list_color(ColorManagementForMapping::default_color());
		std::cout << "\n";
	}
};

template<typename T>
class NameColorizerForPymol : public NameColorizer<T>
{
public:
	NameColorizerForPymol()
	{
	}

	std::string color_string(const T& name) const
	{
		return ColorManagementForPymol::color_to_string_id(color(name));
	}

	void list_colors() const
	{
		ColorManagementForPymol::list_colors_from_map(NameColorizer<T>::map_of_colors());
	}
};

class ResidueNameColorizerByResidueType : public NameColorizerForPymol<std::string>
{
public:
	ResidueNameColorizerByResidueType()
	{
		set_map_of_colors(create_map_of_residue_colors_by_type());
	}

private:
	static std::map<std::string, auxiliaries::Color> create_map_of_residue_colors_by_type()
	{
		const auxiliaries::Color nonpolar(255, 255, 0);
		const auxiliaries::Color acidic(255, 0, 0);
		const auxiliaries::Color basic(0, 0, 255);
		const auxiliaries::Color uncharged(0, 255, 0);

		std::map<std::string, auxiliaries::Color> m;

		m["LEU"]=nonpolar;
		m["VAL"]=nonpolar;
		m["ILE"]=nonpolar;
		m["ALA"]=nonpolar;
		m["PHE"]=nonpolar;
		m["TRP"]=nonpolar;
		m["MET"]=nonpolar;
		m["PRO"]=nonpolar;

		m["ASP"]=acidic;
		m["GLU"]=acidic;

		m["LYS"]=basic;
		m["ARG"]=basic;
		m["HIS"]=basic;

		m["CYS"]=uncharged;
		m["SER"]=uncharged;
		m["THR"]=uncharged;
		m["TYR"]=uncharged;
		m["ASN"]=uncharged;
		m["GLN"]=uncharged;
		m["GLY"]=uncharged;

		return m;
	}
};

class ResidueNameColorizerByResidueHydrophobicity : public NameColorizerForPymol<std::string>
{
public:
	ResidueNameColorizerByResidueHydrophobicity()
	{
		set_map_of_colors(create_map_of_residue_colors_by_hydropathy_indices());
	}

private:
	static auxiliaries::Color color_from_hydropathy_index(const double hi)
	{
		return auxiliaries::Color::from_temperature_to_blue_white_red((1+hi/4.5)/2);
	}

	static std::map<std::string, auxiliaries::Color> create_map_of_residue_colors_by_hydropathy_indices()
	{
		std::map<std::string, auxiliaries::Color> m;
		m["ASP"]=color_from_hydropathy_index(-3.5);
		m["GLU"]=color_from_hydropathy_index(-3.5);
		m["CYS"]=color_from_hydropathy_index(2.5);
		m["MET"]=color_from_hydropathy_index(1.9);
		m["LYS"]=color_from_hydropathy_index(-3.9);
		m["ARG"]=color_from_hydropathy_index(-4.5);
		m["SER"]=color_from_hydropathy_index(-0.8);
		m["THR"]=color_from_hydropathy_index(-0.7);
		m["PHE"]=color_from_hydropathy_index(2.8);
		m["TYR"]=color_from_hydropathy_index(-1.3);
		m["ASN"]=color_from_hydropathy_index(-3.5);
		m["GLN"]=color_from_hydropathy_index(-3.5);
		m["GLY"]=color_from_hydropathy_index(-0.4);
		m["LEU"]=color_from_hydropathy_index(3.8);
		m["VAL"]=color_from_hydropathy_index(4.2);
		m["ILE"]=color_from_hydropathy_index(4.5);
		m["ALA"]=color_from_hydropathy_index(1.8);
		m["TRP"]=color_from_hydropathy_index(-0.9);
		m["HIS"]=color_from_hydropathy_index(-3.2);
		m["PRO"]=color_from_hydropathy_index(-1.6);
		return m;
	}
};

class AtomNameColorizerByAtomType : public NameColorizerForPymol<std::string>
{
public:
	AtomNameColorizerByAtomType()
	{
		set_map_of_colors(create_map_of_atom_colors_by_type());
	}

private:
	static std::map<std::string, auxiliaries::Color> create_map_of_atom_colors_by_type()
	{
		std::map<std::string, auxiliaries::Color> m;

		m["C"]=auxiliaries::Color(0, 255, 0);
		m["N"]=auxiliaries::Color(0, 0, 255);
		m["O"]=auxiliaries::Color(255, 0, 0);
		m["S"]=auxiliaries::Color(255, 255, 0);
		m["P"]=auxiliaries::Color(255, 0, 255);

		return m;
	}
};

class ValueColorizer : public NameColorizerForPymol<int>
{
public:
	ValueColorizer()
	{
		set_map_of_colors(create_map_of_colors_by_values());
	}

private:
	static std::map<int, auxiliaries::Color> create_map_of_colors_by_values()
	{
		std::map<int, auxiliaries::Color> m;
		for(int i=0;i<=100;i++)
		{
			m[i]=auxiliaries::Color::from_temperature_to_blue_white_red(static_cast<double>(i)/100.0);
		}
		m[-1]=auxiliaries::Color(255, 0, 255);
		return m;
	}
};

class ContactColorizerInterface
{
public:
	virtual auxiliaries::Color color(const protein::Atom& a, const protein::Atom& b) const = 0;

	std::string color_string(const protein::Atom& a, const protein::Atom& b) const
	{
		return ColorManagementForPymol::color_to_string_id(color(a, b));
	}

	virtual void list_colors() const = 0;
};

class ContactColorizerByCustomSingleColor : public ContactColorizerInterface
{
public:
	ContactColorizerByCustomSingleColor(const std::string& color_code_string) : color_(255, 255, 255)
	{
		std::istringstream input(color_code_string);
		if(input.good())
		{
			int color_code=0xFFFFFF;
			input >> std::hex >> color_code;
			if(!input.fail() && color_code>0)
			{
				color_=auxiliaries::Color::from_code(color_code);
			}
		}
	}

	auxiliaries::Color color(const protein::Atom& /*a*/, const protein::Atom& /*b*/) const
	{
		return color_;
	}

	virtual void list_colors() const
	{
		ColorManagementForPymol::list_color(color_);
	}

private:
	auxiliaries::Color color_;
};

template<class ResidueNameColorizerType>
class ContactColorizerByFirstResidueName : public ContactColorizerInterface
{
public:
	ContactColorizerByFirstResidueName()
	{
	}

	auxiliaries::Color color(const protein::Atom& a, const protein::Atom& /*b*/) const
	{
		return name_colorizer_.color(a.residue_name);
	}

	virtual void list_colors() const
	{
		name_colorizer_.list_colors();
	}

private:
	ResidueNameColorizerType name_colorizer_;
};

class ContactColorizerByFirstAtomName : public ContactColorizerInterface
{
public:
	ContactColorizerByFirstAtomName()
	{
	}

	auxiliaries::Color color(const protein::Atom& a, const protein::Atom& /*b*/) const
	{
		return name_colorizer_.color(a.atom_name.substr(0, 1));
	}

	virtual void list_colors() const
	{
		name_colorizer_.list_colors();
	}

private:
	AtomNameColorizerByAtomType name_colorizer_;
};

class ContactColorizerByInterResidueContactScore : public ContactColorizerInterface
{
public:
	ContactColorizerByInterResidueContactScore(
			const std::map< contacto::ContactID<protein::ResidueID>, contacto::InterResidueContactDualAreas >& combined_inter_residue_contacts,
			const std::string& specific_contact_type,
			const bool binary_mode) :
				combined_inter_residue_contacts_(combined_inter_residue_contacts),
				specific_contact_type_(specific_contact_type),
				binary_mode_(binary_mode)
	{
	}

	auxiliaries::Color color(const protein::Atom& a, const protein::Atom& b) const
	{
		const contacto::ContactID<protein::ResidueID> irc_id(protein::ResidueID::from_atom(a), protein::ResidueID::from_atom(b));
		std::map< contacto::ContactID<protein::ResidueID>, contacto::InterResidueContactDualAreas >::const_iterator it=combined_inter_residue_contacts_.find(irc_id);
		if(it!=combined_inter_residue_contacts_.end())
		{
			std::string contact_type=specific_contact_type_;
			if(contact_type.empty())
			{
				const std::vector<std::string> contact_types=contacto::ContactClassification::classify_atoms_contact<protein::Atom, protein::ResidueID>(a, b);
				if(!contact_types.empty())
				{
					contact_type=contact_types.front();
				}
			}
			std::pair<double, double> area=it->second.area(contact_type);
			if(binary_mode_)
			{
				return value_colorizer_.color((area.first>0.0) ? 0 : 100);
			}
			else
			{
				if(area.first>0.0)
				{
					return value_colorizer_.color(static_cast<int>(floor((std::min(fabs(area.first-area.second), area.first)/area.first)*100.0)));
				}
			}
		}
		return value_colorizer_.color(-1);
	}

	virtual void list_colors() const
	{
		value_colorizer_.list_colors();
	}

private:
	std::map< contacto::ContactID<protein::ResidueID>, contacto::InterResidueContactDualAreas > combined_inter_residue_contacts_;
	std::string specific_contact_type_;
	bool binary_mode_;
	ValueColorizer value_colorizer_;
};

class ContactColorizerByFirstResidueID : public ContactColorizerInterface
{
public:
	ContactColorizerByFirstResidueID(const std::vector<protein::Atom>& atoms)
	{
		for(std::size_t i=0;i<atoms.size();i++)
		{
			const long numeric_rid=numeric_residue_id(protein::ResidueID::from_atom(atoms[i]));
			name_colorizer_.add_name_color(numeric_rid, auxiliaries::Color::from_id(numeric_rid));
		}
	}

	auxiliaries::Color color(const protein::Atom& a, const protein::Atom& /*b*/) const
	{
		return name_colorizer_.color(numeric_residue_id(protein::ResidueID::from_atom(a)));
	}

	virtual void list_colors() const
	{
		name_colorizer_.list_colors();
	}

private:
	static long numeric_residue_id(const protein::ResidueID& residue_id)
	{
		return static_cast<long>(residue_id.residue_number)+(residue_id.chain_id.empty() ? 0 : static_cast<long>(residue_id.chain_id.c_str()[0])*1000000);
	}

	NameColorizerForPymol<long> name_colorizer_;
};

class ContactColorizerByFirstAtomID : public ContactColorizerInterface
{
public:
	ContactColorizerByFirstAtomID(const std::vector<protein::Atom>& atoms)
	{
		for(std::size_t i=0;i<atoms.size();i++)
		{
			const protein::Atom& a=atoms[i];
			name_colorizer_.add_name_color(a.atom_number, auxiliaries::Color::from_id(a.atom_number));
		}
	}

	auxiliaries::Color color(const protein::Atom& a, const protein::Atom& /*b*/) const
	{
		return name_colorizer_.color(a.atom_number);
	}

	virtual void list_colors() const
	{
		name_colorizer_.list_colors();
	}

private:
	NameColorizerForPymol<int> name_colorizer_;
};

class ContactAccepterInterface
{
public:
	virtual bool accept(const protein::Atom& a, const protein::Atom& b) const = 0;
	virtual std::string assign_group_name(const protein::Atom& a) const = 0;
};

class ContactAccepterForInterChain : public ContactAccepterInterface
{
public:
	bool accept(const protein::Atom& a, const protein::Atom& b) const
	{
		return (a.chain_id!=b.chain_id && a.chain_id!="?" && b.chain_id!="?");
	}

	std::string assign_group_name(const protein::Atom& a) const
	{
		return a.chain_id;
	}
};

class ContactAccepterForInsideSet : public ContactAccepterInterface
{
public:
	ContactAccepterForInsideSet(const std::set<protein::ResidueID>& allowed_residue_ids) : allowed_residue_ids_(allowed_residue_ids)
	{
	}

	bool accept(const protein::Atom& a, const protein::Atom& b) const
	{
		return (protein::ResidueID::from_atom(a)!=protein::ResidueID::from_atom(b)&&
				allowed_residue_ids_.count(protein::ResidueID::from_atom(a))>0 &&
				allowed_residue_ids_.count(protein::ResidueID::from_atom(b))>0);
	}

	std::string assign_group_name(const protein::Atom& a) const
	{
		return a.chain_id;
	}

private:
	std::set<protein::ResidueID> allowed_residue_ids_;
};

class ContactAccepterForInterResidue : public ContactAccepterInterface
{
public:
	ContactAccepterForInterResidue(bool only_side_chain_contacts) : only_side_chain_contacts_(only_side_chain_contacts)
	{
	}

	bool accept(const protein::Atom& a, const protein::Atom& b) const
	{
		return (protein::ResidueID::from_atom(a)!=protein::ResidueID::from_atom(b) &&
				(!only_side_chain_contacts_ || (a.location_class==protein::Atom::side_chain && b.location_class==protein::Atom::side_chain)));
	}

	std::string assign_group_name(const protein::Atom& /*a*/) const
	{
		return std::string("residues");
	}

private:
	bool only_side_chain_contacts_;
};

class ContactAccepterForInterInterval : public ContactAccepterInterface
{
public:
	ContactAccepterForInterInterval(
			const std::vector< std::vector< std::pair<protein::ResidueID, protein::ResidueID> > >& intervals,
			bool only_side_chain_contacts,
			bool accept_contacts_inside_residues) :
				intervals_(intervals),
				only_side_chain_contacts_(only_side_chain_contacts),
				accept_contacts_inside_residues_(accept_contacts_inside_residues)
	{
	}

	bool accept(const protein::Atom& a, const protein::Atom& b) const
	{
		if((!only_side_chain_contacts_ || (a.location_class==protein::Atom::side_chain && b.location_class==protein::Atom::side_chain)))
		{
			const int a_iid=interval_id(a);
			const int b_iid=interval_id(b);
			if(a_iid!=b_iid)
			{
				if(a_iid>=0 && b_iid>=0)
				{
					return true;
				}
				else if(intervals_.size()==1)
				{
					return true;
				}
			}
			else if(a_iid>=0 && accept_contacts_inside_residues_)
			{
				return true;
			}
		}
		return false;
	}

	std::string assign_group_name(const protein::Atom& a) const
	{
		const int iid=interval_id(a);
		if(iid>=0)
		{
			std::ostringstream output;
			output << iid;
			return output.str();
		}
		else
		{
			return "all";
		}
	}

private:
	int interval_id(const protein::Atom& a) const
	{
		for(std::size_t i=0;i<intervals_.size();i++)
		{
			for(std::size_t j=0;j<intervals_[i].size();j++)
			{
				const protein::ResidueID& r1=intervals_[i][j].first;
				const protein::ResidueID& r2=intervals_[i][j].second;
				if(a.chain_id==r1.chain_id && a.residue_number>=r1.residue_number && a.residue_number<=r2.residue_number)
				{
					return i;
				}
			}
		}
		return (-1);
	}

	std::vector< std::vector< std::pair<protein::ResidueID, protein::ResidueID> > > intervals_;
	bool only_side_chain_contacts_;
	bool accept_contacts_inside_residues_;
};

std::string atom_name_without_single_quote(const std::string& full_atom_name)
{
	const std::size_t quote_pos=full_atom_name.find('\'', 0);
	if(quote_pos==std::string::npos)
	{
		return full_atom_name;
	}
	else if(quote_pos==0)
	{
		return (full_atom_name.size()>1 ? full_atom_name.substr(1) : std::string("Q"));
	}
	else
	{
		return full_atom_name.substr(0, quote_pos);
	}
}

}

void print_inter_chain_interface_graphics(const auxiliaries::CommandLineOptions& clo)
{
	typedef apollota::InterSphereContactFaceOnHyperboloid CellFace;

	clo.check_allowed_options("--probe: --step: --projections: --face-coloring: --selection-coloring: --groups: --output-names-prefix: --outline --insides --specific-contact-type: --transparent-magenta --binary-coloring");

	const double probe_radius=clo.isopt("--probe") ? clo.arg_with_min_value<double>("--probe", 0) : 1.4;
	const double step_length=clo.isopt("--step") ? clo.arg_with_min_value<double>("--step", 0.1) : 0.5;
	const int projections_count=clo.isopt("--projections") ? clo.arg_with_min_value<int>("--projections", 5) : 5;
	const std::string face_coloring_mode=clo.isopt("--face-coloring") ? clo.arg<std::string>("--face-coloring") : std::string("");
	const std::string selection_coloring_mode=clo.isopt("--selection-coloring") ? clo.arg<std::string>("--selection-coloring") : std::string("");
	const std::string groups_option=clo.isopt("--groups") ? clo.arg<std::string>("--groups") : std::string("");
	const std::string output_names_prefix=clo.isopt("--output-names-prefix") ? clo.arg<std::string>("--output-names-prefix") : std::string("");
	const bool outline=clo.isopt("--outline");
	const bool insides=clo.isopt("--insides");
	const std::string specific_contact_type=clo.isopt("--specific-contact-type") ? clo.arg<std::string>("--specific-contact-type") : std::string("");
	const bool tansparent_magenta=clo.isopt("--transparent-magenta");
	const bool binary_coloring=clo.isopt("--binary-coloring");

	const std::vector<protein::Atom> atoms=auxiliaries::STDContainersIO::read_vector<protein::Atom>(std::cin, "atoms", "atoms", false);

	const apollota::UtilitiesForTriangulation::PairsNeighborsMap pairs_neighbours_map=apollota::UtilitiesForTriangulation::collect_pairs_neighbors_map_from_quadruples_map(apollota::Triangulation::construct_result(apollota::UtilitiesForTriangulation::collect_simple_spheres(atoms), 3.5, false, false).quadruples_map);

#ifdef FOR_OLDER_COMPILERS
	typedef std::auto_ptr<ContactAccepterInterface> AutoPtrToContactAccepterInterface;
#else
	typedef std::unique_ptr<ContactAccepterInterface> AutoPtrToContactAccepterInterface;
#endif

	AutoPtrToContactAccepterInterface contact_accepter;
	if(groups_option.substr(0,1)=="(")
	{
		std::string intervals_string=groups_option;
		bool only_sidechains=false;
		std::size_t spec_pos=groups_option.find("_SS");
		if(spec_pos!=std::string::npos)
		{
			intervals_string=groups_option.substr(0, spec_pos);
			only_sidechains=true;
		}
		std::vector< std::vector< std::pair<protein::ResidueID, protein::ResidueID> > > intervals;
		if(!protein::ResidueIDsIntervalsReader::read_residue_ids_intervals(intervals_string, intervals) || intervals.empty())
		{
			throw std::runtime_error(std::string("Invalid intervals string: ")+intervals_string);
		}
		contact_accepter.reset(new ContactAccepterForInterInterval(intervals, only_sidechains, insides));
	}
	else if(groups_option=="inter_residue" || groups_option=="inter_residue_SS")
	{
		contact_accepter.reset(new ContactAccepterForInterResidue(groups_option=="inter_residue_SS"));
	}
	else if(groups_option=="inter_chain_region")
	{
		std::set<protein::ResidueID> region_residue_ids;
		for(apollota::UtilitiesForTriangulation::PairsNeighborsMap::const_iterator it=pairs_neighbours_map.begin();it!=pairs_neighbours_map.end();++it)
		{
			const std::pair<std::size_t, std::size_t> atoms_ids_pair=std::make_pair(it->first.get(0), it->first.get(1));
			const protein::Atom& a=atoms[atoms_ids_pair.first];
			const protein::Atom& b=atoms[atoms_ids_pair.second];
			if(a.chain_id!=b.chain_id && a.chain_id!="?" && b.chain_id!="?" &&
					apollota::minimal_distance_from_sphere_to_sphere(a, b)<probe_radius*2)
			{
				region_residue_ids.insert(protein::ResidueID::from_atom(a));
				region_residue_ids.insert(protein::ResidueID::from_atom(b));
			}
		}
		contact_accepter.reset(new ContactAccepterForInsideSet(region_residue_ids));
	}
	else
	{
		contact_accepter.reset(new ContactAccepterForInterChain());
	}

	std::vector<CellFace> faces_vector;
	std::map< std::pair<std::size_t, std::size_t>, std::size_t > faces_vector_map;
	typedef std::map< std::pair<std::string, std::string>, std::vector< std::pair<std::size_t, std::size_t> > > InterfacesMap;
	InterfacesMap inter_chain_interfaces;

	for(apollota::UtilitiesForTriangulation::PairsNeighborsMap::const_iterator it=pairs_neighbours_map.begin();it!=pairs_neighbours_map.end();++it)
	{
		const std::pair<std::size_t, std::size_t> atoms_ids_pair=std::make_pair(it->first.get(0), it->first.get(1));

		const protein::Atom& a=atoms[atoms_ids_pair.first];
		const protein::Atom& b=atoms[atoms_ids_pair.second];

		if(contact_accepter->accept(a,b))
		{
			std::vector<const protein::Atom*> cs;
			cs.reserve(it->second.size());
			for(apollota::UtilitiesForTriangulation::PairsNeighborsMap::mapped_type::const_iterator jt=it->second.begin();jt!=it->second.end();++jt)
			{
				cs.push_back(&(atoms[*jt]));
			}

			const CellFace cell_face=CellFace::construct(a, b, cs, probe_radius, step_length, projections_count);
			if(!cell_face.mesh_vertices().empty())
			{
				const std::pair<std::size_t, std::size_t> reversed_atoms_ids_pair=std::make_pair(atoms_ids_pair.second, atoms_ids_pair.first);
				faces_vector.push_back(cell_face);
				faces_vector_map[atoms_ids_pair]=faces_vector.size()-1;
				faces_vector_map[reversed_atoms_ids_pair]=faces_vector.size()-1;
				inter_chain_interfaces[std::make_pair(contact_accepter->assign_group_name(a), contact_accepter->assign_group_name(b))].push_back(atoms_ids_pair);
				inter_chain_interfaces[std::make_pair(contact_accepter->assign_group_name(b), contact_accepter->assign_group_name(a))].push_back(reversed_atoms_ids_pair);
			}
		}
	}

	if(inter_chain_interfaces.empty())
	{
		throw std::runtime_error("No interfaces found");
	}

#ifdef FOR_OLDER_COMPILERS
	typedef std::auto_ptr<const ContactColorizerInterface> AutoPtrToContactColorizerInterface;
#else
	typedef std::unique_ptr<const ContactColorizerInterface> AutoPtrToContactColorizerInterface;
#endif

	AutoPtrToContactColorizerInterface face_colorizer;
	if(face_coloring_mode=="residue_type")
	{
		face_colorizer.reset(new ContactColorizerByFirstResidueName<ResidueNameColorizerByResidueType>());
	}
	else if(face_coloring_mode=="residue_hydrophobicity")
	{
		face_colorizer.reset(new ContactColorizerByFirstResidueName<ResidueNameColorizerByResidueHydrophobicity>());
	}
	else if(face_coloring_mode=="atom_type")
	{
		face_colorizer.reset(new ContactColorizerByFirstAtomName());
	}
	else if(face_coloring_mode=="inter_residue_contact_scores")
	{
		const std::map< contacto::ContactID<protein::ResidueID>, contacto::InterResidueContactDualAreas > combined_inter_residue_contacts=
				auxiliaries::STDContainersIO::read_map< contacto::ContactID<protein::ResidueID>, contacto::InterResidueContactDualAreas >(std::cin, "combined inter-residue contacts", "combined_residue_contacts", true);
		face_colorizer.reset(new ContactColorizerByInterResidueContactScore(combined_inter_residue_contacts, specific_contact_type, binary_coloring));
	}
	else if(face_coloring_mode=="residue_id")
	{
		face_colorizer.reset(new ContactColorizerByFirstResidueID(atoms));
	}
	else if(face_coloring_mode=="atom_id")
	{
		face_colorizer.reset(new ContactColorizerByFirstAtomID(atoms));
	}
	else if(face_coloring_mode.substr(0,2)=="0x")
	{
		face_colorizer.reset(new ContactColorizerByCustomSingleColor(face_coloring_mode));
	}
	else
	{
		face_colorizer.reset(new ContactColorizerByFirstResidueName< NameColorizerForPymol<std::string> >());
	}

	apollota::OpenGLPrinter::print_setup(std::cout);

	for(InterfacesMap::const_iterator it=inter_chain_interfaces.begin();it!=inter_chain_interfaces.end();++it)
	{
		const std::string obj_name=output_names_prefix+std::string("obj_")+it->first.first+"_"+it->first.second;
		const std::string cgo_name=output_names_prefix+std::string("iface_")+it->first.first+"_"+it->first.second;
		apollota::OpenGLPrinter opengl_printer(std::cout, obj_name, cgo_name);
		for(std::size_t i=0;i<it->second.size();++i)
		{
			const std::pair<std::size_t, std::size_t> atoms_ids_pair=it->second[i];
			const protein::Atom& a=atoms[atoms_ids_pair.first];
			const protein::Atom& b=atoms[atoms_ids_pair.second];
			const CellFace& cell_face=faces_vector[faces_vector_map.find(atoms_ids_pair)->second];
			const auxiliaries::Color face_color=face_colorizer->color(a, b);
			if(!(tansparent_magenta && face_color.r==255 && face_color.g==0 && face_color.b==255))
			{
				opengl_printer.print_color(face_color.r_double(), face_color.g_double(), face_color.b_double());
				if(outline)
				{
					std::vector<apollota::SimplePoint> loop_vertices=cell_face.mesh_vertices();
					loop_vertices.pop_back();
					opengl_printer.print_line_strip(loop_vertices, true);
				}
				else
				{
					const apollota::SimplePoint normal=apollota::sub_of_points<apollota::SimplePoint>(b, a).unit();
					const apollota::SimplePoint shift=normal*0.001;
					std::vector<apollota::SimplePoint> fan_vertices;
					fan_vertices.reserve(cell_face.mesh_vertices().size()+1);
					fan_vertices.push_back(cell_face.mesh_vertices().back()+shift);
					for(std::size_t j=0;j+1<cell_face.mesh_vertices().size();j++)
					{
						fan_vertices.push_back(cell_face.mesh_vertices()[j]+shift);
					}
					fan_vertices.push_back(cell_face.mesh_vertices().front()+shift);
					opengl_printer.print_triangle_strip(fan_vertices, std::vector<apollota::SimplePoint>(fan_vertices.size(), normal), true);
				}
			}
		}
	}

	if(!selection_coloring_mode.empty())
	{
		bool color_pymol_selection_at_atomic_level=true;
		AutoPtrToContactColorizerInterface selection_colorizer;
		if(selection_coloring_mode=="residue_type")
		{
			selection_colorizer.reset(new ContactColorizerByFirstResidueName<ResidueNameColorizerByResidueType>());
			color_pymol_selection_at_atomic_level=false;
		}
		else if(selection_coloring_mode=="residue_hydrophobicity")
		{
			selection_colorizer.reset(new ContactColorizerByFirstResidueName<ResidueNameColorizerByResidueHydrophobicity>());
			color_pymol_selection_at_atomic_level=false;
		}
		else if(selection_coloring_mode=="atom_type")
		{
			selection_colorizer.reset(new ContactColorizerByFirstAtomName());
		}
		else if(selection_coloring_mode=="residue_id")
		{
			selection_colorizer.reset(new ContactColorizerByFirstResidueID(atoms));
			color_pymol_selection_at_atomic_level=false;
		}
		else if(selection_coloring_mode=="atom_id")
		{
			selection_colorizer.reset(new ContactColorizerByFirstAtomID(atoms));
		}
		else if(selection_coloring_mode.substr(0,2)=="0x")
		{
			selection_colorizer.reset(new ContactColorizerByCustomSingleColor(selection_coloring_mode));
		}
		else
		{
			selection_colorizer.reset(new ContactColorizerByFirstResidueName< NameColorizerForPymol<std::string> >());
			color_pymol_selection_at_atomic_level=false;
		}

		std::cout << "cmd.color('gray')\n\n";
		std::cout << "cmd.hide('nonbonded')\n\n";
		std::cout << "cmd.hide('lines')\n\n";

		selection_colorizer->list_colors();

		for(InterfacesMap::const_iterator it=inter_chain_interfaces.begin();it!=inter_chain_interfaces.end();++it)
		{
			const std::string selection_name=output_names_prefix+std::string("iface_sel_")+it->first.first+"_"+it->first.second;

			std::map< protein::ResidueID, std::vector< std::pair<std::size_t, std::size_t> > > selectable_residue_ids;
			for(std::size_t i=0;i<it->second.size();++i)
			{
				const std::pair<std::size_t, std::size_t>& atoms_ids_pair=it->second[i];
				const protein::Atom& a=atoms[atoms_ids_pair.first];
				selectable_residue_ids[protein::ResidueID::from_atom(a)].push_back(atoms_ids_pair);
			}

			std::cout << "cmd.select('" << selection_name << "', '";
			for(std::map< protein::ResidueID, std::vector< std::pair<std::size_t, std::size_t> > >::const_iterator jt=selectable_residue_ids.begin();jt!=selectable_residue_ids.end();++jt)
			{
				const protein::ResidueID& rid=jt->first;
				if(jt!=selectable_residue_ids.begin())
				{
					std::cout << " or ";
				}
				std::cout << "resi " << rid.residue_number << " and chain " << rid.chain_id;
			}
			std::cout << "')\n\n";

			for(std::map< protein::ResidueID, std::vector< std::pair<std::size_t, std::size_t> > >::const_iterator jt=selectable_residue_ids.begin();jt!=selectable_residue_ids.end();++jt)
			{
				const std::vector< std::pair<std::size_t, std::size_t> >& atoms_ids_pairs=jt->second;
				if(color_pymol_selection_at_atomic_level)
				{
					for(std::size_t i=0;i<atoms_ids_pairs.size();i++)
					{
						const std::pair<std::size_t, std::size_t>& atoms_ids_pair=atoms_ids_pairs[i];
						const protein::Atom& a=atoms[atoms_ids_pair.first];
						const protein::Atom& b=atoms[atoms_ids_pair.second];
						std::cout << "cmd.color('" << selection_colorizer->color_string(a, b) << "', 'resi " << a.residue_number << " and name " << atom_name_without_single_quote(a.atom_name) << " and chain " << a.chain_id << "')\n";
					}
				}
				else if(!atoms_ids_pairs.empty())
				{
					const std::pair<std::size_t, std::size_t>& atoms_ids_pair=atoms_ids_pairs.front();
					const protein::Atom& a=atoms[atoms_ids_pair.first];
					const protein::Atom& b=atoms[atoms_ids_pair.second];
					std::cout << "cmd.color('" << selection_colorizer->color_string(a, b) << "', 'resi " << a.residue_number << " and chain " << a.chain_id << "')\n";
				}
			}
			std::cout << "\n";

			std::cout << "cmd.show('sticks', '" << selection_name << "')\n\n";
		}

		std::cout << "cmd.deselect()\n\n";
	}
}
