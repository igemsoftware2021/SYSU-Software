#ifndef CONTACTO_INTER_ATOM_CONTACT_H_
#define CONTACTO_INTER_ATOM_CONTACT_H_

#include <iostream>

namespace contacto
{

struct InterAtomContact
{
	int a;
	int b;
	double area;

	InterAtomContact() : a(0), b(0), area(0)
	{
	}

	InterAtomContact(const int a, const int b, const double area) : a(a), b(b), area(area)
	{
	}

	bool operator==(const InterAtomContact& c) const
	{
		return (a==c.a && b==c.b && area==c.area);
	}

	bool operator<(const InterAtomContact& c) const
	{
		if(a<c.a) return true;
		else if(a>c.a) return false;

		if(b<c.b) return true;
		else if(b>c.b) return false;

		if(area<c.area) return true;
		else if(area>c.area) return false;

		return false;
	}

	friend std::ostream& operator<<(std::ostream &output, const InterAtomContact &contact)
	{
		output << contact.a << " ";
		output << contact.b << " ";
		output << contact.area;
		return output;
	}

	friend std::istream& operator>>(std::istream &input, InterAtomContact &contact)
	{
		input >> contact.a;
		input >> contact.b;
		input >> contact.area;
		return input;
	}
};

}

#endif /* CONTACTO_INTER_ATOM_CONTACT_H_ */
