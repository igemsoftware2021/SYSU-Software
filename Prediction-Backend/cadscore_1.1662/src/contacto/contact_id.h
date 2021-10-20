#ifndef CONTACTO_CONTACT_ID_H_
#define CONTACTO_CONTACT_ID_H_

namespace contacto
{

template<typename SingleIDType>
struct ContactID
{
	typedef SingleIDType SingleID;

	SingleID a;
	SingleID b;

	ContactID()
	{
	}

	ContactID(const SingleID& a, const SingleID& b) : a(a), b(b)
	{
	}

	bool operator==(const ContactID& cid) const
	{
		return (a==cid.a && b==cid.b);
	}

	bool operator<(const ContactID& cid) const
	{
		return (a<cid.a || (a==cid.a && b<cid.b));
	}

	friend std::ostream& operator<<(std::ostream& output, const ContactID& contact_id)
	{
		output << contact_id.a << " " << contact_id.b;
		return output;
	}

	friend std::istream& operator>>(std::istream& input, ContactID& contact_id)
	{
		input >> contact_id.a;
		input >> contact_id.b;
		return input;
	}
};

}

#endif /* CONTACTO_CONTACT_ID_H_ */
