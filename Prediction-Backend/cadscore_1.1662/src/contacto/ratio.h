#ifndef CONTACTO_RATIO_H_
#define CONTACTO_RATIO_H_

namespace contacto
{

struct Ratio
{
	double difference;
	double reference;

	Ratio() : difference(0), reference(0)
	{
	}

	friend std::ostream& operator<<(std::ostream &output, const Ratio &ratio)
	{
		output << ratio.difference << " ";
		output << ratio.reference;
		return output;
	}

	friend std::istream& operator>>(std::istream &input, Ratio &ratio)
	{
		input >> ratio.difference;
		input >> ratio.reference;
		return input;
	}
};

}

#endif /* CONTACTO_RATIO_H_ */
