#ifndef CONTACTO_RESIDUE_CONTACT_AREA_DIFFERENCE_BASIC_SCORING_FUNCTORS_H_
#define CONTACTO_RESIDUE_CONTACT_AREA_DIFFERENCE_BASIC_SCORING_FUNCTORS_H_

#include <cmath>

namespace contacto
{

struct RawDifferenceProducer
{
	double operator()(const double target_area, const double model_area) const
	{
		return (target_area-model_area);
	}
};

struct SimpleDifferenceProducer
{
	double operator()(const double target_area, const double model_area) const
	{
		return fabs(target_area-model_area);
	}
};

struct BoundedDifferenceProducer
{
	double operator()(const double target_area, const double model_area) const
	{
		return std::min(fabs(target_area-model_area), target_area);
	}
};

struct SimpleReferenceProducer
{
	double operator()(const double target_area, const double /*model_area*/) const
	{
		return target_area;
	}
};

struct SummingReferenceProducer
{
	double operator()(const double target_area, const double model_area) const
	{
		return (target_area+model_area);
	}
};

}

#endif /* CONTACTO_RESIDUE_CONTACT_AREA_DIFFERENCE_BASIC_SCORING_FUNCTORS_H_ */
