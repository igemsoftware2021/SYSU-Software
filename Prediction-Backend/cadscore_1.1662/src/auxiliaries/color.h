#ifndef AUXILIARIES_COLOR_H_
#define AUXILIARIES_COLOR_H_

#include <algorithm>

namespace auxiliaries
{

struct Color
{
	unsigned char r;
	unsigned char g;
	unsigned char b;

	Color() : r(0), g(0), b(0)
	{
	}

	Color(const unsigned char r, const unsigned char g, const unsigned char b) : r(r), g(g), b(b)
	{
	}

	static Color from_code(const unsigned int rgb)
	{
		return Color(((rgb&0xFF0000) >> 16), ((rgb&0x00FF00) >> 8), (rgb&0x0000FF));
	}

	static Color from_temperature_to_blue_white_red(const double t)
	{
		Color c;
		if(t<0)
		{
			c.b=255;
		}
		else if(t>1)
		{
			c.r=255;
		}
		else if(t<=0.5)
		{
			c.b=255;
			c.r=static_cast<unsigned char>(255*(t/0.5));
			c.g=c.r;
		}
		else if(t>0.5)
		{
			c.r=255;
			c.b=static_cast<unsigned char>(255*(1-(t-0.5)/0.5));
			c.g=c.b;
		}
		return c;
	}

	static Color from_temperature_to_green_yellow_red(const double t)
	{
		Color c;
		if(t<0)
		{
			c.g=255;
		}
		else if(t>1)
		{
			c.r=255;
		}
		else if(t<=0.5)
		{
			c.g=255;
			c.r=static_cast<unsigned char>(255*(t/0.5));
		}
		else if(t>0.5)
		{
			c.r=255;
			c.g=static_cast<unsigned char>(255*(1-(t-0.5)/0.5));
		}
		return c;
	}

	static Color from_two_values_to_green_yellow_red(const double a, const double b)
	{
		Color c;
		const double max_val=std::max(a, b);
		if(max_val>0)
		{
			c.r=static_cast<unsigned char>(a/max_val*255);
			c.g=static_cast<unsigned char>(b/max_val*255);
		}
		return c;
	}

	static Color from_id(const long id)
	{
		const long generator=1234567;
		const long limiter=0xFFFFFF;
		return from_code(static_cast<unsigned int>(((id+1)*generator)%limiter));
	}


	double r_double() const
	{
		return (static_cast<double>(r)/255.0);
	}

	double g_double() const
	{
		return (static_cast<double>(g)/255.0);
	}

	double b_double() const
	{
		return (static_cast<double>(b)/255.0);
	}
};

}

#endif /* AUXILIARIES_COLOR_H_ */
