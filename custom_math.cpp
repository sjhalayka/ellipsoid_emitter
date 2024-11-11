
#include "custom_math.h"


custom_math::vector_3::vector_3(const long double &src_x, const long double &src_y, const long double &src_z)
{
	x = src_x;
	y = src_y;
	z = src_z;
}

bool custom_math::vector_3::operator==(const vector_3 &rhs)
{
	if(x == rhs.x && y == rhs.y && z == rhs.z)
		return true;

	return false;
}

bool custom_math::vector_3::operator!=(const vector_3 &rhs)
{
	if(x == rhs.x && y == rhs.y && z == rhs.z)
		return false;

	return true;
}

void custom_math::vector_3::zero(void)
{
	x = y = z = 0;
}

void custom_math::vector_3::rotate_x(const long double &radians)
{
	long double t_y = y;

	y = t_y*cos(radians) + z*sin(radians);
	z = t_y*-sin(radians) + z*cos(radians);
}

void custom_math::vector_3::rotate_y(const long double &radians)
{
	long double t_x = x;

	x = t_x*cos(radians) + z*-sin(radians);
	z = t_x*sin(radians) + z*cos(radians);
}

void custom_math::vector_3::rotate_z(const long double &radians)
{
    long double t_x = x;
    
    x = t_x*cos(radians) + y*-sin(radians);
    y = t_x*sin(radians) + y*cos(radians);
}

custom_math::vector_3 custom_math::vector_3::operator+(const vector_3 &rhs)
{
	return vector_3(x + rhs.x, y + rhs.y, z + rhs.z);
}

custom_math::vector_3 custom_math::vector_3::operator-(const vector_3 &rhs) const
{
	return vector_3(x - rhs.x, y - rhs.y, z - rhs.z);
}

custom_math::vector_3 custom_math::vector_3::operator*(const vector_3 &rhs)
{
	return vector_3(x*rhs.x, y*rhs.y, z*rhs.z);
}

custom_math::vector_3 custom_math::vector_3::operator/(const vector_3& rhs)
{
	return vector_3(x / rhs.x, y / rhs.y, z / rhs.z);
}
custom_math::vector_3 custom_math::vector_3::operator*(const long double &rhs)
{
	return vector_3(x*rhs, y*rhs, z*rhs);
}

custom_math::vector_3 custom_math::vector_3::operator/(const long double &rhs)
{
	return vector_3(x/rhs, y/rhs, z/rhs);
}

custom_math::vector_3 &custom_math::vector_3::operator=(const vector_3 &rhs)
{
	x = rhs.x;
	y = rhs.y;
	z = rhs.z;
	return *this;
}

custom_math::vector_3 &custom_math::vector_3::operator+=(const vector_3 &rhs)
{
	x += rhs.x; y += rhs.y; z += rhs.z;
	return *this;
}

custom_math::vector_3 &custom_math::vector_3::operator*=(const vector_3 &rhs)
{
	x *= rhs.x; y *= rhs.y; z *= rhs.z;
	return *this;
}

custom_math::vector_3 &custom_math::vector_3::operator*=(const long double &rhs)
{
	x *= rhs; y *= rhs; z *= rhs;
	return *this;
}

custom_math::vector_3 custom_math::vector_3::operator-(void)
{
	vector_3 temp;
	temp.x = -x;
	temp.y = -y;
	temp.z = -z;

	return temp;
}

long double custom_math::vector_3::length(void) const
{
	return sqrt(self_dot());
}

custom_math::vector_3 &custom_math::vector_3::normalize(void)
{
	long double len = length();

	if(len != 1)
	{
		x /= len;
		y /= len;
		z /= len;
	}

	return *this;
}

long double custom_math::vector_3::dot(const vector_3 &rhs) const
{
	return x*rhs.x + y*rhs.y + z*rhs.z;
}

long double custom_math::vector_3::self_dot(void) const
{
	return x*x + y*y + z*z;
}

custom_math::vector_3 custom_math::vector_3::cross(const vector_3 &rhs) const
{
	vector_3 cross;
	cross.x = y*rhs.z - rhs.y*z;
	cross.y = z*rhs.x - rhs.z*x;
	cross.z = x*rhs.y - rhs.x*y;

	return cross;
}


custom_math::vector_4::vector_4(const long double &src_x, const long double &src_y, const long double &src_z, const long double &src_w)
{
	x = src_x;
	y = src_y;
	z = src_z;
	w = src_w;
}

void custom_math::vector_4::zero(void)
{
	x = y = z = w = 0;
}

custom_math::vector_4 custom_math::vector_4::operator+(const vector_4 &rhs)
{
	return vector_4(x + rhs.x, y + rhs.y, z + rhs.z, w + rhs.w);
}

custom_math::vector_4 custom_math::vector_4::operator-(const vector_4 &rhs)
{
	return vector_4(x - rhs.x, y - rhs.y, z - rhs.z, w - rhs.w);
}

custom_math::vector_4 custom_math::vector_4::operator*(const vector_4 &rhs)
{
	return vector_4(x*rhs.x, y*rhs.y, z*rhs.z, w*rhs.w);
}

custom_math::vector_4 custom_math::vector_4::operator*(const long double &rhs)
{
	return vector_4(x*rhs, y*rhs, z*rhs, w*rhs);
}

custom_math::vector_4 custom_math::vector_4::operator/(const long double &rhs)
{
	return vector_4(x/rhs, y/rhs, z/rhs, w/rhs);
}

custom_math::vector_4 &custom_math::vector_4::operator=(const vector_4 &rhs)
{
	x = rhs.x;
	y = rhs.y;
	z = rhs.z;
	w = rhs.w;
	return *this;
}

custom_math::vector_4 &custom_math::vector_4::operator+=(const vector_4 &rhs)
{
	x += rhs.x; y += rhs.y; z += rhs.z; w += rhs.w;
	return *this;
}

custom_math::vector_4 &custom_math::vector_4::operator*=(const vector_4 &rhs)
{
	x *= rhs.x; y *= rhs.y; z *= rhs.z; w *= rhs.w;
	return *this;
}

custom_math::vector_4 &custom_math::vector_4::operator*=(const long double &rhs)
{
	x *= rhs; y *= rhs; z *= rhs; w *= rhs;
	return *this;
}

custom_math::vector_4 custom_math::vector_4::operator-(void)
{
	vector_4 temp;
	temp.x = -x;
	temp.y = -y;
	temp.z = -z;
	temp.w = -w;

	return temp;
}

long double custom_math::vector_4::length(void) const
{
	return sqrt(self_dot());
}

custom_math::vector_4 &custom_math::vector_4::normalize(void)
{
	long double len = length();

	if(len != 1)
	{
		x /= len;
		y /= len;
		z /= len;
		w /= len;
	}

	return *this;
}

long double custom_math::vector_4::dot(const vector_4 &rhs) const
{
	return x*rhs.x + y*rhs.y + z*rhs.z + w*rhs.w;
}

long double custom_math::vector_4::self_dot(void) const
{
	return x*x + y*y + z*z + w*w;
}


long double custom_math::d(const long double &a, const long double &b)
{
	return fabs(a - b);
}

long double custom_math::d_3(const vector_3 &a, const vector_3 &b)
{
	return sqrt((a.x - b.x)*(a.x - b.x) + (a.y - b.y)*(a.y - b.y) + (a.z - b.z)*(a.z - b.z));
}

long double custom_math::d_3_sq(const vector_3 &a, const vector_3 &b)
{
	return (a.x - b.x)*(a.x - b.x) + (a.y - b.y)*(a.y - b.y) + (a.z - b.z)*(a.z - b.z);
}

long double custom_math::d_4(const vector_4 &a, const vector_4 &b)
{
	return sqrt((a.x - b.x)*(a.x - b.x) + (a.y - b.y)*(a.y - b.y) + (a.z - b.z)*(a.z - b.z) + (a.w - b.w)*(a.w - b.w));
}

