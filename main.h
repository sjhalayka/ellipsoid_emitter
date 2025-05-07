#ifndef main_H
#define main_H


//#define USE_OPENGL





#include "uv_camera.h"
#include "custom_math.h"

using custom_math::vector_3;
using custom_math::vector_4;
using custom_math::line_segment_3;



#ifdef USE_OPENGL

#include <cstdlib>
#include <GL/glut.h>       //GLUT Library

#endif


#include <iostream>
using std::cout;
using std::endl;

#include <iomanip>
using std::setprecision;

#include <vector>
using std::vector;

#include <string>
using std::string;
using std::to_string;

#include <sstream>
using std::ostringstream;
using std::istringstream;

#include <fstream>
using std::ofstream;
using std::ifstream;

#include <set>
using std::set;

#include <map>
using std::map;

#include <utility>
using std::pair;

#include <mutex>
using std::mutex;

#include <thread>
using std::thread;

#include <utility>
using std::swap;



#ifdef USE_OPENGL
void idle_func(void);
void init_opengl(const int& width, const int& height);
void reshape_func(int width, int height);
void display_func(void);
void keyboard_func(unsigned char key, int x, int y);
void mouse_func(int button, int state, int x, int y);
void motion_func(int x, int y);
void passive_motion_func(int x, int y);
void render_string(int x, const int y, void* font, const string& text);
void draw_objects(void);
#endif


const MyBig G = 6.67430e-11;
const MyBig c = 299792458;
const MyBig c2 = c * c;
const MyBig c3 = c * c * c;
const MyBig c4 = c * c * c * c;

const MyBig pi = 4.0 * atan(1.0);
const MyBig h = 6.62607015e-34; // Planck’s constant
const MyBig hbar = h / (2.0 * pi);
const MyBig k = 1.380649e-23;

vector<vector_3> threeD_oscillators;
vector<vector_3> normals;
vector<line_segment_3> threeD_line_segments;
vector<line_segment_3> threeD_line_segments_intersected;





custom_math::vector_3 background_colour(0.0f, 0.0f, 0.0f);
custom_math::vector_3 control_list_colour(1.0f, 1.0f, 1.0f);

bool draw_axis = true;
bool draw_control_list = true;


#ifdef USE_OPENGL
uv_camera main_camera;


GLint win_id = 0;
GLint win_x = 800, win_y = 600;
MyBig camera_w = 10;// 0.01;

MyBig camera_fov = 45.0f;
MyBig camera_x_transform = 0;
MyBig camera_y_transform = 0;
MyBig u_spacer = 0.01f;
MyBig v_spacer = 0.5f * u_spacer;
MyBig w_spacer = camera_w * 0.01f;
MyBig camera_near = 0.01;
MyBig camera_far = 10000000.0;

bool lmb_down = false;
bool mmb_down = false;
bool rmb_down = false;
int mouse_x = 0;
int mouse_y = 0;

#endif


MyBig hit_sphere(vector_3 ray_origin, vector_3 ray_direction, const vector_3& center, MyBig radius) {
	
	const vector_3 oc = ray_origin - center;
	const MyBig a = ray_direction.self_dot();
	const MyBig b = 2.0 * oc.dot(ray_direction);
	const MyBig c = oc.self_dot()  - radius * radius;
	const MyBig discriminant = b * b - 4.0 * a * c;
	if (discriminant < 0) {
		return -1.0;
	}
	else {
		return (-b - sqrt(discriminant)) / (2.0 * a);
	}
}


bool intersect(vector_3 ray_origin, vector_3 ray_direction, vector_3 sc, MyBig r, MyBig* t1, MyBig* t2)
{
	//solve for tc
	vector_3 L = sc - ray_origin;
	MyBig tc = L.dot(ray_direction);

	if (tc < 0.0) return false;

	MyBig L2 = L.self_dot();

	MyBig d2 = (tc * tc) - L2;

	MyBig radius2 = r * r;
	if (d2 > radius2) return false;

	//solve for t1c
	MyBig t1c = sqrt(radius2 - d2);

	//solve for intersection points
	*t1 = tc - t1c;
	*t2 = tc + t1c;

	return true;
}
// https://paulbourke.net/geometry/circlesphere/raysphere.c
int RaySphere(vector_3 p1, vector_3 p2, vector_3 sc, MyBig r, MyBig* mu1, MyBig* mu2)
{
	MyBig a, b, c;
	MyBig bb4ac;
	vector_3 dp;

	dp.x = p2.x - p1.x;
	dp.y = p2.y - p1.y;
	dp.z = p2.z - p1.z;
	a = dp.x * dp.x + dp.y * dp.y + dp.z * dp.z;
	b = 2 * (dp.x * (p1.x - sc.x) + dp.y * (p1.y - sc.y) + dp.z * (p1.z - sc.z));
	c = sc.x * sc.x + sc.y * sc.y + sc.z * sc.z;
	c += p1.x * p1.x + p1.y * p1.y + p1.z * p1.z;
	c -= 2 * (sc.x * p1.x + sc.y * p1.y + sc.z * p1.z);
	c -= r * r;
	bb4ac = b * b - 4 * a * c;
	if (fabs(a) < 1e-5 || bb4ac < 0) {
		*mu1 = 0;
		*mu2 = 0;
		return(false);
	}

	*mu1 = (-b + sqrt(bb4ac)) / (2 * a);
	*mu2 = (-b - sqrt(bb4ac)) / (2 * a);

	return(true);
}



// https://www.shadertoy.com/view/7d3BRl
// https://math.stackexchange.com/questions/2931909/normal-of-a-point-on-the-surface-of-an-ellipsoid

vector_4 RayEllipsoid(vector_3 ro, vector_3 rd, vector_3 r)
{
	vector_3 r2 = r * r;
	MyBig a = rd.dot(rd / r2);
	MyBig b = ro.dot(rd / r2);
	MyBig c = ro.dot(ro / r2);
	MyBig h = b * b - a * (c - 1.0);

	if (h < 0.0)
		return vector_4(-1, 0, 0, 0);

	MyBig t = (-b - sqrt(h)) / a;

	vector_3 pos = ro + rd * t;

	return vector_4(t, pos.x, pos.y, pos.z);
}

vector_3 EllipsoidNormal(vector_3 pos, vector_3 ra)
{
	vector_3 normal = (pos / (ra * ra));
	normal.normalize();

	return -normal;
}

//
//void thread_func(
//	size_t& count,
//	const size_t first_index,
//	const size_t last_index,
//	vector_3 sphere_location)
//{
//
//	count = 0;
//
//	for (size_t i = first_index; i <= last_index /*threeD_line_segments.size()*/; i++)
//	{
//		const vector_3 dir = (threeD_line_segments[i].end - threeD_line_segments[i].start).normalize();
//
//		if (1)//dir.dot(sphere_location) > 0)
//		{
//			MyBig mu1 = 0, mu2 = 0;
//
//			if (RaySphere(threeD_line_segments[i].start, threeD_line_segments[i].end, sphere_location, 1.0, &mu1, &mu2))
//			{
//				count++;
//
//				//if (skip_saving_intersected_segments)
//				//    continue;
//
//				line_segment_3 ls_;
//				ls_.start = threeD_line_segments[i].start;
//				ls_.end = threeD_line_segments[i].start + threeD_line_segments[i].end * mu2;
//
//				threeD_line_segments_intersected.push_back(ls_);
//			}
//		}
//	}
//}
//
//size_t get_intersecting_line_count(const vector_3 sphere_location,
//	const MyBig sphere_radius,
//	const MyBig dimension)
//{
//	vector<thread> threads;
//	int num_cpu_threads = 1;// std::thread::hardware_concurrency();
//
//	vector<size_t> counts(num_cpu_threads, 0);
//
//	const size_t div = threeD_line_segments.size() / num_cpu_threads;
//	const size_t modulus = threeD_line_segments.size() % num_cpu_threads;
//
//	vector<size_t> first_index(num_cpu_threads, 0);
//	vector<size_t> last_index(num_cpu_threads, 0);
//
//	first_index[0] = 0;
//	last_index[0] = div - 1;
//
//	for (size_t i = 1; i < (num_cpu_threads - 1); i++)
//	{
//		first_index[i] = first_index[i - 1] + div;
//		last_index[i] = last_index[i - 1] + div;
//	}
//
//	first_index[num_cpu_threads - 1] = first_index[num_cpu_threads - 2] + div;
//	last_index[num_cpu_threads - 1]  = threeD_line_segments.size() - 1;
//
//	for (int i = 0; i < num_cpu_threads; i++)
//		threads.push_back(thread(thread_func, std::ref(counts[i]), first_index[i], last_index[i], sphere_location));
//
//	for (int i = 0; i < num_cpu_threads; i++)
//		threads[i].join();
//
//	size_t count = 0;
//
//	for (size_t i = 0; i < counts.size(); i++)
//		count += counts[i];
//
//	return count;
//}



MyBig get_intersecting_line_count(
	const MyBig n,
	const vector_3 sphere_location,
	const MyBig sphere_radius)
{
	const MyBig big_area =
		4 * pi
		* sphere_location.x * sphere_location.x;

	const MyBig small_area =
		pi
		* sphere_radius * sphere_radius;

	const MyBig ratio =
		small_area
		/ big_area;

	return n * ratio;
}


bool circle_intersect(vector_3 location, const vector_3 normal, const MyBig circle_location, const MyBig circle_radius)
{
	if (normal.dot(circle_location) <= 0) 
		return false;
	
	vector_3 v = normal;

	const MyBig ratio = v.x / circle_location;

	const vector_3 circle_origin(circle_location, 0, 0);

	v.y = v.y / ratio;
	v.z = v.z / ratio;
	v.x = circle_location;

	vector_3 v2;
	v2.x = circle_origin.x - v.x;
	v2.y = circle_origin.y - v.y;
	v2.z = circle_origin.z - v.z;

	if (v2.length() > circle_radius)
		return false;

	return true;
}



size_t get_intersecting_line_count(const vector_3 sphere_location,
	const MyBig sphere_radius,
	const MyBig dimension,
	const bool skip_saving_intersected_segments)
{
	threeD_line_segments_intersected.clear();

	size_t count = 0;

	for (size_t i = 0; i < threeD_oscillators.size(); i++)
	{
		if (circle_intersect(threeD_oscillators[i], normals[i], sphere_location.x, sphere_radius))
		{
			count++;

			if (skip_saving_intersected_segments)
				continue;

			line_segment_3 ls_;
			ls_.start = threeD_oscillators[i];
			ls_.end = threeD_oscillators[i] + normals[i]*3;// +threeD_line_segments[i].end * mu1;

			threeD_line_segments_intersected.push_back(ls_);
		}
	}

	return count;
}




size_t get_intersecting_line_count_sphere(const vector_3 sphere_location,
	const MyBig sphere_radius,
	const MyBig dimension,
	const bool skip_saving_intersected_segments)
{
	threeD_line_segments_intersected.clear();
	size_t count = 0;

	for (size_t i = 0; i < threeD_line_segments.size(); i++)
	{
		const vector_3 dir = (threeD_line_segments[i].end - threeD_line_segments[i].start).normalize();

		if (dir.dot(sphere_location) > 0)
		{
			MyBig mu1 = hit_sphere(threeD_oscillators[i], normals[i], sphere_location, sphere_radius);


			if(mu1 > 0.0)
			//MyBig mu1 = 0, mu2 = 0;

			//if(intersect(threeD_oscillators[i], normals[i], sphere_location, sphere_radius, &mu1, &mu2))

			//if (RaySphere(threeD_line_segments[i].start, threeD_line_segments[i].end, sphere_location, sphere_radius, &mu1, &mu2))
			{
				//if (mu1 < 0 || mu2 < 0)
				//{
				//	//cout << mu1 << "  " << mu2 << endl;
				//	continue;
				//}	

				count++;



				if (skip_saving_intersected_segments)
				    continue;

				line_segment_3 ls_;
				ls_.start = threeD_oscillators[i];
				ls_.end = threeD_oscillators[i] + normals[i]*mu1;// +threeD_line_segments[i].end * mu1;

				threeD_line_segments_intersected.push_back(ls_);
			}
		}
	}

	return count;
}

#include <random>

std::mt19937 generator(0);
std::uniform_real_distribution<long double> dis(0.0, 1.0);


vector_3 RandomUnitVector(void)
{
	MyBig z = dis(generator) * 2.0 - 1.0;
	MyBig a = dis(generator) * 2.0 * pi;

	//MyBig z = static_cast<MyBig>(rand() % RAND_MAX) / static_cast<MyBig>(RAND_MAX) * 2 - 1;
	//MyBig a = static_cast<MyBig>(rand() % RAND_MAX) / static_cast<MyBig>(RAND_MAX) * 2 * pi;
	MyBig r = sqrt(1.0f - z * z);
	MyBig x = r * cos(a);
	MyBig y = r * sin(a);
	return vector_3(x, y, z).normalize();
}

vector_3 slerp(vector_3 s0, vector_3 s1, const MyBig t)
{
	vector_3 s0_norm = s0;
	s0_norm.normalize();

	vector_3 s1_norm = s1;
	s1_norm.normalize();

	const MyBig cos_angle = s0_norm.dot(s1_norm);
	const MyBig angle = acos(cos_angle);

	const MyBig p0_factor = sin((1 - t) * angle) / sin(angle);
	const MyBig p1_factor = sin(t * angle) / sin(angle);

	return s0 * p0_factor + s1 * p1_factor;
}


#endif
