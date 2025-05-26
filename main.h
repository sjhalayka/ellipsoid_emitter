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

#include <iomanip>

#include <vector>

#include <string>

#include <sstream>

#include <fstream>

#include <set>

#include <map>

#include <utility>

#include <mutex>

#include <thread>

#include <utility>

#include <chrono>

#include <cmath>

using namespace std;


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
vector<line_segment_3> threeD_line_segments_intersected2;




custom_math::vector_3 background_colour(0.0f, 0.0f, 0.0f);
custom_math::vector_3 control_list_colour(1.0f, 1.0f, 1.0f);

bool draw_axis = true;
bool draw_control_list = true;


#ifdef USE_OPENGL
uv_camera main_camera;


GLint win_id = 0;
GLint win_x = 800, win_y = 600;
MyBig camera_w = 100;// 0.01;

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




#endif
