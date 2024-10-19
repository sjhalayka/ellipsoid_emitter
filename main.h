#ifndef main_H
#define main_H

#include "uv_camera.h"
#include "custom_math.h"

using custom_math::vector_3;
using custom_math::vector_4;
using custom_math::line_segment_3;


#include <cstdlib>
#include <GL/glut.h>       //GLUT Library

#include <iostream>
using std::cout;
using std::endl;

#include <iomanip>
using std::setprecision;

#include <vector>
using std::vector;

#include <string>
using std::string;

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


void idle_func(void);
void init_opengl(const int &width, const int &height);
void reshape_func(int width, int height);
void display_func(void);
void keyboard_func(unsigned char key, int x, int y);
void mouse_func(int button, int state, int x, int y);
void motion_func(int x, int y);
void passive_motion_func(int x, int y);

void render_string(int x, const int y, void *font, const string &text);
void draw_objects(void);




size_t n = 100000;

const double receiver_pos = 10;// 10 kiloparsecs
const double receiver_radius = 1.0;// receiver_pos / 10.0;

    


const double G = 6.67430e-11;


vector<vector_3> threeD_oscillators;
vector<vector_3> normals;

const double pi = 4.0f * atanf(1.0f);
const double c_meters = 299792458;


vector<line_segment_3> threeD_line_segments;
vector<line_segment_3> threeD_line_segments_intersected;



custom_math::vector_3 background_colour(0.0f, 0.0f, 0.0f);
custom_math::vector_3 control_list_colour(1.0f, 1.0f, 1.0f);

bool draw_axis = true;
bool draw_control_list = true;

uv_camera main_camera;


GLint win_id = 0;
GLint win_x = 800, win_y = 600;
double camera_w = 10.0;

double camera_fov = 45.0f;
double camera_x_transform = 0;
double camera_y_transform = 0;
double u_spacer = 0.01f;
double v_spacer = 0.5f*u_spacer;
double w_spacer = camera_w*0.01f;
double camera_near = 1.0f;
double camera_far = 10000.0;

bool lmb_down = false;
bool mmb_down = false;
bool rmb_down = false;
int mouse_x = 0;
int mouse_y = 0;






// https://paulbourke.net/geometry/circlesphere/raysphere.c
int RaySphere(vector_3 p1, vector_3 p2, vector_3 sc, double r, double* mu1, double* mu2)
{
    double a, b, c;
    double bb4ac;
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
        return(FALSE);
    }

    *mu1 = (-b + sqrt(bb4ac)) / (2 * a);
    *mu2 = (-b - sqrt(bb4ac)) / (2 * a);

    return(TRUE);
}



// https://www.shadertoy.com/view/7d3BRl
// https://math.stackexchange.com/questions/2931909/normal-of-a-point-on-the-surface-of-an-ellipsoid

vector_4 RayEllipsoid(vector_3 ro, vector_3 rd, vector_3 r)
{
    vector_3 r2 = r * r;
    double a = rd.dot(rd / r2);
    double b = ro.dot(rd / r2);
    double c = ro.dot(ro / r2);
    double h = b * b - a * (c - 1.0);

    if (h < 0.0)
        return vector_4(-1, 0, 0, 0);

    double t = (-b - sqrt(h)) / a;

    vector_3 pos = ro + rd * t;

    return vector_4(t, pos.x, pos.y, pos.z);
}

vector_3 EllipsoidNormal(vector_3 pos, vector_3 ra)
{
    vector_3 normal = (pos / (ra * ra));
    normal.normalize();
    
    return -normal;
}








#endif
