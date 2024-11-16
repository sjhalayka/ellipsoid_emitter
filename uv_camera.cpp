// Shawn Halayka
// This code and data is in the public domain.


#include "uv_camera.h"

const double pi = 4.0f*atan(1.0f);
const double pi_half = 0.5f*pi;
const double pi_2 = 2.0f*pi;
const double epsilon = 1e-6f;



uv_camera::uv_camera(void)
{
	u = v = 0;
	w = 4;
	fov = 45;
	near_plane = 1.0;
	far_plane = 4.0;
	win_x = win_y = 0;
}

void uv_camera::Set(const MyBig u_rad, const MyBig v_rad, const MyBig w_metres, const MyBig fov_deg, const int width_px, const int height_px, MyBig src_near, MyBig src_far)
{
	u = u_rad;
	v = v_rad;
	w = w_metres;
	near_plane = src_near;
	far_plane = src_far;

	static const double lock = epsilon * 1000.0f;

	if(u < -pi_half + lock)
		u = -pi_half + lock;
	else if(u > pi_half - lock)
		u = pi_half - lock;

	while(v < 0)
		v += pi_2;

	while(v > pi_2)
		v -= pi_2;

	if(w < 0)
		w = 0;
//	else if(w > 10000)
//		w = 10000;

	fov = fov_deg;
	win_x = width_px;
	win_y = height_px;

	Transform();
}


void uv_camera::Transform(void)
{
	Reset();
	Rotate();
	Translate();

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	gluPerspective(
		fov,
		static_cast<GLdouble>(win_x)/static_cast<GLdouble>(win_y), 
		near_plane, far_plane);

	gluLookAt(
		eye.x, eye.y, eye.z, // Eye position.
		eye.x + look_at.x, eye.y + look_at.y, eye.z + look_at.z, // Look at position (not direction).
		up.x, up.y, up.z); // Up direction vector.
}

void uv_camera::Set(void)
{
	// Force a recalculation of the camera vectors and frustum.
	Set(u, v, w, fov, win_x, win_y, near_plane, far_plane);
}
//
//void uv_camera::Set_Large_Screenshot(size_t num_cams, size_t cam_index_x, size_t cam_index_y)
//{
//	// No guarantees about the behaviour of this functionality. It wasn't tested a lot.
//
//	glMatrixMode(GL_PROJECTION);
//	glLoadIdentity();
//
//	// Image plane reference:
//	// http://www.songho.ca/opengl/gl_transform.html
//    const double deg_to_rad = (1.0f/360.0f)*2.0f*pi;
//	double aspect = static_cast<double>(win_x) / static_cast<double>(win_y);
//    double tangent = tan((fov/2.0f)*deg_to_rad);
//    double height = near_plane * tangent; // Half height of near_plane plane.
//    double width = height * aspect; // Half width of near_plane plane.
//
//	double cam_width = 2.0f*width/num_cams;
//	double cam_height = 2.0f*height/num_cams;
//
//	double left = -width + cam_index_x*cam_width;
//	double right = -width + (cam_index_x + 1)*cam_width;
//	double bottom = -height + cam_index_y*cam_height;
//	double top = -height + (cam_index_y + 1)*cam_height;
//
//	// Instead of gluPerspective...
//    glFrustum(left, right, bottom, top, near_plane.convert_to<double>(), far_plane.convert_to<double>());
//
//	gluLookAt(
//		eye.x.convert_to<double>(), eye.y.convert_to<double>(), eye.z.convert_to<double>(), // Eye position.
//		eye.x.convert_to<double>() + look_at.x, eye.y.convert_to<double>() + look_at.y, eye.z.convert_to<double>() + look_at.z, // Look at position (not direction).
//		up.x.convert_to<double>(), up.y.convert_to<double>(), up.z.convert_to<double>()); // Up direction vector.
//}

void uv_camera::Reset(void)
{
	eye.zero();
	look_at.zero();
	up.zero();
	right.zero();

	look_at.z = -1;
	up.y = 1;
	right.x = 1;
}

void uv_camera::Rotate(void)
{
	// Rotate about the world x axis
	look_at.rotate_x(u);
	up.rotate_x(u);
	// Right only rotates on the x axis

	// Rotate about the world y axis
	look_at.rotate_y(v);
	up.rotate_y(v);
	right.rotate_y(v);
}

void uv_camera::Translate(void)
{
	// Place the eye directly across the sphere from the look-at vector's "tip",
	// Then scale the sphere radius by w
	eye.x = -look_at.x*w;
	eye.y = -look_at.y*w;
	eye.z = -look_at.z*w;
}
