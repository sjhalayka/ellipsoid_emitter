#ifndef uv_camera_h
#define uv_camera_h

#ifdef USE_OPENGL
#include <cstdlib>
#include <GL/glut.h>       //GLUT Library
#endif


#include "custom_math.h"




// UV camera
//
// latitude:     | longitude:    | radius:       |
//       *_*_    |        ___    |        ___    |
//      */   \   |       /   \   |       /   \   |
// u:  *|  x  |  |  v:  |**x**|  |  w:  |  x**|  |
//      *\___/   |       \___/   |       \___/   |
//       * *     |               |               |
// 

class uv_camera
{
public:
	// Use as read-only
	MyBig u, v, w, fov;
	int win_x, win_y;
    custom_math::vector_3 eye, look_at, up, right;
	MyBig near_plane;
	MyBig far_plane;

public:
	uv_camera(void);

	// Must initialize or change camera settings through these two functions
	void Set(const MyBig u_rad, const MyBig v_rad, const MyBig w_metres, const MyBig fov_deg, const int width_px, const int height_px, MyBig src_near, MyBig src_far);
	void Set(void);
	void Set_Large_Screenshot(size_t num_cams, size_t cam_num_x, size_t cam_num_y);
protected:
	void Transform(void);
	void Reset(void);
	void Rotate(void);
	void Translate(void);
};


#endif
