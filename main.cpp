#include "main.h"


int main(int argc, char** argv)
{
	cout << setprecision(20) << endl;
	srand(0);

	double dimension = 2.5;

	size_t n = 10000000;

	if (dimension <= 2)
		dimension = 2.01;
	else if (dimension > 3)
		dimension = 3;

	double disk_like = 3 - dimension;

	// Start with pseudorandom oscillator locations
	for (size_t i = 0; i < n; i++)
		threeD_oscillators.push_back(RandomUnitVector());

	// Spread the oscillators out, so that they are distributed evenly across
	// the surface of the ellipsoid emitter
	for (size_t i = 0; i < n; i++)
	{
		vector_3 ring;

		ring.x = threeD_oscillators[i].x;
		ring.y = 0;
		ring.z = threeD_oscillators[i].z;

		const double disk_like = 3 - dimension;

		vector_3 s = slerp(threeD_oscillators[i], ring, disk_like);

		threeD_oscillators[i] = s;
	}

	// Get position on oblate ellipsoid
	for (size_t i = 0; i < n; i++)
	{
		vector_3 vec = threeD_oscillators[i];

		const vector_4 rv = RayEllipsoid(vector_3(0, 0, 0), vec, vector_3(1.0, 1.0 - disk_like, 1.0));

		threeD_oscillators[i] = vector_3(rv.y, rv.z, rv.w);
	}

	// Get position and normal on prolate ellipsoid
	for (size_t i = 0; i < n; i++)
	{
		vector_3 vec = threeD_oscillators[i];

		const vector_4 rv = RayEllipsoid(vector_3(0, 0, 0), vec, vector_3(1.0 - disk_like, 1.0, 1.0 - disk_like));

		vector_3 collision_point = vector_3(rv.y, rv.z, rv.w);

		vector_3 normal = EllipsoidNormal(collision_point, vector_3(1.0 - disk_like, 1.0, 1.0 - disk_like));

		normals.push_back(normal);

		line_segment_3 ls;
		ls.start = threeD_oscillators[i];
		ls.end = threeD_oscillators[i] + normals[i] * 1e30;

		threeD_line_segments.push_back(ls);
	}

	// Get intersecting lines
	vector_3 reciever_pos(10, 0, 0);
	size_t collision_count = get_intersecting_line_segments(reciever_pos, 1.0, dimension);

	double ratio = static_cast<double>(collision_count) / static_cast<double>(n);

	cout << ratio << endl;



	glutInit(&argc, argv);
	init_opengl(win_x, win_y);
	glutReshapeFunc(reshape_func);
	glutIdleFunc(idle_func);
	glutDisplayFunc(display_func);
	glutKeyboardFunc(keyboard_func);
	glutMouseFunc(mouse_func);
	glutMotionFunc(motion_func);
	glutPassiveMotionFunc(passive_motion_func);
	//glutIgnoreKeyRepeat(1);
	glutMainLoop();
	glutDestroyWindow(win_id);

	return 0;
}


void idle_func(void)
{
	glutPostRedisplay();
}

void init_opengl(const int& width, const int& height)
{
	win_x = width;
	win_y = height;

	if (win_x < 1)
		win_x = 1;

	if (win_y < 1)
		win_y = 1;

	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
	glutInitWindowPosition(0, 0);
	glutInitWindowSize(win_x, win_y);
	win_id = glutCreateWindow("orbit");

	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LEQUAL);
	glDepthMask(GL_TRUE);
	glShadeModel(GL_SMOOTH);
	glEnable(GL_CULL_FACE);
	glCullFace(GL_BACK);

	glClearColor(static_cast<float>(background_colour.x), static_cast<float>(background_colour.y), static_cast<float>(background_colour.z), 1);
	glClearDepth(1.0f);

	main_camera.Set(0, 0, camera_w, camera_fov, win_x, win_y, camera_near, camera_far);
}

void reshape_func(int width, int height)
{
	win_x = width;
	win_y = height;

	if (win_x < 1)
		win_x = 1;

	if (win_y < 1)
		win_y = 1;

	glutSetWindow(win_id);
	glutReshapeWindow(win_x, win_y);
	glViewport(0, 0, win_x, win_y);

	main_camera.Set(main_camera.u, main_camera.v, main_camera.w, main_camera.fov, win_x, win_y, camera_near, camera_far);
}

// Text drawing code originally from "GLUT Tutorial -- Bitmap Fonts and Orthogonal Projections" by A R Fernandes
void render_string(int x, const int y, void* font, const string& text)
{
	for (size_t i = 0; i < text.length(); i++)
	{
		glRasterPos2i(x, y);
		glutBitmapCharacter(font, text[i]);
		x += glutBitmapWidth(font, text[i]) + 1;
	}
}
// End text drawing code.

void draw_objects(void)
{
	glDisable(GL_LIGHTING);

	glPushMatrix();

	//glScaled(1.0 / receiver_radius, 1.0 / receiver_radius, 1.0 / receiver_radius);
	glPointSize(1.0);
	glLineWidth(1.0f);


	glBegin(GL_POINTS);

	glColor3f(1, 1, 1);

	for (size_t i = 0; i < threeD_oscillators.size(); i++)
		glVertex3d(threeD_oscillators[i].x, threeD_oscillators[i].y, threeD_oscillators[i].z);

	glEnd();




	glEnable(GL_ALPHA);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);


	//glBegin(GL_LINES);

	//glColor4f(0, 1, 0, 0.1f);

	//for (size_t i = 0; i < threeD_oscillators.size(); i++)
	//{
	//	//if (threeD_line_segments[i].start.z > 0 || threeD_line_segments[i].end.z > 0)
	//	//	continue;

	//	if (threeD_oscillators[i].z > 0)
	//		continue;


	//	glVertex3d(threeD_oscillators[i].x, threeD_oscillators[i].y, threeD_oscillators[i].z);
	//	glVertex3d(threeD_oscillators[i].x + normals[i].x, threeD_oscillators[i].y + normals[i].y, threeD_oscillators[i].z + normals[i].z);
	//}

	//glEnd();







	//glBegin(GL_LINES);

	//glColor4f(0, 1, 0, 0.1f);

	//for (size_t i = 0; i < threeD_line_segments.size(); i++)
	//{
	//	if(threeD_line_segments[i].start.z > 0 || threeD_line_segments[i].end.z > 0)
	//		continue;

	//	glVertex3d(threeD_line_segments[i].start.x, threeD_line_segments[i].start.y, threeD_line_segments[i].start.z);
	//	glVertex3d(threeD_line_segments[i].end.x, threeD_line_segments[i].end.y, threeD_line_segments[i].end.z);
	//}

	//glEnd();






	glBegin(GL_LINES);

	glColor4f(0, 0, 1, 0.1f);

	for (size_t i = 0; i < threeD_line_segments_intersected.size(); i++)
	{
		glVertex3d(threeD_line_segments_intersected[i].start.x, threeD_line_segments_intersected[i].start.y, threeD_line_segments_intersected[i].start.z);
		glVertex3d(threeD_line_segments_intersected[i].end.x, threeD_line_segments_intersected[i].end.y, threeD_line_segments_intersected[i].end.z);
	}

	glEnd();






	glDisable(GL_BLEND);



	//// If we do draw the axis at all, make sure not to draw its outline.
	//if(true == draw_axis)
	//{
	//	glBegin(GL_LINES);

	//	glColor3f(1, 0, 0);
	//	glVertex3f(0, 0, 0);
	//	glVertex3f(static_cast<float>(emitter_radius), 0, 0);
	//	glColor3f(0, 1, 0);
	//	glVertex3f(0, 0, 0);
	//	glVertex3f(0, static_cast<float>(emitter_radius), 0);
	//	glColor3f(0, 0, 1);
	//	glVertex3f(0, 0, 0);
	//	glVertex3f(0, 0, static_cast<float>(emitter_radius));

	//	glEnd();
	//}

	glPopMatrix();
}




void display_func(void)
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// Draw the model's components using OpenGL/GLUT primitives.
	draw_objects();

	if (true == draw_control_list)
	{
		// Text drawing code originally from "GLUT Tutorial -- Bitmap Fonts and Orthogonal Projections" by A R Fernandes
		// http://www.lighthouse3d.com/opengl/glut/index.php?bmpfontortho
		glMatrixMode(GL_PROJECTION);
		glPushMatrix();
		glLoadIdentity();
		gluOrtho2D(0, static_cast<float>(win_x), 0, static_cast<float>(win_y));
		glScalef(1, -1, 1); // Neat. :)
		glTranslatef(0, -static_cast<float>(win_y), 0); // Neat. :)
		glMatrixMode(GL_MODELVIEW);
		glPushMatrix();
		glLoadIdentity();

		glColor3d(control_list_colour.x, control_list_colour.y, control_list_colour.z);

		int break_size = 22;
		int start = 20;
		ostringstream oss;

		render_string(10, start, GLUT_BITMAP_HELVETICA_18, string("Mouse controls:"));
		render_string(10, start + 1 * break_size, GLUT_BITMAP_HELVETICA_18, string("  LMB + drag: Rotate camera"));
		render_string(10, start + 2 * break_size, GLUT_BITMAP_HELVETICA_18, string("  RMB + drag: Zoom camera"));

		render_string(10, start + 4 * break_size, GLUT_BITMAP_HELVETICA_18, string("Keyboard controls:"));
		render_string(10, start + 5 * break_size, GLUT_BITMAP_HELVETICA_18, string("  w: Draw axis"));
		render_string(10, start + 6 * break_size, GLUT_BITMAP_HELVETICA_18, string("  e: Draw text"));
		render_string(10, start + 7 * break_size, GLUT_BITMAP_HELVETICA_18, string("  u: Rotate camera +u"));
		render_string(10, start + 8 * break_size, GLUT_BITMAP_HELVETICA_18, string("  i: Rotate camera -u"));
		render_string(10, start + 9 * break_size, GLUT_BITMAP_HELVETICA_18, string("  o: Rotate camera +v"));
		render_string(10, start + 10 * break_size, GLUT_BITMAP_HELVETICA_18, string("  p: Rotate camera -v"));



		custom_math::vector_3 eye = main_camera.eye;
		custom_math::vector_3 eye_norm = eye;
		eye_norm.normalize();

		oss.clear();
		oss.str("");
		oss << "Camera position: " << eye.x << ' ' << eye.y << ' ' << eye.z;
		render_string(10, win_y - 2 * break_size, GLUT_BITMAP_HELVETICA_18, oss.str());

		oss.clear();
		oss.str("");
		oss << "Camera position (normalized): " << eye_norm.x << ' ' << eye_norm.y << ' ' << eye_norm.z;
		render_string(10, win_y - break_size, GLUT_BITMAP_HELVETICA_18, oss.str());

		glPopMatrix();
		glMatrixMode(GL_PROJECTION);
		glPopMatrix();
		// End text drawing code.
	}

	glFlush();
	glutSwapBuffers();
}

void keyboard_func(unsigned char key, int x, int y)
{
	switch (tolower(key))
	{
	case 'w':
	{
		draw_axis = !draw_axis;
		break;
	}
	case 'e':
	{
		draw_control_list = !draw_control_list;
		break;
	}
	case 'u':
	{
		main_camera.u -= u_spacer;
		main_camera.Set();
		break;
	}
	case 'i':
	{
		main_camera.u += u_spacer;
		main_camera.Set();
		break;
	}
	case 'o':
	{
		main_camera.v -= v_spacer;
		main_camera.Set();
		break;
	}
	case 'p':
	{
		main_camera.v += v_spacer;
		main_camera.Set();
		break;
	}

	case ' ':
	{
		//repulse();

		break;
	}




	default:
		break;
	}
}

void mouse_func(int button, int state, int x, int y)
{
	if (GLUT_LEFT_BUTTON == button)
	{
		if (GLUT_DOWN == state)
			lmb_down = true;
		else
			lmb_down = false;
	}
	else if (GLUT_MIDDLE_BUTTON == button)
	{
		if (GLUT_DOWN == state)
			mmb_down = true;
		else
			mmb_down = false;
	}
	else if (GLUT_RIGHT_BUTTON == button)
	{
		if (GLUT_DOWN == state)
			rmb_down = true;
		else
			rmb_down = false;
	}
}

void motion_func(int x, int y)
{
	int prev_mouse_x = mouse_x;
	int prev_mouse_y = mouse_y;

	mouse_x = x;
	mouse_y = y;

	int mouse_delta_x = mouse_x - prev_mouse_x;
	int mouse_delta_y = prev_mouse_y - mouse_y;

	if (true == lmb_down && (0 != mouse_delta_x || 0 != mouse_delta_y))
	{
		main_camera.u -= static_cast<float>(mouse_delta_y) * u_spacer;
		main_camera.v += static_cast<float>(mouse_delta_x) * v_spacer;
	}
	else if (true == rmb_down && (0 != mouse_delta_y))
	{
		main_camera.w -= static_cast<float>(mouse_delta_y) * w_spacer;

		if (main_camera.w < 2.0f)
			main_camera.w = 2.0f;

	}

	main_camera.Set(); // Calculate new camera vectors.
}

void passive_motion_func(int x, int y)
{
	mouse_x = x;
	mouse_y = y;
}




