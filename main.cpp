#include "main.h"







bool intersect_AABB(const vector_3 min_location, const vector_3 max_location, const vector_3 ray_origin, const vector_3 ray_dir, MyBig& tmin, MyBig& tmax)
{
	tmin = (min_location.x - ray_origin.x) / ray_dir.x;
	tmax = (max_location.x - ray_origin.x) / ray_dir.x;

	if (tmin > tmax) swap(tmin, tmax);

	MyBig tymin = (min_location.y - ray_origin.y) / ray_dir.y;
	MyBig tymax = (max_location.y - ray_origin.y) / ray_dir.y;

	if (tymin > tymax) swap(tymin, tymax);

	if ((tmin > tymax) || (tymin > tmax))
		return false;

	if (tymin > tmin)
		tmin = tymin;

	if (tymax < tmax)
		tmax = tymax;

	MyBig tzmin = (min_location.z - ray_origin.z) / ray_dir.z;
	MyBig tzmax = (max_location.z - ray_origin.z) / ray_dir.z;

	if (tzmin > tzmax) swap(tzmin, tzmax);

	if ((tmin > tzmax) || (tzmin > tmax))
		return false;

	if (tzmin > tmin)
		tmin = tzmin;

	if (tzmax < tmax)
		tmax = tzmax;

	return true;
}


// beta is density
// alpha is gradient
void get_density_and_gradient(MyBig& beta, MyBig& alpha)
{
	const size_t n = 100000000;

	const MyBig emitter_radius = sqrt((n * G * hbar * log(2.0)) / (k * c3 * pi));
	const MyBig emitter_area = 4 * pi * emitter_radius * emitter_radius;
	const MyBig mass = c2 * emitter_radius / (2.0 * G);
	const MyBig mass2 = sqrt((n * c * hbar * log(2.0)) / (4 * G * k * pi));



	cout << "emitter_radius: " << emitter_radius << endl;
	cout << "mass: " << mass << endl;
	//cout << "mass2 : " << mass2 << endl;
	cout << "n: " << n << endl;
	cout << endl;

	const MyBig start_dim = 2.01; // Minimum 2.000001
	const MyBig end_dim = 3.0; // Maximum 3
	const size_t dim_res = 2; // Larger than 1
	const MyBig dim_step_size = (end_dim - start_dim) / (dim_res - 1);

	//for (MyBig D = start_dim; D <= end_dim; D += dim_step_size)
	{
		threeD_oscillators.clear();
		normals.clear();
		threeD_line_segments.clear();
		threeD_line_segments_intersected.clear();
		threeD_line_segments_intersected2.clear();
		//threeD_oscillators.resize(n);
		normals.resize(n);
		//threeD_line_segments.resize(n);	

		const MyBig D = 2.001;

		const MyBig disk_like = 3 - D;

		//const MyBig fractionality = 1.0 - 2.0 * (0.5 - fmod(D, 1.0));

		// Get normal on prolate ellipsoid
		for (size_t i = 0; i < n; i++)
		{
			vector_3 oscillator = RandomUnitVector() * 0.01; // Something much smaller than unit vectors

			const vector_4 rv = RayEllipsoid(vector_3(0, 0, 0), oscillator, vector_3(1.0 - disk_like, 1.0, 1.0 - disk_like));

			normals[i] = EllipsoidNormal(vector_3(rv.y, rv.z, rv.w), vector_3(1.0 - disk_like, 1.0, 1.0 - disk_like));
			//threeD_oscillators[i] = vector_3(0, 0, 0);

			//line_segment_3 ls;
			//ls.start = threeD_oscillators[i];
			//ls.end = threeD_oscillators[i] + normals[i];

			//threeD_line_segments[i] = ls;
		}


		const MyBig start_distance = 20.0;
		const MyBig end_distance = 100.0;
		const size_t distance_res = 10;

		const MyBig distance_step_size =
			(end_distance - start_distance)
			/ (distance_res - 1);

		//for (size_t step_index = 0; step_index < distance_res; step_index++)
		{
			//const MyBig r =
			//	start_distance + step_index * distance_step_size;

			const MyBig r =
				start_distance;




			MyBig epsilon = 0.000001;

			MyBig count0 = 0;
			MyBig density0 = 0;

			// Unit box
			vector_3 min_location(-0.5 + r, -0.5, -0.5);
			vector_3 max_location(0.5 + r, 0.5, 0.5);

			for (size_t i = 0; i < n; i++)
			{
				vector_3 ray_origin = vector_3(0, 0, 0);//threeD_oscillators[i];
				vector_3 ray_dir = normals[i];

				MyBig tmin = 0, tmax = 0;

				if (intersect_AABB(min_location, max_location, ray_origin, ray_dir, tmin, tmax))
				{
					// If pointing in the wrong direction
					if (tmin < 0 || tmax < 0)
						continue;

					vector_3 ray_hit_start = ray_origin + ray_dir * tmin;
					vector_3 ray_hit_end = ray_origin + ray_dir * tmax;

					MyBig l = (ray_hit_end - ray_hit_start).length();

					line_segment_3 ls;
					ls.start = ray_hit_start;
					ls.end = ray_hit_end;

					//threeD_line_segments_intersected.push_back(ls);
					count0 += 1;
					density0 += l;
				}
			}

			MyBig count1 = 0;
			MyBig density1 = 0;

			// Unit box
			min_location = vector_3(-0.5 + r + epsilon, -0.5, -0.5);
			max_location = vector_3(0.5 + r + epsilon, 0.5, 0.5);

			for (size_t i = 0; i < n; i++)
			{
				vector_3 ray_origin = vector_3(0, 0, 0);//threeD_oscillators[i];
				vector_3 ray_dir = normals[i];

				MyBig tmin = 0, tmax = 0;

				if (intersect_AABB(min_location, max_location, ray_origin, ray_dir, tmin, tmax))
				{
					// If pointing in the wrong direction
					if (tmin < 0 || tmax < 0)
						continue;

					vector_3 ray_hit_start = ray_origin + ray_dir * tmin;
					vector_3 ray_hit_end = ray_origin + ray_dir * tmax;

					MyBig l = (ray_hit_end - ray_hit_start).length();

					line_segment_3 ls;
					ls.start = ray_hit_start;
					ls.end = ray_hit_end;

					//threeD_line_segments_intersected2.push_back(ls);
					count1 += 1;
					density1 += l;
				}
			}

			density0 /= (max_location.x - min_location.x) * (max_location.y - min_location.y) * (max_location.z - min_location.z);
			density1 /= (max_location.x - min_location.x) * (max_location.y - min_location.y) * (max_location.z - min_location.z);

			alpha = (density1 - density0) / epsilon;
			beta = density0;

			MyBig g = -alpha * pi;
			MyBig g_ = n / (2.0 * pow(r, D));

			cout << g << " " << g_ << endl;

			cout << "g " <<  g / g_ << endl;

		//	MyBig g_N = g * r * c * hbar * log(2.0) / (k * pi * 2.0 * mass);
			

			MyBig g_N = G * mass / (r * r);

			MyBig g_N_flat = sqrt((g * G * c * hbar * log(2.0)) / (2 * k * pi * pow(r, 1.0 - disk_like)));

			cout << g_N_flat / g_N << endl;
			


			////MyBig g_N_3 = sqrt((G * g * c * hbar * log(2.0)) / (2*r * k * pi));
			//MyBig g_N_3 = sqrt((g * G * c * hbar * log(2.0)) / (2 * k * r * r * pi));

			////cout << g_N_ / g_N_3 << endl;

			//MyBig v_3 = pow((g * G * r * c * hbar * log(2.0)) / (2 * k * pi), 1.0f / 4.0f);
			//MyBig v_ = sqrt(G * mass / r);
			////cout << v_3 << endl;


			//MyBig v = sqrt(g * r * r * c * hbar * log(2.0) / (k * pi * 2.0 * mass));

			//cout << v / v_3 << endl;

			//MyBig a_ = v * v / r;
			//MyBig a__ = g * r * c * hbar * log(2.0) / (k * pi * 2.0 * mass);
			//MyBig a = -alpha * r * c * hbar * log(2.0) / (k * 2.0 * mass);


			////cout << a << " " << a_ << endl;

			////cout << g_N / g_N_3 << endl;


		}
	}
}

int main(int argc, char** argv)
{

	std::chrono::high_resolution_clock::time_point global_time_start = std::chrono::high_resolution_clock::now();

	MyBig density = 0, gradient = 0;
	get_density_and_gradient(density, gradient);

	std::chrono::high_resolution_clock::time_point global_time_end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<float, std::milli> elapsed = global_time_end - global_time_start;

	cout << elapsed.count() / 1000.0f << " seconds elapsed" << endl;









#ifdef USE_OPENGL

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

#endif

}




#ifdef USE_OPENGL

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



	glPointSize(2.0);
	glLineWidth(1.0f);


	//glBegin(GL_POINTS);

	//glColor3f(1, 1, 1);

	//for (size_t i = 0; i < threeD_oscillators.size(); i++)
	//{
	//	if (i % 1000 != 0)
	//		continue;

	//	glVertex3d(threeD_oscillators[i].x, threeD_oscillators[i].y, threeD_oscillators[i].z);
	//}

	//glEnd();




	glEnable(GL_ALPHA);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);


	//glBegin(GL_LINES);

	//glColor4f(1.0f, 0.5f, 0, 0.2f);

	//for (size_t i = 0; i < threeD_line_segments.size(); i++)
	//{
	//	//if (threeD_line_segments[i].start.z > 0 || threeD_line_segments[i].end.z > 0)
	//	//	continue;

	//	//if (threeD_oscillators[i].z > 0)
	//	//	continue;


	//	glVertex3d(threeD_line_segments[i].start.x, threeD_line_segments[i].start.y, threeD_line_segments[i].start.z);
	//	glVertex3d(threeD_line_segments[i].end.x, threeD_line_segments[i].end.y, threeD_line_segments[i].end.z);
	//}

	//glEnd();




	//glBegin(GL_LINES);

	//glColor4f(1.0f, 0.5f, 0, 0.2f);

	//for (size_t i = 0; i < threeD_line_segments_intersected.size(); i++)
	//{
	//	//if (threeD_line_segments[i].start.z > 0 || threeD_line_segments[i].end.z > 0)
	//	//	continue;

	//	//if (threeD_oscillators[i].z > 0)
	//	//	continue;


	//	glVertex3d(threeD_line_segments_intersected[i].start.x, threeD_line_segments_intersected[i].start.y, threeD_line_segments_intersected[i].start.z);
	//	glVertex3d(threeD_line_segments_intersected[i].end.x, threeD_line_segments_intersected[i].end.y, threeD_line_segments_intersected[i].end.z);
	//}



//	glEnd();


	//glColor4f(1.0f, 1.0f, 1.0f, 1.0f);

	//glutSolidSphere(1.0, 20, 20);





	glBegin(GL_LINES);

	glColor4f(0, 1, 0, 0.1f);

	for (size_t i = 0; i < threeD_line_segments.size(); i++)
	{
		glVertex3d(threeD_line_segments[i].start.x, threeD_line_segments[i].start.y, threeD_line_segments[i].start.z);
		glVertex3d(threeD_line_segments[i].end.x, threeD_line_segments[i].end.y, threeD_line_segments[i].end.z);
	}

	glEnd();



	glBegin(GL_LINES);

	glColor4f(0, 0.5, 1.0, 0.2f);

	//cout << threeD_line_segments_intersected.size() << endl;

	for (size_t i = 0; i < threeD_line_segments_intersected2.size(); i++)
	{
		glVertex3d(threeD_line_segments_intersected2[i].start.x, threeD_line_segments_intersected2[i].start.y, threeD_line_segments_intersected2[i].start.z);
		glVertex3d(threeD_line_segments_intersected2[i].end.x, threeD_line_segments_intersected2[i].end.y, threeD_line_segments_intersected2[i].end.z);
	}

	glEnd();


	glBegin(GL_LINES);

	glColor4f(1, 0.5, 0.0, 0.2f);

	//cout << threeD_line_segments_intersected.size() << endl;

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

		if (main_camera.w < 0.00001)
			main_camera.w = 0.00001;

	}

	main_camera.Set(); // Calculate new camera vectors.
}

void passive_motion_func(int x, int y)
{
	mouse_x = x;
	mouse_y = y;
}

#endif
