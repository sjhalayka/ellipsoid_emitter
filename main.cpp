#include "main.h"







int main(int argc, char** argv)
{
	cout << setprecision(10);

	const MyBig start_dim = 2.999;
	const MyBig end_dim = 3;

	const size_t dim_res = 2;
	const size_t n = 1000000;

	const size_t output_mod = 10000;

	const MyBig dim_step_size = (end_dim - start_dim) / (dim_res - 1);

	//for (MyBig D = start_dim; D <= end_dim; D += dim_step_size)
	{
		const MyBig D = 3;

		threeD_oscillators.clear();
		normals.clear();
		threeD_line_segments.clear();


		cout << "Allocating memory for oscillators" << endl;
		threeD_oscillators.resize(n);

		cout << "Allocating memory for normals" << endl;
		normals.resize(n);

		cout << "Allocating memory for line segments" << endl;
		threeD_line_segments.resize(n);

		//if (D <= 2)
		//	D = 2.001;
		//else if (D > 3)
		//	D = 3;

		const MyBig disk_like = 3 - D;
		//const MyBig falloff_exponent = D;
		//const MyBig fractionality = 1.0 - 2.0 * (0.5 - fmod(D, 1.0));


		// Start with pseudorandom oscillator locations
		for (size_t i = 0; i < n; i++)
		{
			vector_3 r = RandomUnitVector();
			threeD_oscillators[i] = r;

			if (i % output_mod == 0)
				cout << "Getting pseudorandom locations: " << i << " of " << n << endl;

		}

		//// Spread the oscillators out, so that they are distributed evenly across
		//// the surface of the ellipsoid emitter
		//for (size_t i = 0; i < n; i++)
		//{
		//	vector_3 ring;

		//	ring.x = threeD_oscillators[i].x;
		//	ring.y = 0;
		//	ring.z = threeD_oscillators[i].z;

		//	vector_3 s = slerp(threeD_oscillators[i], ring, disk_like);

		//	threeD_oscillators[i] = s;

		//	if (i % output_mod == 0)
		//		cout << "SLERPing: " << i << " of " << n << endl;
		//}

		//// Get position on oblate ellipsoid
		//for (size_t i = 0; i < n; i++)
		//{
		//	vector_3 vec = threeD_oscillators[i];

		//	const vector_4 rv = RayEllipsoid(vector_3(0, 0, 0), vec, vector_3(1.0, 1.0 - disk_like, 1.0));

		//	threeD_oscillators[i] = vector_3(rv.y, rv.z, rv.w);

		//	if (i % output_mod == 0)
		//		cout << "Getting ellipsoid locations: " << i << " of " << n << endl;
		//}

			



		// Get position and normal on prolate ellipsoid
		for (size_t i = 0; i < n; i++)
		{
			//const vector_4 rv = RayEllipsoid(vector_3(0, 0, 0), threeD_oscillators[i], vector_3(1.0 - disk_like, 1.0, 1.0 - disk_like));

			//vector_3 collision_point = vector_3(rv.y, rv.z, rv.w);

			//vector_3 normal = EllipsoidNormal(collision_point, vector_3(1.0 - disk_like, 1.0, 1.0 - disk_like));

			vector_3 normal = RandomUnitVector();
			if (threeD_oscillators[i].dot(normal) < 0)
				normal = -normal;

			//vector_3 normal = threeD_oscillators[i];
			//normal.normalize();

			normals[i] = normal;

			line_segment_3 ls;
			ls.start = threeD_oscillators[i];
			ls.end = threeD_oscillators[i] + normals[i];// *1e30;

			threeD_line_segments[i] = ls;

			if (i % output_mod == 0)
				cout << "Getting elipsoid normals: " << i << " of " << n << endl;
		}


		string filename = to_string(D) + ".txt";

		ofstream out_file(filename.c_str());
		out_file << setprecision(30);

		const MyBig receiver_radius = 1.0;

		const MyBig start_distance = (1 + receiver_radius);
		const MyBig end_distance = 100.0;

		const size_t distance_res = 100;

		const MyBig distance_step_size = (end_distance - start_distance) / (distance_res - 1);

		for (MyBig r = start_distance; r <= end_distance; r += distance_step_size)
		{
			//const vector_3 receiver_pos(r, 0, 0);

			//const MyBig epsilon = 1.0;

			//vector_3 receiver_pos_plus = receiver_pos;
			//receiver_pos_plus.x += epsilon;

			//const long long signed int collision_count_plus = get_intersecting_line_count(receiver_pos_plus, receiver_radius, D, true);

			//const long long signed int collision_count = get_intersecting_line_count(receiver_pos, receiver_radius, D, false);

			//vector_3 gradient;
			//gradient.x = static_cast<MyBig>(collision_count_plus - collision_count) / epsilon;

			//const MyBig gradient_strength =
			//	-gradient.x
			//	/ (receiver_radius * receiver_radius);

			//const MyBig gradient_strength_ =
			//	n / (2 * pow(receiver_pos.x, 3.0));

			const MyBig emitter_mass = c2 * 1.0 / (2.0 * G);

			//const MyBig time_dilation = sqrt(1.0 - 1.0 / receiver_pos.x);

			//const MyBig newton_strength = gradient_strength * receiver_pos.x * c * hbar * log(2.0) / (k * 2.0 * pi * emitter_mass);




			//cout << gradient_strength / gradient_strength_ << endl;

			//const MyBig newton_strength_ =
			//	n *
			//	c *
			//	hbar *
			//	log(2.0)
			//	/
			//	(k *
			//		pow(receiver_pos.x, 2.0)
			//		*
			//		emitter_mass *
			//		4.0 *
			//		pi);

			//const MyBig newton_strength__ = G *  emitter_mass / ( pow(receiver_pos.x, 3.0));




			//cout << newton_strength__ / (newton_strength) << endl;


			////	const MyBig newton_strength =
			////		G * emitter_mass / (/*time_dilation * */ pow(receiver_pos.x, 2.0));

			//	//cout << gradient_strength / gradient_strength_ << endl;






			const vector_3 receiver_pos(r, 0, 0);

			const MyBig epsilon = 1.0;

			vector_3 receiver_pos_plus = receiver_pos;
			receiver_pos_plus.x += epsilon;



			//const long long signed int collision_count_plus =
			//	get_intersecting_line_count(
			//		receiver_pos_plus, receiver_radius, 3.0, true);

			//const long long signed int collision_count =
			//	get_intersecting_line_count(
			//		receiver_pos, receiver_radius, 3.0, false);

			const MyBig collision_count_plus = 
				get_intersecting_line_count(n, 
					receiver_pos_plus, receiver_radius);

			const MyBig collision_count =
				get_intersecting_line_count(n,
					receiver_pos, receiver_radius);

			const MyBig gradient =
				(MyBig(collision_count_plus) - MyBig(collision_count))
				/ epsilon;

			const MyBig gradient_strength =
				-gradient
				/ (receiver_radius * receiver_radius);


			//const MyBig newton_strength_ =
			//	 G * k * emitter_mass / (4 * c2 * pow(receiver_pos.x, 2.0));

			const MyBig newton_strength_ =
				G  * emitter_mass  / (pow(receiver_pos.x, 2.0));

			const MyBig newton_strength__ =
				gradient_strength * receiver_pos.x * c * hbar * log(2)
				/ (k * 2 * pi * emitter_mass);


			//const MyBig newton_strength =
			//	n *
			//	c *
			//	hbar *
			//	log(2.0)
			//	/
			//	(k *
			//		pow(receiver_pos.x, 2.0)
			//		*
			//		emitter_mass *
			//		4.0 *
			//		pi);

			cout << newton_strength__ / newton_strength_ << endl;

			cout << "r: " << r << " gradient strength: "
				<< gradient_strength << endl;

			out_file << r << " " << gradient_strength << endl;








			//cout << "D: " << D << " r: " << r << " " << gradient_strength << endl;

			//out_file << r << " " << gradient_strength << endl;
		}

		out_file.close();
	}

	//return 0;






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


	glBegin(GL_LINES);

	glColor4f(1.0f, 0.5f, 0, 0.2f);

	for (size_t i = 0; i < threeD_oscillators.size(); i++)
	{
		//if (threeD_line_segments[i].start.z > 0 || threeD_line_segments[i].end.z > 0)
		//	continue;

		//if (threeD_oscillators[i].z > 0)
		//	continue;


		glVertex3d(threeD_oscillators[i].x, threeD_oscillators[i].y, threeD_oscillators[i].z);
		glVertex3d(threeD_oscillators[i].x + normals[i].x, threeD_oscillators[i].y + normals[i].y, threeD_oscillators[i].z + normals[i].z);
	}

	glEnd();


	//glColor4f(1.0f, 1.0f, 1.0f, 1.0f);

	//glutSolidSphere(1.0, 20, 20);





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

	glColor4f(0, 0, 1, 1.0f);

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
