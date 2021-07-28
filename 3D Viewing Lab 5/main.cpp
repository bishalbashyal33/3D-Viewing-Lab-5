#include"graphics.h"
#include<chrono>
#include<string>
#include "mesh.h"

poly3d p1;
poly3d p2;
cam  vcam;
ptlight l1;
float fps = 0.0f;
int angle = 0, state = 0;
void draw_line(int x1, int x2, int y1, int y2) {
	int dx, dy, i, e;
	int incx, incy, inc1, inc2;
	int x, y;

	dx = x2 - x1;
	dy = y2 - y1;

	if (dx < 0) dx = -dx;
	if (dy < 0) dy = -dy;
	incx = 1;
	if (x2 < x1) incx = -1;
	incy = 1;
	if (y2 < y1) incy = -1;
	x = x1; y = y1;
	if (dx > dy) {
		putpixel(x, y, 20);
		e = 2 * dy - dx;
		inc1 = 2 * (dy - dx);
		inc2 = 2 * dy;
		for (i = 0; i < dx; i++) {
			if (e >= 0) {
				y += incy;
				e += inc1;
			}
			else
				e += inc2;
			x += incx;
			putpixel(x, y, 20);
		}

	}
	else {
		putpixel(x, y, 20);
		e = 2 * dx - dy;
		inc1 = 2 * (dx - dy);
		inc2 = 2 * dx;
		for (i = 0; i < dy; i++) {
			if (e >= 0) {
				x += incx;
				e += inc1;
			}
			else
				e += inc2;
			y += incy;
			putpixel(x, y, 20);
		}
	}
}

void myDisplay(int x1, int x2, int y1, int y2) {
	draw_line(x1, x2, y1, y2);
	glFlush();
}

void renderer() {



}

void readKeyboard(unsigned char key, int x, int y) {

	vertex temp = vcam.lookdir();
	vertex forward = temp * 0.1f;

	//cout << "here";
	if (key == 'w') {
		angle++;
		state = 0;
		update(0);
		display();

	}
	if (key == 'a') {


		angle++;
		state = 1;
		update(0);
		display();
	}

	if (key == 's') {
		angle++;
		state = 2;
		update(0);
		display();
	}

	if (key == 'd') {
		angle++;
		state = 3;
		update(0);
		display();

	}

	if (key == 't') {
		vcam.vCamera = vcam.vCamera + forward;
		//vcam.vCamera.y -= 0.1f;
		update(0);
		display();


	}
	if (key == 'f') {

		vcam.yaw -= 0.1f;
		update(0);
		display();

	}
	if (key == 'h') {
		vcam.yaw += 0.1f;
		update(0);
		display();

	}
	if (key == 'g') {
		vcam.vCamera = vcam.vCamera - forward;
		//vcam.vCamera.y += 0.1f;
		update(0);
		display();

	}
	if (key == 'i') {
		l1.pointlight.x += +1.0f;
		update(0);
		display();

	}
	if (key == 'k') {
		l1.pointlight.x -= +1.0f;
		update(0);
		display();

	}
	if (key == 'j') {
		l1.pointlight.y += +1.0f;
		update(0);
		display();

	}
	if (key == 'j') {
		l1.pointlight.y -= +1.0f;
		update(0);
		display();

	}
	if (key == 'o') {
		l1.pointlight.z += +1.0f;
		update(0);
		display();

	}

	if (key == 'p') {
		l1.pointlight.z -= +1.0f;
		update(0);
		display();

	}


	/*vcam.vCamera = vcam.vCamera + temp;*/


}

void readarrows(int key, int x, int y) {

	if (key == GLUT_KEY_UP) {
		//p1.vCam.vCamera.x += 10.0f;
		vcam.vCamera.y -= 0.1f;

		update(0);
		display();
	}
	if (key == GLUT_KEY_DOWN) {


		vcam.vCamera.y += 0.1f;

		update(0);
		display();

	}
	if (key == GLUT_KEY_LEFT) {


		vcam.vCamera.x -= 0.1f;

		update(0);
		display();

	}
	if (key == GLUT_KEY_RIGHT) {


		vcam.vCamera.x += 1.0f;

		update(0);
		display();

	}

}



void mouse(int button, int state, int x, int y) {



}




//int main(int argc, char** argv) {
//
//	glutInit(&argc, argv);
//	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
//	glutInitWindowSize(getmaxx(), getmaxy());
//	glutInitWindowPosition(0, 0);
//	glutCreateWindow("3D Renderer");
//	myInit();
//	glutDisplayFunc(renderer);	
//	glutKeyboardFunc(readKeyboard);
//	glutMouseFunc(mouse);
//	glutMainLoop();
//	return 0;
//}
void update(int j) {

	cleargrid();
	int i = 60;
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	//DDAlgorithm(0, 0, 200, 200, -1.0f);
	auto lastframe = std::chrono::high_resolution_clock::now();
	//p2.render(angle, state, vcam, l1, 6.0f);
	p1.render(angle, state, vcam, l1);
	float deltatime = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - lastframe).count();
	lastframe = std::chrono::high_resolution_clock::now();
	fps = (1e6 / deltatime);
	glutSetWindowTitle(to_string(fps).c_str());
	glFlush();
	i = i % 360;
}


int main(int argc, char** argv) {

	glutInit(&argc, argv);
	glutInitWindowSize(window_width, window_height);
	glutInitWindowPosition(0, 0);
	//const char* str = to_string(fps).c_str();
	glutCreateWindow("Renderer");
	glutReshapeFunc(reshape);
	glutKeyboardFunc(readKeyboard);
	glutSpecialFunc(readarrows);
	glutDisplayFunc(display);
	initcanvas(argc, argv);
	update(0);
	glutMainLoop();
	return 0;
}

