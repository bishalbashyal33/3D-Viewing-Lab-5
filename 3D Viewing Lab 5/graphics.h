#pragma once
#pragma once
#ifdef __APPLE__
#include <GLUT/glut.h>
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#else
#include "structures.h"
#include <stdio.h>
#include<stdlib.h>
#include<windows.h>
#include <GL.H>
#include<GLU.H>
#include<GLAux.h>
#include<glut.h>
#include <iostream>

#include "maths.h"
#endif

#include <iostream>




#define WHITE 1.0, 1.0, 1.0
#define BLACK 0.0, 0.0, 0.0

// values are read from "game.config"
GLint FPS = 60;
GLint window_width = 800;
GLint window_height = 600;

bool* grid;
vec3* color;
vertex rgb;
float* zBuffer;

float t = 0;
void display() {
    glClear(GL_COLOR_BUFFER_BIT);
    glLoadIdentity();
    //     glPointSize(5.5f);

    auto alpla = sin(((t += 0.1f) + 1) / 2.f);

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glBegin(GL_POINTS);
    for (GLint x = 0; x < window_width; ++x) {
        for (GLint y = 0; y < window_height; ++y) {
            if (grid[x + y * window_width]) {
                vec3& c = color[x + y * window_width];
                // glColor3f(c.x, c.y, c.z);
                glColor4f(static_cast<float>(c.x), static_cast<float>(c.y), static_cast<float>(c.z), 1);

                // glColor4f(rgb.r,rgb.g,rgb.b,-c.z);
                // std::cout << c.z<<std::endl;
                 // glColor3f(static_cast<float>(rand()) / static_cast<float>(RAND_MAX), static_cast<float>(rand()) / static_cast<float>(RAND_MAX), static_cast<float>(rand()) / static_cast<float>(RAND_MAX));
                glVertex2i(x, y);
            }
        }
    }
    glEnd();

    glFlush();
    glutSwapBuffers();
}

void cleargrid() {
    for (GLint x = 0; x < window_width * window_height; ++x) {
        grid[x] = false;
        color[x] = vec3(0);
        zBuffer[x] = -99999999999.9999999999;
    }
}

void getPixel(int x, int y, bool& pixel, vec3& col) {
    if (x < window_width && x >= 0 && y < window_height && y >= 0) {
        col = color[x + y * window_width];
        pixel = grid[x + y * window_width];
    }
}

void reshape(int w, int h) {
    auto oldWidth = window_width;
    auto oldHeight = window_height;
    window_width = w;
    window_height = h;

    bool* newGrid = new bool[window_height * window_width];
    vec3* newcolor = new vec3[window_height * window_width];
    float* newzBuffer = new float[window_height * window_width];

    for (GLint x = 0; x < window_width; ++x) {
        for (GLint y = 0; y < window_height; ++y) {
            if (x < oldWidth && x >= 0 && y < oldHeight && y >= 0) {
                newGrid[x + y * window_width] = grid[x + y * oldWidth];
                newcolor[x + y * window_width] = color[x + y * oldWidth];
                newzBuffer[x + y * window_width] = zBuffer[x + y * oldWidth];
            }
            else {
                newGrid[x + y * window_width] = false;
                newcolor[x + y * window_width] = vec3(0);
                newzBuffer[x + y * window_width] = -99999999999.9999999999;
            }
        }
    }

    delete[] grid;
    delete[] color;
    delete[] zBuffer;
    grid = newGrid;
    color = newcolor;
    zBuffer = newzBuffer;

    glViewport(0, 0, window_width, window_height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(0.0, window_width, 0.0, window_height);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glutPostRedisplay();
}

void update(int);

int Round(float a) {
    return (int)(a + 0.5);
}

void initcanvas(int argc, char** argv) {

    // Do drawing here
    zBuffer = new float[window_height * window_width];
    grid = new bool[window_height * window_width];
    color = new vec3[window_height * window_width];
    for (GLint x = 0; x < window_height * window_width; ++x) {
        grid[x] = false;
        color[x] = vec3(0);
        zBuffer[x] = -(99999999999.9999999999);
    }
}

void putpixel(int x, int y, const vec3& col = 1) {
    if (x < window_width && x >= 0 && y < window_height && y >= 0) {
        color[x + y * window_width] = col;
        grid[x + y * window_width] = true;
    }
}

void zbufferputpixel(int x, int y, float zBuf, const vec3& col = 1) {
    window_width = (int)window_width;
    if (x < window_width && x >= 0 && y < window_height && y >= 0) {
        if (zBuffer[x + y * window_width] <= zBuf)
        {
            color[x + y * window_width] = col;
            grid[x + y * window_width] = true;
            zBuffer[x + y * window_width] = zBuf;
        }

    }

}

void putpixel_adjusted(int x, int y, const vec3_T<float>& col = 0) {
    putpixel(x + window_width / 2, y + window_height / 2, col);
}

void putpixel_adjusted2(int x, int y, const vec3_T<float>& col = 0) {
    putpixel(x, -y + window_height, col);
}

void putpixel_adjusted2(vertex P, const vec3_T<float>& col = 0) {
    vec3_T<float> color = { P.r,P.g,P.b };
    putpixel(P.x + 0.5f, -P.y + 0.5f + window_height, color);
    //cout << "Color" << col;
}

void putpixel_adjusted3(int x, int y, int zbuf, const vec3_T<float>& col = 0) {

    zbufferputpixel(x, -y + window_height, zbuf, col);
    //cout << "Color" << col;
}

int getmaxx() {
    return 800;

}

int getmaxy() {
    return 600;
}






