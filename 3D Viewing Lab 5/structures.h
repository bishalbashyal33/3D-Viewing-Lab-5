#pragma once
#pragma once
#include<iostream>
using namespace std;
#include<vector>
#include<fstream>
#include<strstream>
#include"graphics.h"
#include"maths.h"

struct vertex
{
	float x, y, z;
	float dp = -1.0f;
	float w = 0;
	float r = 1.0f;
	float g = 1.0f;
	float b = 1.0f;
	float intensity = 0.0f;
	vec3 normal = vec3(0);

	vertex operator-(vertex& vert) {
		vertex temp;
		temp.x = x - vert.x;
		temp.y = y - vert.y;
		temp.z = z - vert.z;
		return temp;
	}

	vertex operator* (float k) {
		vertex temp;
		temp.x = x * k;
		temp.y = y * k;
		temp.z = z * k;
		return temp;
	}

	vertex operator+ (vertex& vert) {
		vertex temp;
		temp.x = x + vert.x;
		temp.y = y + vert.y;
		temp.z = z + vert.z;
		return temp;
	}

	

	bool operator==(vertex& vert) {
		if (vert.x == x && vert.y == y && vert.z == z)return true;
		return false;
	}

	friend ostream& operator<<(ostream& output,vertex& v) {
		output << "(" << v.x << " , " << v.y<<","<<v.z<<')';
		return output;
	}




};

struct triangle
{
	vertex p[3];
};




struct mesh
{

	vector<triangle> triangles;

	bool LoadFromObjectFile(string sFilename)
	{
		ifstream f;

		f.open(sFilename, ios::out | ios::in);
		if (!f.is_open()) {
			cout << "File not opened";
			return false;

		}
		// Local cache of verts
		vector<vertex> verts;

		while (!f.eof())
		{
			char line[128];
			f.getline(line, 128);

			strstream s;
			s << line;

			char junk;

			if (line[0] == 'v')
			{
				vertex v;
				s >> junk >> v.x >> v.y >> v.z;
				verts.push_back(v);
			}

			if (line[0] == 'f')
			{
				int f[3];
				s >> junk >> f[0] >> f[1] >> f[2];

				triangles.push_back({ verts[f[0] - 1], verts[f[1] - 1], verts[f[2] - 1] });
			}
			if (line[0] == 'v' && line[1] == 'n') {

			}
		}

		return true;
	}

};

struct mat4d
{

	float m[4][4] = { 0 };




};

//Scanline algorithm







