#pragma once


#include"graphics.h"

#include<iostream>
#include<fstream>
#include<strstream>
#include<chrono>
#include<thread>
#include<math.h>
#include"structures.h"
#include"matrices.h"
#include"mathtools.h"
#include<thread>
#include <vector>
#include<algorithm>
#include <time.h>
#include"camera.h"
#include "modelparser.h"

//using namespace std;

void displaymatrix3(float mat[3][3], int n) {
	cout << "Matrix Display" << endl;
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < n; j++) {
			cout << mat[i][j] << "\t";

		}
		cout << endl;
	}
	cout << endl;

}

class poly3d {
private:
	mesh mesh;
	mat4d mat;
	ModelParse* object;



public:
	/*vertex vCamera = { 0,0,-10 };
	cam vCam;
	*/

	poly3d() {
		//mat = projectionMatrix();
		mat = projectionMatrix();


		//-------------------------------Assignment 2-------------------------//
		//mat = obliquematrix(90,90,1);
		//---------------------------------------------------------//


		//-------------------------------Assignment 3-------------------------//
		//mat = prespectivematrix(0,0,10,0,0);

		//---------------------------------------------------------//






		/*mat.m[0][0] = aspectRatio * fieldofView;
		mat.m[1][1] = fieldofView;
		mat.m[2][2] = zFar / (zFar - zNear);
		mat.m[3][2] = (-zFar * zNear) / (zFar - zNear);
		mat.m[2][3] = 1.0f;
		mat.m[3][3] = 0.0f;
		//*/
		//	mesh.triangles = {

		//	// south
		//	{ 0.0f, 0.0f, 0.0f,    0.0f, 1.0f, 0.0f,    1.0f, 1.0f, 0.0f },
		//	{ 0.0f, 0.0f, 0.0f,    1.0f, 1.0f, 0.0f,    1.0f, 0.0f, 0.0f },

		//	// east                                                      
		//	{ 1.0f, 0.0f, 0.0f,    1.0f, 1.0f, 0.0f,    1.0f, 1.0f, 1.0f },
		//	{ 1.0f, 0.0f, 0.0f,    1.0f, 1.0f, 1.0f,    1.0f, 0.0f, 1.0f },

		//	// north                                                     
		//	{ 1.0f, 0.0f, 1.0f,    1.0f, 1.0f, 1.0f,    0.0f, 1.0f, 1.0f },
		//	{ 1.0f, 0.0f, 1.0f,    0.0f, 1.0f, 1.0f,    0.0f, 0.0f, 1.0f },

		//	// west                                                      
		//	{ 0.0f, 0.0f, 1.0f,    0.0f, 1.0f, 1.0f,    0.0f, 1.0f, 0.0f },
		//	{ 0.0f, 0.0f, 1.0f,    0.0f, 1.0f, 0.0f,    0.0f, 0.0f, 0.0f },

		//	// top                                                       
		//	{ 0.0f, 1.0f, 0.0f,    0.0f, 1.0f, 1.0f,    1.0f, 1.0f, 1.0f },
		//	{ 0.0f, 1.0f, 0.0f,    1.0f, 1.0f, 1.0f,    1.0f, 1.0f, 0.0f },

		//	// bottom                                                    
		//	{ 1.0f, 0.0f, 1.0f,    0.0f, 0.0f, 1.0f,    0.0f, 0.0f, 0.0f },
		//	{ 1.0f, 0.0f, 1.0f,    0.0f, 0.0f, 0.0f,    1.0f, 0.0f, 0.0f },

		//};


		//mesh.LoadFromObjectFile("pumpkin.txt");
		object = new ModelParse("TestObj1.txt");




	}

	float fTheta;
	void pipeline(int i, int state, cam& vcam, ptlight& l1, float tranvalue) {
		float angle = i;
		/*mat4d matRotX = { 1,0,0,0,0,cosf(angle),-sinf(angle),0,0,sinf(angle),cosf(angle),0,0,0,0,1 };
		mat4d matRotZ = { cos(angle),-sinf(angle),0,0,sinf(angle),cosf(angle),0,0,0,0,1,0,0,0,0,1 };
		mat4d matTran = { 1,0,0,400,0,0,0,300,0,0,0,0,0,0,0,1 };
		mat4d matScale = { 50,0,0,0,0,50,0,0,0,0,50,0,0,0,0,1 };*/


		// Set up rotation matrices
		mat4d matRotZ, matRotX, matRotY, matTran, matTranXY, scaleMatrix, vscaleMatrix, vtranMatrix, flipymatrix,halfscale,fiftyscale,translatetocenter,scale20,shearZ;

		fTheta = 1.0 * (i / 10.0f);

		// Rotation Z
		matRotZ = RotationMatrixZ(fTheta);
		
		// Rotation X
		matRotX = RotationMatrixX(fTheta);

		matRotY = RotationMatrixY(fTheta);

		//Translation
		matTran = TranslationMatrix(0.0f, 0.0f, tranvalue + 5.0f);

		//Scale by 1/2
		halfscale = ScaleMatrix(0.5f, 0.5f, 0.5f);
		fiftyscale = ScaleMatrix(50, 50, 50);
		scale20 = ScaleMatrix(20, 20, 20);
		shearZ = shearz(2, 2);

		matTranXY = TranslationMatrix(1.0f, 1.0f, 0.0f);

		scaleMatrix = ScaleMatrix(0.5f * (float)getmaxx(), 0.5f * (float)getmaxy(), 1);

		//
		translatetocenter = TranslationMatrix(0.5f * (float)getmaxx(), 0.5f * (float)getmaxy(),0);
		vscaleMatrix = ScaleMatrix(1, 1, 1);
		flipymatrix = flipy();


		vector<triangle> sortedTriangles;

		float dp;

		//camera matrix view


		mat4d viewmat = vcam.matview();
		mat4d modalmatrix=Matrix_MakeIdentity();
		/*cout << "inverse" << endl;
		displaymat(viewmat);*/

		//------------------------------------------Modoel Parse Object Implementation----------------------------------------------//
		
		mesh = getsmesh(object);
		//-------------------------------------------Implementation Ends-------------------------------------------------------------//
		int tricount = 0;
		bool flag = true;
		for (auto& tri : mesh.triangles)
		{
			//cout << i;


			triangle triProjected, triTranslated, triRotatedZ, triRotatedZX, triScaled, triCamviewed;

			// Rotate in Z-Axis
			//---------------for Assignment - 1----------------------------//
			
			MultiplyMatrixVector(tri, tri, halfscale);
		//	MultiplyMatrixVector(tri, tri, shearZ);

			/*cout << "Before:";
			cout << triTranslated.p[0] << "," << triTranslated.p[1] << "," << triTranslated.p[2]<<endl;
			MultiplyMatrixVector(triTranslated, triTranslated, fiftyscale);
			cout << "After:";
			cout << triTranslated.p[0] << "," << triTranslated.p[1] << "," << triTranslated.p[2] << endl;*/

		



			//--------------------------------------------------------------//


			if (state == 0)
				MultiplyMatrixVector(tri, triTranslated, matRotZ);
				modalmatrix = Matrix_MultiplyMatrix(modalmatrix, matRotZ);



			// Rotate in X-Axis
			if (state == 1)
				MultiplyMatrixVector(tri, triTranslated, matRotX);
				modalmatrix = Matrix_MultiplyMatrix(modalmatrix, matRotX);


			if (state == 2)
				MultiplyMatrixVector(tri, triTranslated, matRotY);
				modalmatrix = Matrix_MultiplyMatrix(modalmatrix, matRotY);


			if (state == 3) {

				MultiplyMatrixVector(tri, triRotatedZ, matRotZ);
				modalmatrix = Matrix_MultiplyMatrix(modalmatrix, matRotZ);
				MultiplyMatrixVector(triRotatedZ, triTranslated, matRotX);
				modalmatrix = Matrix_MultiplyMatrix(modalmatrix, matRotX);
			}


			// Offset into the screen

			//transformMultiplication(triTranslated, triTranslated, matTran);
			


			vertex normal, line1, line2;
			line1.x = (triTranslated.p[1].x - triTranslated.p[0].x);
			line1.y = triTranslated.p[1].y - triTranslated.p[0].y;
			line1.z = triTranslated.p[1].z - triTranslated.p[0].z;
			line2.x = triTranslated.p[2].x - triTranslated.p[0].x;
			line2.y = triTranslated.p[2].y - triTranslated.p[0].y;
			line2.z = triTranslated.p[2].z - triTranslated.p[0].z;

			normal = crossProduct(line1, line2);
			// normalise the normal
			float l = sqrtf(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z);
			normal.x /= l; normal.y /= l; normal.z /= l;

	

			if ((normal.x * (triTranslated.p[0].x - vcam.vCamera.x) +
				normal.y * (triTranslated.p[0].y - vcam.vCamera.y) +
				normal.z * (triTranslated.p[0].z - vcam.vCamera.z) < 0.0f))
			{
				//cout << "V" << vCamera.x;
				// Illumination
				vertex light_direction = { 1.0f, 1.0f, 1.0f };
				/*vertex light_direction = vcam.lookdir();*/
				float l = sqrtf(light_direction.x * light_direction.x + light_direction.y * light_direction.y + light_direction.z * light_direction.z);
				light_direction.x /= l; light_direction.y /= l; light_direction.z /= l;


				//-----------------------Point Source Implementation--------------------------------//

					//l1.caclulateintensity(tri,triTranslated);
					l1.calculatephongintensity(tri, triTranslated,vcam,modalmatrix);
			/*	if (flag)
					l1.calculatefaceintensity(tri, triTranslated,matRotZ);
				flag = false;*/
				
				//tricount++;


				//----------------------Point Source  Implementation Ends---------------------------//

				// How similar is normal to light direction
				dp = normal.x * light_direction.x + normal.y * light_direction.y + normal.z * light_direction.z;


				MultiplyMatrixVector(triTranslated, triCamviewed, viewmat);

				// Project triangles from 3D --> 2D

				MultiplyMatrixVector(triCamviewed, triProjected, mat);

				//transformMultiplication(triProjected, triProjected, vscaleMatrix);
				// Scale into view
				/*triProjected.p[0].x += 1.0f; triProjected.p[0].y += 1.0f;
				triProjected.p[1].x += 1.0f; triProjected.p[1].y += 1.0f;
				triProjected.p[2].x += 1.0f; triProjected.p[2].y += 1.0f;*/

				//translates it by 1.0f.
				transformMultiplication(triProjected, triProjected, matTranXY);

				//Scales the cube to bring it to the  origin

				//---------------Pay Attention to this piece of code---------------------------//

				transformMultiplication(triProjected, triProjected, scaleMatrix);
				//transformMultiplication(triProjected, triProjected, translatetocenter);

				//transformMultiplication(triProjected, triProjected, flipymatrix);



				//transformMultiplication(triProjected, triProjected, vtranMatrix);
				triProjected.p->dp = dp;

				//intensity implementation for point source

				l1.getintensity(triProjected); //sets intensity value in triprojected triangles.


				sortedTriangles.push_back(triProjected);

				/*drawTriangle(triProjected.p[0].x, triProjected.p[0].y,
					triProjected.p[1].x, triProjected.p[1].y,
					triProjected.p[2].x, triProjected.p[2].y);*/
					//	this_thread::sleep_for(chrono::milliseconds(20));



			}




			i = i % 360;

			//mesh.triangles.clear();
		}
		//Painter's algorithm
		sort(sortedTriangles.begin(), sortedTriangles.end(), [](triangle& t1, triangle& t2)
			{
				float z1 = (t1.p[0].z + t1.p[1].z + t1.p[2].z) / 3.0f;
				float z2 = (t2.p[0].z + t2.p[1].z + t2.p[2].z) / 3.0f;
				return z1 > z2;
			});

		for (auto& triProjected : sortedTriangles) {
			float intensity[3] = { triProjected.p[0].intensity ,triProjected.p[1].intensity ,triProjected.p[2].intensity };
			rasterTriangle(triProjected.p[0].x, triProjected.p[0].y, triProjected.p[0].z,
				triProjected.p[1].x, triProjected.p[1].y, triProjected.p[1].z,
				triProjected.p[2].x, triProjected.p[2].y, triProjected.p[2].z, triProjected.p[0].dp, intensity);


		}
		std::this_thread::sleep_for(std::chrono::milliseconds(0));

		//}

	}

	void render(int& i, int state, cam& vcam, ptlight& l1, float tranvalue = 1.0f) {

		pipeline(i, state, vcam, l1, tranvalue);


	}



	





	void drawTriangle(float x1, float y1, float x2, float y2, float x3, float y3) {
		DDAlgorithm(x1, y1, x2, y2, 1.0f);
		DDAlgorithm(x2, y2, x3, y3, 1.0f);
		DDAlgorithm(x3, y3, x1, y1, 1.0f);

	}

	void rasterTriangle(float& x1, float& y1, float& z1, float& x2, float& y2, float& z2, float& x3, float& y3, float& z3, float& dp, float* intensity) {
		//vector<vertex> vert[3];

			//drawTriangle(x1, y1, x2, y2, x3, y3);

		vertex v1{ x1,y1,z1 }; vertex v2{ x2,y2,z2 }; vertex v3{ x3,y3,z3 };
		v1.r = intensity[0] * 1.0f;//((float)rand() / RAND_MAX);
		v1.g = intensity[0] * 1.0f;//((float)rand() / RAND_MAX);
		v1.b = intensity[0] * 1.0f;//((float)rand() / RAND_MAX);
		v2.r = intensity[1] * 1.0f;// ((float)rand() / RAND_MAX);
		v2.g = intensity[1] * 1.0f;//((float)rand() / RAND_MAX);
		v2.b = intensity[1] * 1.0f;//((float)rand() / RAND_MAX);
		v3.r = intensity[2] * 1.0f;//((float)rand() / RAND_MAX);
		v3.g = intensity[2] * 1.0f;//((float)rand() / RAND_MAX);
		v3.b = intensity[1] * 1.0f;//((float)rand() / RAND_MAX);
		v1.dp = dp;
		v2.dp = dp;
		v3.dp = dp;

		//fillTriangle(v1, v2, v3, dp);
		//	grouradFiller(v1, v2, v3, dp);
			//brasterize(v1, v2, v3, dp);

		vertex arr[3] = { v1,v2,v3 };
		vec3 col = { 0.0f ,1,1 };
		//float intensity[3] = { dp,dp,dp };
	//	float inten[3] = { 0.5f+intensity[0],0.5f+intensity[1],0.5f+intensity[3] };
		/*cout <<"I1: "<< intensity[0];
		cout << "I2: " << intensity[1];
		cout << "I3: " << intensity[2];*/

		//-------------------------------------for Assignment 2-------------------------//

		//drawTriangle(v1.x, v1.y, v2.x, v2.y, v3.x, v3.y);

		//------------------------------------------------------------------------------//


		bellatriangle(arr, 0, col, intensity);
		dp = intensity[0];
		//fillTriangle(v1, v2, v3, dp);
		//trasterize(v1, v2, v3, dp);
		//fillTriangle(v1.x, v1.y, v2.x, v2.y,v3.x,v3.y);

	}

};
//vertex poly3d::vCamera = { 0,0,0 };