#pragma once

#pragma once
#include<math.h>
#include "mathtools.h"
#include"structures.h"


float zNear = 0.1f;
float zFar = 100.0f;
float angleofview = 90.0f;
float aspectRatio = 1;//float(float(8) / 6);
float fieldofView = 1.0f / tanf(angleofview * 0.5f / 180.0f * 2*acosf(0.0f));
//mat4d projectionMatrix({aspectRatio * fieldofView,fieldofView,zFar / (zFar - zNear),(-zFar * zNear) / (zFar - zNear),1.0f,0.0f });


//  float angle = 120;


// namespace mt

auto projectionMatrix() {
	mat4d mat;
	mat.m[0][0] = aspectRatio * fieldofView;
	mat.m[1][1] = fieldofView;
	mat.m[2][2] = zFar / (zFar - zNear);
	mat.m[3][2] = (-zFar * zNear) / (zFar - zNear);
	mat.m[2][3] = 1.0f;
	mat.m[3][3] = 0.0f;
	return mat;
}

auto flipy() {
	mat4d mat;
	mat.m[0][0] = 1;
	mat.m[1][1] = 1;
	mat.m[2][2] = -1;
	mat.m[3][3] = 1;
	return mat;
}


auto shearz(float a,float b) {
	mat4d mat;
	mat.m[0][0] = 1;
	mat.m[1][1] = 1;
	mat.m[2][0] = a;
	mat.m[2][1] = b;
	mat.m[2][2] = 1;
	mat.m[3][3] = 1;

	return mat;
}






auto RotationMatrixZ(float fTheta) {
	mat4d matRotZ;
	matRotZ.m[0][0] = cosf(fTheta);
	matRotZ.m[0][1] = sinf(fTheta);
	matRotZ.m[1][0] = -sinf(fTheta);
	matRotZ.m[1][1] = cosf(fTheta);
	matRotZ.m[2][2] = 1.0f;
	matRotZ.m[3][3] = 1.0f;
	return matRotZ;
}

auto RotationMatrixX(float fTheta) {
	mat4d matRotX;
	matRotX.m[0][0] = 1.0f;
	matRotX.m[1][1] = cosf(fTheta * 0.5f);
	matRotX.m[1][2] = sinf(fTheta * 0.5f);
	matRotX.m[2][1] = -sinf(fTheta * 0.5f);
	matRotX.m[2][2] = cosf(fTheta * 0.5f);
	matRotX.m[3][3] = 1.0f;
	return matRotX;
}

auto RotationMatrixY(float fTheta) {

	mat4d matRotX;
	matRotX.m[0][0] = cosf(fTheta);
	matRotX.m[1][1] = 1.0f;
	matRotX.m[0][2] = sinf(fTheta);
	matRotX.m[2][0] = -sinf(fTheta);
	matRotX.m[2][2] = cosf(fTheta);
	matRotX.m[3][3] = 1.0f;
	return matRotX;

}

/*
mat.m[0][0] =
mat.m[0][1] =
mat.m[0][2] =
mat.m[0][3] =
mat.m[1][0] =
mat.m[1][0] =

*/

auto obliquematrix(float alpha, float phi, int col = 0) {
	
	alpha = (2 * acosf(0.0f) / 180) * (alpha);
	phi = (2 * acosf(0.0f) / 180) * (phi);
	float L1 = 1 / tan(alpha);

	mat4d mat;
	mat.m[0][0] = 1;
	mat.m[1][1] = 1;
	mat.m[2][0] = L1 * cos(phi);
	mat.m[2][1] = L1 * sin(phi);
	mat.m[3][3] = 1;

	return mat;


}

auto prespectivematrix(float xprp, float yprp, float zprp, float zvp, int col = 0) {
	float dp = zprp - zvp;
	mat4d mat;
	mat.m[0][0] = 1;
	mat.m[1][1] = 1;
	mat.m[2][0] = xprp / dp;
	mat.m[2][1] = yprp / dp;
	mat.m[2][2] = -zvp / dp;
	mat.m[2][3] = -1/dp;
	mat.m[3][0] = -xprp * zvp / dp;
	mat.m[3][1] = -yprp * zvp / dp;
	mat.m[3][2] = zprp * zvp / dp;
	mat.m[3][3] = zprp / dp;
	return mat;


}






auto TranslationMatrix(float tx, float ty, float tz) {
	mat4d matTran;
	matTran.m[0][0] = 1;
	matTran.m[0][1] = 0;
	matTran.m[0][2] = 0;
	matTran.m[1][1] = 1;
	matTran.m[1][2] = 0;
	matTran.m[2][2] = 1;
	matTran.m[0][3] = tx;
	matTran.m[1][3] = ty;
	matTran.m[2][3] = tz;
	matTran.m[3][3] = 1;
	return matTran;
}


auto ScaleMatrix(float sx, float sy, float sz) {
	mat4d mat;
	mat.m[0][0] = sx;
	mat.m[1][1] = sy;
	mat.m[2][2] = sz;
	mat.m[3][3] = 1.0f;
	return mat;

}

auto OrthoMatrix() {
	mat4d mat;
	mat.m[0][0] = 1;
	mat.m[1][1] = 1;
	mat.m[2][2] = 1;
	mat.m[3][3] = 1;
	return mat;

}

//auto lookAtMatrix(vertex target, vertex pos, vertex up) {
//	
//
//	mat4d mat;
//	mat.m[0][0] =A.x ;
//	mat.m[0][1] =B.x ;
//	mat.m[0][2] =C.x ;
//	mat.m[1][0] = A.y;
//	mat.m[1][1] = B.y;
//	mat.m[1][2] = C.y;
//	mat.m[2][0] = A.z ;
//	mat.m[2][1] = B.z;
//	mat.m[2][2] = C.z;
//	mat.m[3][0] = -dotProduct(T,A);
//	mat.m[3][1] = -dotProduct(T, B);
//	mat.m[3][2] = -dotProduct(T, C);
//	mat.m[3][3] = 1;
//	return mat;
//
//}

auto pointAtMatrix(vertex& pos, vertex& target, vertex& up) {
	// Calculate new forward direction
	vertex newForward = target - pos;
	newForward = normalize(newForward);

	// Calculate new Up direction
	vertex a = newForward * dotProduct(up, newForward);
	vertex newUp = up - a;
	newUp = normalize(newUp);

	// New Right direction is easy, its just cross product
	vertex newRight = crossProduct(newUp, newForward);
	/*mat4d matrix;
	matrix.m[0][0] = newRight.x;	matrix.m[0][1] = newRight.y;	matrix.m[0][2] = newRight.z;	matrix.m[0][3] = 0.0f;
	matrix.m[1][0] = newUp.x;		matrix.m[1][1] = newUp.y;		matrix.m[1][2] = newUp.z;		matrix.m[1][3] = 0.0f;
	matrix.m[2][0] = newForward.x;	matrix.m[2][1] = newForward.y;	matrix.m[2][2] = newForward.z;	matrix.m[2][3] = 0.0f;
	matrix.m[3][0] = pos.x;			matrix.m[3][1] = pos.y;			matrix.m[3][2] = pos.z;			matrix.m[3][3] = 1.0f;
	return matrix;*/

	mat4d matrix;
	matrix.m[0][0] = newRight.x;	matrix.m[0][1] = newRight.y;	matrix.m[0][2] = newRight.z;	matrix.m[0][3] = 0.0f;
	matrix.m[1][0] = newUp.x;		matrix.m[1][1] = newUp.y;		matrix.m[1][2] = newUp.z;		matrix.m[1][3] = 0.0f;
	matrix.m[2][0] = newForward.x;	matrix.m[2][1] = newForward.y;	matrix.m[2][2] = newForward.z;	matrix.m[2][3] = 0.0f;
	matrix.m[3][0] = pos.x;			matrix.m[3][1] = pos.y;			matrix.m[3][2] = pos.z;			matrix.m[3][3] = 1.0f;
	return matrix;



}

mat4d Matrix_Inverse(mat4d& m) // Only for Rotation/Translation Matrices
{
	mat4d matrix;
	matrix.m[0][0] = m.m[0][0]; matrix.m[0][1] = m.m[1][0]; matrix.m[0][2] = m.m[2][0]; matrix.m[0][3] = 0.0f;
	matrix.m[1][0] = m.m[0][1]; matrix.m[1][1] = m.m[1][1]; matrix.m[1][2] = m.m[2][1]; matrix.m[1][3] = 0.0f;
	matrix.m[2][0] = m.m[0][2]; matrix.m[2][1] = m.m[1][2]; matrix.m[2][2] = m.m[2][2]; matrix.m[2][3] = 0.0f;
	matrix.m[3][0] = -(m.m[3][0] * matrix.m[0][0] + m.m[3][1] * matrix.m[1][0] + m.m[3][2] * matrix.m[2][0]);
	matrix.m[3][1] = -(m.m[3][0] * matrix.m[0][1] + m.m[3][1] * matrix.m[1][1] + m.m[3][2] * matrix.m[2][1]);
	matrix.m[3][2] = -(m.m[3][0] * matrix.m[0][2] + m.m[3][1] * matrix.m[1][2] + m.m[3][2] * matrix.m[2][2]);
	matrix.m[3][3] = 1.0f;
	return matrix;
}







// Rotation X
