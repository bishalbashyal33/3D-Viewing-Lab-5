#pragma once

#pragma once
#include<math.h>
#include"graphics.h"
#include"structures.h"
#include<iostream>
#include<limits>
#include"tools.h"
#include"modelparser.h"
using namespace std;

void MultiplyMatrixVector(vertex& in, vertex& out, mat4d& m)
{


    out.x = in.x * m.m[0][0] + in.y * m.m[1][0] + in.z * m.m[2][0] + m.m[3][0];
    out.y = in.x * m.m[0][1] + in.y * m.m[1][1] + in.z * m.m[2][1] + m.m[3][1];
    out.z = in.x * m.m[0][2] + in.y * m.m[1][2] + in.z * m.m[2][2] + m.m[3][2];
    float w = in.x * m.m[0][3] + in.y * m.m[1][3] + in.z * m.m[2][3] + m.m[3][3];
    if (w != 0.0f)
    {
        out.x /= w; out.y /= w; out.z /= w;
    }
}

void MultiplyMatrixVector(triangle& in, triangle& out, mat4d& m)
{
    for (int i = 0; i < 3; i++) {
        out.p[i].x = in.p[i].x * m.m[0][0] + in.p[i].y * m.m[1][0] + in.p[i].z * m.m[2][0] + m.m[3][0];
        out.p[i].y = in.p[i].x * m.m[0][1] + in.p[i].y * m.m[1][1] + in.p[i].z * m.m[2][1] + m.m[3][1];
        out.p[i].z = in.p[i].x * m.m[0][2] + in.p[i].y * m.m[1][2] + in.p[i].z * m.m[2][2] + m.m[3][2];
        float w = in.p[i].x * m.m[0][3] + in.p[i].y * m.m[1][3] + in.p[i].z * m.m[2][3] + m.m[3][3];
        if (w != 0.0f)
        {
            out.p[i].x /= w; out.p[i].y /= w; out.p[i].z /= w;
        }
    }
}

vertex Matrix_MultiplyVector(mat4d& m, vertex& i)
{
    vertex v;
    v.x = i.x * m.m[0][0] + i.y * m.m[1][0] + i.z * m.m[2][0] + i.w * m.m[3][0];
    v.y = i.x * m.m[0][1] + i.y * m.m[1][1] + i.z * m.m[2][1] + i.w * m.m[3][1];
    v.z = i.x * m.m[0][2] + i.y * m.m[1][2] + i.z * m.m[2][2] + i.w * m.m[3][2];
    v.w = i.x * m.m[0][3] + i.y * m.m[1][3] + i.z * m.m[2][3] + i.w * m.m[3][3];
    return v;
}

vertex Vector_Div(vertex& v1, float k)
{
    return { v1.x / k, v1.y / k, v1.z / k };
}

void transformMultiplication(triangle& in, triangle& out, mat4d& m) {
    for (int i = 0; i < 3; i++) {
        out.p[i].x = in.p[i].x * m.m[0][0] + in.p[i].y * m.m[0][1] + in.p[i].z * m.m[0][1] + m.m[0][3];
        out.p[i].y = in.p[i].x * m.m[1][0] + in.p[i].y * m.m[1][1] + in.p[i].z * m.m[1][2] + m.m[1][3];
        out.p[i].z = in.p[i].x * m.m[2][0] + in.p[i].y * m.m[2][1] + in.p[i].z * m.m[2][2] + m.m[2][3];
        // float w = in.p[i].x * m.m[0][3] + in.p[i].y * m.m[1][3] + in.p[i].z * m.m[2][3] + m.m[3][3];

    }

}


void MultiplyMatrixVector(triangle& in, triangle& out, mat4d& m, int n)
{
    for (int i = 0; i < 3; i++) {
        out.p[i].x = in.p[i].x * m.m[0][0] + in.p[i].y * m.m[1][0] + in.p[i].z * m.m[2][0] + m.m[3][0];
        out.p[i].y = in.p[i].x * m.m[0][1] + in.p[i].y * m.m[1][1] + in.p[i].z * m.m[2][1] + m.m[3][1];
        out.p[i].z = in.p[i].x * m.m[0][2] + in.p[i].y * m.m[1][2] + in.p[i].z * m.m[2][2] + m.m[3][2];


    }
}

mat4d Matrix_MultiplyMatrix(mat4d& m1, mat4d& m2)
{
    mat4d matrix;
    for (int c = 0; c < 4; c++)
        for (int r = 0; r < 4; r++)
            matrix.m[r][c] = m1.m[r][0] * m2.m[0][c] + m1.m[r][1] * m2.m[1][c] + m1.m[r][2] * m2.m[2][c] + m1.m[r][3] * m2.m[3][c];
    return matrix;
}

void DDAlgorithm(float x1, float y1, float x2, float y2, float dp) {
    float dx, dy, xinc, yinc, step;
    dx = x2 - x1;
    dy = y2 - y1;
    if (abs(dx) > abs(dy))
    {
        step = abs(dx);
    }
    else {
        step = abs(dy);
    }
    xinc = dx / step;
    yinc = dy / step;
    int count = 0;
    //cout << "From DDA" << endl;
    while (count <= step) {
        /*  putpixel_adjusted(round(x1), round(y1),vec3(0,0,dp));
          putpixel(round(x1), round(y1), vec3(0, 0, dp));*/
        putpixel_adjusted2(round(x1), round(y1), vec3(1, 1, dp));


        // cout << "(" << round(x1) << "," << round(y1) << ")" << endl;
        x1 = x1 + xinc;
        y1 = y1 + yinc;
        count++;
    }
}


void Bresenham(int x1, int y1, int x2, int y2) {

    int lx, ly;
    int dx = abs(x2 - x1);
    int dy = abs(y2 - y1);
    int x = x1;
    int y = y1;
    x2 > x1 ? lx = 1 : lx = -1;
    y2 > y1 ? ly = 1 : ly = -1;
    int p = 2 * dy - dx;
    int p0 = 2 * dy - dx;
    int pk = p0;

    if (dx > dy) {


        int count = 0;
        while (count <= dx) {
            if (pk < 0) {

                putpixel(x, y, -1);
                x = x + lx;
                pk = pk + 2 * dy;
            }

            else {
                putpixel(x, y, -1);
                x = x + lx;
                y = y + ly;
                pk = pk + 2 * dy - 2 * dx;
            }
            count++;
        }


    }
    else {
        p0 = 2 * dx - dy;
        pk = p0;
        int count = 0;
        while (count <= dy) {
            if (pk < 0) {

                putpixel(x, y, 255);
                y = y + ly;
                pk = pk + 2 * dx;

            }
            else {

                putpixel(x, y, 255);
                x = x + lx;
                y = y + ly;
                pk = pk + 2 * dx - 2 * dy;
            }
            count++;
        }
    }
}

//


//void fillBottomFlatTriangle(vertex v1, vertex v2, vertex v3, float& dp)
//{
//
//    float invslope1 = (v2.x - v1.x) / (v2.y - v1.y);
//    // cout <<"a"<< (v2.x - v1.x) / (v2.y - v1.y)<<endl;
//    float invslope2 = (v3.x - v1.x) / (v3.y - v1.y);
//
//    float curx1 = v1.x;
//    float curx2 = v1.x;
//    const int ystart = (int)ceil(v1.y - 0.5f);
//    const int yend = (int)ceil(v3.y - 0.5f);
//
//    for (int y = ystart; y < yend; y++)
//    {
//        const float px0 = invslope1 * (float(y) + 0.5f - v1.y) + v1.x;
//        const float px1 = invslope2 * (float(y) + 0.5f - v2.y) + v2.x;
//
//        const int xstart = (int)ceil(px0 - 0.5f);
//        const int xend = (int)ceil(px1 - 0.5f);
//
//        for (int x = xstart; x < xend; x++) {
//            putpixel_adjusted2(x, y, dp);
//        }
//    }
//
//  
//}
//
//void fillTopFlatTriangle(vertex& v1, vertex& v2, vertex& v3, float& dp)
//{
//
//    float invslope1 = (v3.x - v1.x) / (v3.y - v1.y);
//    float invslope2 = (v3.x - v2.x) / (v3.y - v2.y);
//
//    float curx1 = v3.x;
//    float curx2 = v3.x;
//    const int ystart = (int)ceil(v1.y - 0.5f);
//    const int yend = (int)ceil(v3.y - 0.5f);
//
//    for (int y = ystart; y <yend; y++)
//    {
//        const float px0 = invslope1 * (float(y) + 0.5f - v1.y) + v1.x;
//        const float px1 = invslope2 * (float(y) + 0.5f - v2.y) + v2.x;
//
//        const int xstart = (int)ceil(px0 - 0.5f);
//        const int xend = (int)ceil(px1 - 0.5f);
//
//        for (int x = xstart; x < xend; x++) {
//            putpixel_adjusted2(x, y, dp);
//        }
//    }
//}
//
//
//
//





void fillBottomFlatTriangle(vertex v1, vertex v2, vertex v3, float& dp)
{

    float invslope1 = (v2.x - v1.x) / (v2.y - v1.y);
    // cout <<"a"<< (v2.x - v1.x) / (v2.y - v1.y)<<endl;
    float invslope2 = (v3.x - v1.x) / (v3.y - v1.y);

    float curx1 = v1.x;
    float curx2 = v1.x;
    DDAlgorithm(v1.x, v1.y, v2.x, v2.y, dp);
    DDAlgorithm(v2.x, v2.y, v3.x, v3.y, dp);
    DDAlgorithm(v3.x, v3.y, v1.x, v1.y, dp);


    for (int scanlineY = v1.y; scanlineY < v2.y - 0.5f; scanlineY++)
    {
        if (scanlineY == v2.y) {
            cout << "x" << curx1 << "," << curx2 << endl;
            cout << v2.y << endl;
        }
        DDAlgorithm(curx1, scanlineY, curx2, scanlineY, dp);
        curx1 += invslope1;
        curx2 += invslope2;


    }

}

void fillTopFlatTriangle(vertex& v1, vertex& v2, vertex& v3, float& dp)
{

    float invslope1 = (v3.x - v1.x) / (v3.y - v1.y);
    float invslope2 = (v3.x - v2.x) / (v3.y - v2.y);

    float curx1 = v3.x;
    float curx2 = v3.x;
    DDAlgorithm(v1.x, v1.y, v2.x, v2.y, dp);
    DDAlgorithm(v2.x, v2.y, v3.x, v3.y, dp);
    DDAlgorithm(v3.x, v3.y, v1.x, v1.y, dp);

    for (int scanlineY = v3.y; scanlineY > v1.y; scanlineY--)
    {
        DDAlgorithm(curx1, scanlineY, curx2, scanlineY, dp);
        curx1 -= invslope1;
        curx2 -= invslope2;

    }
}




void fillTriangle(vertex v1, vertex v2, vertex v3, float& dp)
{
    vertex arr[3] = { v1,v2,v3 };
    if (arr[0].y > arr[1].y) { swap(arr[1], arr[0]); }
    if (arr[0].y > arr[2].y) { swap(arr[2], arr[0]); }
    if (arr[1].y > arr[2].y) { swap(arr[2], arr[1]); }



    if (int(arr[1].y) == int(arr[2].y))
    {

        fillBottomFlatTriangle(arr[0], arr[1], arr[2], dp);
    }
    /* check for trivial case of top-flat triangle */
    else if (int(arr[0].y) == int(arr[1].y))
    {
        fillTopFlatTriangle(arr[0], arr[1], arr[2], dp);
    }
    else
    {
        /* general case - split the triangle in a topflat and bottom-flat one */
        vertex* v4 = new vertex({
        (arr[0].x + ((float)(arr[1].y - arr[0].y) / (float)(arr[2].y - arr[0].y)) * (arr[2].x - arr[0].x)), arr[1].y,0 });
        fillBottomFlatTriangle(arr[0], arr[1], *v4, dp);
        fillTopFlatTriangle(arr[1], *v4, arr[2], dp);

    }
}

void filledTriangle(vertex A, vertex B, vertex C, float& dp) {
    vertex arr[3] = { A,B,C };
    if (arr[0].y > arr[1].y) { swap(arr[1], arr[0]); }
    if (arr[0].y > arr[2].y) { swap(arr[2], arr[0]); }
    if (arr[1].y > arr[2].y) { swap(arr[2], arr[1]); }
    A = arr[0];
    B = arr[1];
    C = arr[2];
    float dx1, dx2, dx3;
    if (B.y - A.y > 0)  dx1 = (B.x - A.x) / (B.y - A.y); else dx1 = 0;
    if (C.y - A.y > 0)  dx2 = (C.x - A.x) / (C.y - A.y); else dx2 = 0;
    if (C.y - B.y > 0)  dx3 = (C.x - B.x) / (C.y - B.y); else dx3 = 0;

    vertex S, E;
    S = A;
    E = A;
    if (dx1 > dx2) {
        for (; S.y <= B.y; S.y++, E.y++, S.x += dx2, E.x += dx1)
            DDAlgorithm(S.x, S.y, E.x, S.y, dp);
        E = B;
        for (; S.y <= C.y; S.y++, E.y++, S.x += dx2, E.x += dx3)
            DDAlgorithm(S.x, S.y, E.x, S.y, dp);;
    }
    else {
        for (; S.y <= B.y; S.y++, E.y++, S.x += dx1, E.x += dx2)
            DDAlgorithm(S.x, S.y, E.x, S.y, dp);
        S = B;
        for (; S.y <= C.y; S.y++, E.y++, S.x += dx3, E.x += dx2)
            DDAlgorithm(S.x, S.y, E.x, S.y, dp);
    }


}


void grouradFiller(vertex A, vertex B, vertex C, float& dp)

{
    vertex arr[3] = { A,B,C };
    if (arr[0].y > arr[1].y) { swap(arr[1], arr[0]); }
    if (arr[0].y > arr[2].y) { swap(arr[2], arr[0]); }
    if (arr[1].y > arr[2].y) { swap(arr[2], arr[1]); }
    A = arr[0];
    B = arr[1];
    C = arr[2];
    float dx1, dr1, dg1, db1, ddp1, ddp2, ddp3, dx2, dr2, dg2, db2, dx3, dr3, dg3, db3, dr, dg, db, ddp;
    if (B.y - A.y > 0) {
        dx1 = (B.x - A.x) / (B.y - A.y);
        dr1 = (B.r - A.r) / (B.y - A.y);
        dg1 = (B.g - A.g) / (B.y - A.y);
        db1 = (B.b - A.b) / (B.y - A.y);
        ddp1 = (B.dp - A.dp) / (B.y - A.y);
    }
    else
        dx1 = dr1 = dg1 = db1 = 0, ddp1 = 0;;

    if (C.y - A.y > 0) {
        dx2 = (C.x - A.x) / (C.y - A.y);
        dr2 = (C.r - A.r) / (C.y - A.y);
        dg2 = (C.g - A.g) / (C.y - A.y);
        db2 = (C.b - A.b) / (C.y - A.y);
        ddp2 = (C.dp - A.dp) / (C.y - A.y);
    }
    else
        dx2 = dr2 = dg2 = db2 = ddp2 = 0;

    if (C.y - B.y > 0) {
        dx3 = (C.x - B.x) / (C.y - B.y);
        dr3 = (C.r - B.r) / (C.y - B.y);
        dg3 = (C.g - B.g) / (C.y - B.y);
        db3 = (C.b - B.b) / (C.y - B.y);
        ddp3 = (C.dp - B.dp) / (C.y - B.y);
    }
    else
        dx3 = dr3 = dg3 = db3 = ddp3 = 0;
    vertex S, E, P;
    S = A;
    E = A;
    if (dx1 > dx2) {
        for (; S.y <= B.y; S.y++, E.y++) {
            if (E.x - S.x > 0) {
                dr = (E.r - S.r) / (E.x - S.x);
                dg = (E.g - S.g) / (E.x - S.x);
                db = (E.b - S.b) / (E.x - S.x);
                ddp = (E.dp - S.dp) / (E.x - S.x);
            }
            else
                dr = dg = db = ddp = 0;
            P = S;
            for (; P.x < E.x; P.x++) {
                putpixel_adjusted2(P, P.dp);
                P.r += dr; P.g += dg; P.b += db; P.dp += ddp;
            }
            S.x += dx2; S.r += dr2; S.g += dg2; S.b += db2; S.dp += ddp2;
            E.x += dx1; E.r += dr1; E.g += dg1; E.b += db1; E.dp += ddp1;
        }

        E = B;
        for (; S.y <= C.y; S.y++, E.y++) {
            if (E.x - S.x > 0) {
                dr = (E.r - S.r) / (E.x - S.x);
                dg = (E.g - S.g) / (E.x - S.x);
                db = (E.b - S.b) / (E.x - S.x);
                ddp = (E.dp - S.dp) / (E.x - S.x);
            }
            else
                dr = dg = db = 0;
            P = S;
            for (; P.x <= E.x; P.x++) {
                putpixel_adjusted2(P, P.dp);
                P.r += dr; P.g += dg; P.b += db; P.dp += ddp;
            }
            S.x += dx2; S.r += dr2; S.g += dg2; S.b += db2; S.dp += ddp2;
            E.x += dx3; E.r += dr3; E.g += dg3; E.b += db3; E.dp += ddp3;
        }
    }
    else {
        for (; S.y <= B.y; S.y++, E.y++) {
            if (E.x - S.x > 0) {
                dr = (E.r - S.r) / (E.x - S.x);
                dg = (E.g - S.g) / (E.x - S.x);
                db = (E.b - S.b) / (E.x - S.x);
                ddp = (E.dp - S.dp) / (E.x - S.x);
            }
            else
                dr = dg = db = ddp = 0;

            P = S;
            for (; P.x < E.x; P.x++) {
                putpixel_adjusted2(P, P.dp);
                P.r += dr; P.g += dg; P.b += db; P.dp += ddp;
            }
            S.x += dx1; S.r += dr1; S.g += dg1; S.b += db1; S.dp += ddp1;
            E.x += dx2; E.r += dr2; E.g += dg2; E.b += db2; E.dp += ddp2;
        }

        S = B;
        for (; S.y <= C.y; S.y++, E.y++) {
            if (E.x - S.x > 0) {
                dr = (E.r - S.r) / (E.x - S.x);
                dg = (E.g - S.g) / (E.x - S.x);
                db = (E.b - S.b) / (E.x - S.x);
                ddp = (E.dp - S.dp) / (E.x - S.x);
            }
            else
                dr = dg = db = ddp = 0;

            P = S;
            for (; P.x < E.x; P.x++) {
                putpixel_adjusted2(P, P.dp);
                P.r += dr; P.g += dg; P.b += db; P.dp += ddp;
            }
            S.x += dx3; S.r += dr3; S.g += dg3; S.b += db3; S.dp += ddp3;
            E.x += dx2; E.r += dr2; E.g += dg2; E.b += db2; E.dp += ddp2;
        }
    }
}

















mat4d Matrix_MakeIdentity()
{
    mat4d matrix;
    matrix.m[0][0] = 1.0f;
    matrix.m[1][1] = 1.0f;
    matrix.m[2][2] = 1.0f;
    matrix.m[3][3] = 1.0f;
    return matrix;
}



vertex normalize(vertex normal) {
    float l = sqrtf(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z);
    normal.x /= l; normal.y /= l; normal.z /= l;
    return normal;

}

vertex crossProduct(vertex line1, vertex line2) {
    vertex normal;
    normal.x = line1.y * line2.z - line1.z * line2.y;
    normal.y = line1.z * line2.x - line1.x * line2.z;
    normal.z = line1.x * line2.y - line1.y * line2.x;
    return normal;
}

float dotProduct(vertex A, vertex B) {
    float dp = A.x * B.x + A.y * B.y + A.z * B.z;
    return dp;
}











void trasterize(vertex t0, vertex t1, vertex t2, float& color)
{
    // sort the vertices lower?to?upper (bubblesort yay!) 
    vertex A, B, C;
    if (t0.y > t1.y) std::swap(t0, t1);
    if (t0.y > t2.y) std::swap(t0, t2);
    if (t1.y > t2.y) std::swap(t1, t2);
    int total_height = t2.y - t0.y;
    for (int i = 0; i < total_height; i++) {
        bool second_half = i > t1.y - t0.y || t1.y == t0.y;
        int segment_height = second_half ? t2.y - t1.y : t1.y - t0.y;
        float alpha = (float)i / total_height;
        float beta = (float)(i - (second_half ? t1.y - t0.y : 0)) / segment_height;
        A.x = t0.x + (t2.x - t0.x) * alpha;
        A.y = t0.y + (t2.y - t0.y) * alpha;

        B.x = second_half ? t1.x + (t2.x - t1.x) * beta : t0.x + (t1.x - t0.x) * beta;
        B.y = second_half ? t1.y + (t2.y - t1.y) * beta : t0.y + (t1.y - t0.y) * beta;
        if (A.x > B.x) std::swap(A, B);
        for (int j = A.x; j <= B.x; j++)
        {
            putpixel_adjusted2(j, t0.y + i, color);
        }
    }
}

void brasterize(vertex& V1, vertex& V2, vertex& V3, float& color)
{
    if (V1.y == V2.y && V2.y == V3.y) return;
    //Bubble sort on y-position
    if (V1.y > V2.y) { swap(V1, V2); }
    if (V1.y > V3.y) { swap(V1, V3); }
    if (V2.y > V3.y) { swap(V3, V2); }

    //divide triangle into two halves

    int height = V3.y - V1.y;
    if (height == 0)return;

    for (int y = V1.y; y <= V2.y; y++)
    {
        int partialHeight = V2.y - V1.y + 1; // +1 because both upper and lower limit is included

        float alpha = (float)(y - V1.y) / height;// be careful with divisions by zero 
        if (partialHeight != 0)
        {
            float beta = (float)(y - V1.y) / partialHeight;
            int Ax = (V1.x + (V3.x - V1.x) * alpha), Ay = V1.y + (V3.y - V1.y) * alpha;
            int Bx = V1.x + (V2.x - V1.x) * beta, By = V1.y + (V2.y - V1.y) * beta;
            if (Ax > Bx) { swap(Ax, Bx); }
            for (int j = Ax; j <= Bx; j++)
                putpixel_adjusted2(j, y, color);
        }

    }

    for (int y = V2.y; y <= V3.y; y++)
    {
        int partialHeight = V3.y - V2.y + 1; // +1 because both upper and lower limit is included

        float alpha = (float)(y - V1.y) / height;
        if (partialHeight != 0)
        {
            float beta = (float)(y - V2.y) / partialHeight; // be careful with divisions by zero 

            int Ax = V1.x + (V3.x - V1.x) * alpha, Ay = V1.y + (V3.y - V1.y) * alpha;
            int Bx = V2.x + (V3.x - V2.x) * beta, By = V2.y + (V3.y - V2.y) * beta;
            if (Ax > Bx) { swap(Ax, Bx); }
            for (int j = Ax; j <= Bx; j++)
                putpixel_adjusted2(j, y, color);
        }

    }
}



vertex Vector_IntersectPlane(vertex& plane_p, vertex& plane_n, vertex& lineStart, vertex& lineEnd)
{
    plane_n = normalize(plane_n);
    float plane_d = -dotProduct(plane_n, plane_p);
    float ad = dotProduct(lineStart, plane_n);
    float bd = dotProduct(lineEnd, plane_n);
    float t = (-plane_d - ad) / (bd - ad);
    vertex lineStartToEnd = lineEnd - lineStart;
    vertex lineToIntersect = lineStartToEnd * t;
    return lineStart + lineToIntersect;
}

int Triangle_ClipAgainstPlane(vertex plane_p, vertex plane_n, triangle& in_tri, triangle& out_tri1, triangle& out_tri2)
{
    // Make sure plane normal is indeed normal
    plane_n = normalize(plane_n);

    // Return signed shortest distance from point to plane, plane normal must be normalised
    auto dist = [&](vertex& p)
    {
        vertex n = normalize(p);
        return (plane_n.x * p.x + plane_n.y * p.y + plane_n.z * p.z - dotProduct(plane_n, plane_p));
    };

    // Create two temporary storage arrays to classify points either side of plane
    // If distance sign is positive, point lies on "inside" of plane
    vertex* inside_points[3];  int nInsidePointCount = 0;
    vertex* outside_points[3]; int nOutsidePointCount = 0;

    // Get signed distance of each point in triangle to plane
    float d0 = dist(in_tri.p[0]);
    float d1 = dist(in_tri.p[1]);
    float d2 = dist(in_tri.p[2]);

    if (d0 >= 0) { inside_points[nInsidePointCount++] = &in_tri.p[0]; }
    else { outside_points[nOutsidePointCount++] = &in_tri.p[0]; }
    if (d1 >= 0) { inside_points[nInsidePointCount++] = &in_tri.p[1]; }
    else { outside_points[nOutsidePointCount++] = &in_tri.p[1]; }
    if (d2 >= 0) { inside_points[nInsidePointCount++] = &in_tri.p[2]; }
    else { outside_points[nOutsidePointCount++] = &in_tri.p[2]; }

    // Now classify triangle points, and break the input triangle into 
    // smaller output triangles if required. There are four possible
    // outcomes...

    if (nInsidePointCount == 0)
    {
        // All points lie on the outside of plane, so clip whole triangle
        // It ceases to exist

        return 0; // No returned triangles are valid
    }

    if (nInsidePointCount == 3)
    {
        // All points lie on the inside of plane, so do nothing
        // and allow the triangle to simply pass through
        out_tri1 = in_tri;

        return 1; // Just the one returned original triangle is valid
    }

    if (nInsidePointCount == 1 && nOutsidePointCount == 2)
    {
        // Triangle should be clipped. As two points lie outside
        // the plane, the triangle simply becomes a smaller triangle

        // Copy appearance info to new triangle
       /* out_tri1.col = in_tri.col;
        out_tri1.sym = in_tri.sym;*/

        // The inside point is valid, so keep that...
        out_tri1.p[0] = *inside_points[0];

        // but the two new points are at the locations where the 
        // original sides of the triangle (lines) intersect with the plane
        out_tri1.p[1] = Vector_IntersectPlane(plane_p, plane_n, *inside_points[0], *outside_points[0]);
        out_tri1.p[2] = Vector_IntersectPlane(plane_p, plane_n, *inside_points[0], *outside_points[1]);

        return 1; // Return the newly formed single triangle
    }

    if (nInsidePointCount == 2 && nOutsidePointCount == 1)
    {
        // Triangle should be clipped. As two points lie inside the plane,
        // the clipped triangle becomes a "quad". Fortunately, we can
        // represent a quad with two new triangles

        // Copy appearance info to new triangles
        /*out_tri1.col = in_tri.col;
        out_tri1.sym = in_tri.sym;

        out_tri2.col = in_tri.col;
        out_tri2.sym = in_tri.sym;*/

        // The first triangle consists of the two inside points and a new
        // point determined by the location where one side of the triangle
        // intersects with the plane
        out_tri1.p[0] = *inside_points[0];
        out_tri1.p[1] = *inside_points[1];
        out_tri1.p[2] = Vector_IntersectPlane(plane_p, plane_n, *inside_points[0], *outside_points[0]);

        // The second triangle is composed of one of he inside points, a
        // new point determined by the intersection of the other side of the 
        // triangle and the plane, and the newly created point above
        out_tri2.p[0] = *inside_points[1];
        out_tri2.p[1] = out_tri1.p[2];
        out_tri2.p[2] = Vector_IntersectPlane(plane_p, plane_n, *inside_points[1], *outside_points[0]);

        return 2; // Return two newly formed triangles which form a quad
    }
}






























void displaymat(mat4d matrix) {
    cout << "display matrix" << endl;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {

            cout << matrix.m[i][j] << "\t";
        }
        cout << endl;

    }
}

vertex matmultvertex(mat4d& m, vertex& i)
{
    vertex v;
    v.x = i.x * m.m[0][0] + i.y * m.m[1][0] + i.z * m.m[2][0] + i.w * m.m[3][0];
    v.y = i.x * m.m[0][1] + i.y * m.m[1][1] + i.z * m.m[2][1] + i.w * m.m[3][1];
    v.z = i.x * m.m[0][2] + i.y * m.m[1][2] + i.z * m.m[2][2] + i.w * m.m[3][2];
    v.w = i.x * m.m[0][3] + i.y * m.m[1][3] + i.z * m.m[2][3] + i.w * m.m[3][3];
      /* v.x = i.x * m.m[0][0] + i.y * m.m[0][1] + i.z * m.m[0][2];
       v.y = i.x * m.m[1][0] + i.y * m.m[1][1] + i.z * m.m[1][2];
       v.z = i.x * m.m[2][0] + i.y * m.m[2][1] + i.z * m.m[2][2];
        v.w = i.x * m.m[0][3] + i.y * m.m[1][3] + i.z * m.m[2][3] + i.w * m.m[3][3];*/

    return v;
}



//void fillTriangle(vertex v1, vertex v2, vertex v3, float& dp)
//{
//    vertex arr[3] = { v1,v2,v3 };
//    //bubbleSort(arr,3);
//   //arr[0].y < arr[1].y ? swap(arr[1], arr[0]) : arr[0].y < arr[2].y ? swap(arr[2], arr[0]) : arr[1].y < arr[2].y ? swap(arr[2], arr[1]) : swap(arr[1], arr[1]);
//    if (arr[0].y > arr[1].y) { swap(arr[1], arr[0]); }
//    if (arr[1].y > arr[2].y) { swap(arr[2], arr[1]); }
//    if (arr[0].y > arr[1].y) { swap(arr[1], arr[0]); }
//
//
//
//    if (int(arr[1].y) == int(arr[2].y))
//    {
//        if (arr[1].x > arr[2].x)swap(arr[2], arr[1]);
//        fillBottomFlatTriangle(arr[0], arr[1], arr[2], dp);
//    }
//    /* check for trivial case of top-flat triangle */
//    else if (int(arr[0].y) == int(arr[1].y))
//    {
//        if (arr[0].x > arr[1].x)swap(arr[0], arr[1]);
//        fillTopFlatTriangle(arr[0], arr[1], arr[2], dp);
//    }
//    else
//    {
//        /* general case - split the triangle in a topflat and bottom-flat one */
//        vertex* v4 = new vertex({
//        (arr[0].x + ((float)(arr[1].y - arr[0].y) / (float)(arr[2].y - arr[0].y)) * (arr[2].x - arr[0].x)), arr[1].y,0 });
//        if (arr[1].x < (*v4).x) {//major left
//            fillBottomFlatTriangle(arr[0], arr[1], *v4, dp);
//            fillTopFlatTriangle(arr[1], *v4, arr[2], dp);
//        }
//        else {//major right
//            fillBottomFlatTriangle(arr[0], *v4, arr[1], dp);
//            fillTopFlatTriangle(*v4, arr[1], arr[2], dp);
//        }
//
//    }


//}

vertex barycentric(vertex A, vertex B, vertex C, vertex P) {
    vertex s[2];
    s[0] = { C.x - A.x , B.x - A.x , A.x - P.x };
    s[1] = { C.y - A.y , B.y - A.y , A.y - P.y };
    vertex u = crossProduct(s[0], s[1]);
    if (std::abs(u.z) > 0) {
        vertex temp;
        temp.x = 1.f - (u.x + u.y) / u.z;
        temp.y = u.y / u.z;
        temp.z = u.x / u.z;
        return temp;
    }
    vertex temp = { -1,1,1 };
    return temp; // in this case generate negative coordinates, it will be thrown away by the rasterizator
}

void bellatriangle(vertex* pts, float* zbuffer, const vec3_T<float>& color, float* intensity)
{

    vertex bboxmin{ 99999999999.9999999999,99999999999.9999999999 ,0 };
    vertex bboxmax{ -99999999999.9999999999, -99999999999.9999999999,0 };

    for (int i = 0; i <3; i++) {

        //kinda redundant but cannot loop this
        bboxmin.x = min(bboxmin.x, pts[i].x);
        bboxmax.x = max(bboxmax.x, pts[i].x);

        bboxmin.y = min(bboxmin.y, pts[i].y);
        bboxmax.y = max(bboxmax.y, pts[i].y);
    }
    vertex P;
    vec3 clr;
    for (int i = bboxmin.x; i <= bboxmax.x; i++)
    {
        for (int j = bboxmin.y; j <= bboxmax.y; j++)
        {
            P.x = static_cast<float>(i);
            P.y = static_cast<float>(j);
            vertex bc_screen = barycentric(pts[0], pts[1], pts[2], P);
            if (bc_screen.x < 0 || bc_screen.y < 0 || bc_screen.z < 0) continue;
            clr = (color * (bc_screen.x) * intensity[0] + color * (bc_screen.y) * intensity[1] + color * (bc_screen.z) * intensity[2]);
            //std::cout << clr << std::endl;
           // clr /= 3;
            P.z = pts[0].z * bc_screen.x + pts[1].z * bc_screen.y + pts[2].z * bc_screen.z + 1;

            putpixel_adjusted3(int(P.x), int(P.y), P.z, clr); //specially designed for working with zbuffer
            // putpixel_adjusted2(P, clr);
          /*  putpixel_adjusted3(P.x + 0.5f, P.y + 0.5f, P.z, clr);
            putpixel_adjusted3(P.x - 0.5f, P.y - 0.5f, P.z, clr);
            putpixel_adjusted3(P.x - 1.0f, P.y - 1.0f, P.z, clr);
            putpixel_adjusted3(P.x + 1.0f, P.y + 1.0f, P.z, clr);*/
           // putpixel_adjusted3(P.x - 0.5f, P.y - 0.5f, P.z, clr);


        }
      
    }
}


vertex vec3tovertex(const vec3& temp) {
    vertex vert;
    vert.x = temp.x;
    vert.y = temp.y;
    vert.z = temp.z;
    return vert;
}

vec3 vertextovec3(const vertex& temp) {
    vec3 vert;
    vert.x = temp.x;
    vert.y = temp.y;
    vert.z = temp.z;
    return vert;
}


mesh getsmesh(ModelParse* object) {
    triangle tris;
    mesh temp;
    for (int i = 0; i < object->nfaces(); i++) {
        std::vector<vec3i> face = object->face(i);
        //std::cout << face[i];
        vec3 points[3];


        for (int j = 0; j < 3; j++)
        {
            //cout << face[j].x;
            vec3 v = object->vert(face[j].x);
            tris.p[j] = vec3tovertex(v);
            tris.p[j].normal = object->norm(i, j);
        }
        temp.triangles.push_back(tris);
    }

    return temp;

}

vertex subvertex(vertex a, vertex b)
{
    vertex temp;
    temp.x = a.x - b.x;
    temp.y = a.y - b.y;
    temp.z = a.z - b.z;
    return temp;
}



//-------------------------------------------------------------------Experimental functions()-----------------------------------------------------------------------------------//
//vec3 rotation3D(vec3 Old, float angleX, float angleY, float angleZ)
//{
//    vec3 New = Old;
//
//    //Y-rotation
//    if (angleY != 0)
//    {
//        New.x = (1.0f * Old.z * sin(angleY) + 1.0f * Old.x * cos(angleY));
//        New.y = Old.y;
//        New.z = (1.0f * Old.z * cos(angleY) - 1.0f * Old.x * sin(angleY));
//        Old = New;
//    }
//
//    //X-rotation
//    if ((angleX) != 0)
//    {
//        New.x = Old.x;
//        New.y = (1.0f * Old.y * cos(angleX) - 1.0f * Old.z * sin(angleX));
//        New.z = (1.0f * Old.y * sin(angleX) + 1.0f * Old.z * cos(angleX));
//        Old = New;
//    }
//
//
//
//    //Zrotation
//    if (angleZ != 0)
//    {
//        New.x = (-1.0f * Old.y * sin(angleZ) + 1.0f * Old.x * cos(angleZ));
//        New.y = (1.0f * Old.y * cos(angleZ) + 1.0f * Old.x * sin(angleZ));
//        New.z = Old.z;
//    }
//    return New;
//}
//mat4 lookat(vec3 eye, vec3 center, vec3 up) //ModelView
//{
//    //eye   :   position of the eye
//    //center:   origin of the new axes
//    //up    :   vertical vector in final render
//
//    //the z-axis is the vector c-e (centre and eye)
//    vec3 z = (eye - center).normalize();
//    //the x-axis is given by cross product between z and up
//    vec3 x = vec3::cross(up, z).normalize();
//    //the y-axis is given by cross product between z and x
//    vec3 y = vec3::cross(x, z).normalize(); //up vector is not necessarily the new y axis . But why??
//
//    mat4 temp;
//    //translate 
//    temp(0, 3) = -center.x;
//    temp(1, 3) = -center.y;
//    temp(2, 3) = -center.z;
//
//    //inverse
//    temp(0, 0) = x.x;
//    temp(1, 0) = y.x;
//    temp(2, 0) = z.x;
//
//    temp(0, 1) = x.y;
//    temp(1, 1) = y.y;
//    temp(2, 1) = z.y;
//
//    temp(0, 2) = x.z;
//    temp(1, 2) = y.z;
//    temp(2, 2) = z.z;
//
//    return temp;
//}
//
//vec3 world2screen(vec3 v) {
//
//    /*
//    mat4 ModelView = lookat(eye, center, vec3(0, 1, 0));
//    mat4 Projection;
//    mat4 ViewPort = viewport(width*0.1,height*0.1, width * 0.8, height * 0.8);
//    Projection(3, 2) = -1.0f / ((eye - center).normalize()).z;
//    return (ViewPort * Projection * ModelView * (v));
//    */
//    //vec3 camera = { 0,0,3 };
//    vec3 light_dir = vec3(0, 1, 1).normalize();;
//    vec3 eye = { 0,0,3 };
//    vec3 center = { 0,0,0 };
//    vec3 up = { 0,0,1 };
//
//    const vec3_T<float>& WHIT = { 1,1,1 };
//    const vec3_T<float>& RED = { 1,0,0 };
//    const vec3_T<float>& BLUE = { 0,1,0 };
//    const vec3_T<float>& GREEN = { 0,0,1 };
//    float scale = 1.0f;
//    float tempx = (float)((int)(((v.x + 1) * window_width / 2. + .5) * scale));
//    float tempy = (float)((int)(((v.y + 1) * window_height / 2. + .5) * scale));
//    float tempz = (float)((v.z * scale));
//
//    //perspective;
//    float r = -1 / eye.z;
//    tempx = (float)((int)(tempx / (1 + r * tempz) + 0.5));
//    tempy = (float)((int)(tempy / (1 + r * tempz) + 0.5));
//    tempz = tempz / (1 + r * tempz);
//
//    lookat(BLUE, GREEN, RED);
//
//    vec3 temp = { tempx,tempy,tempz };
//    return temp;
//
//}

