#pragma once
#pragma once
#include <vector>
#include<fstream>
#include<sstream>
#include<string>
#include"maths.h"

class ModelParse {
private:
    std::vector<vec3> vertices;				//vertex array
    std::vector<std::vector<vec3i>> faces;	//face array
    std::vector<vec2> textures;	//texture co-ordinates array
    std::vector<vec3> normal;	//normal co-ordinates array for Gouraud Shading
public:
    //ModelParse(std::string filename);
    //~ModelParse();
    //int nverts();	//returns number of vertices
    //int nfaces();	//returns number of faces
    //int ntextures();	//returns number of texture vertices
    //int nnormal();		//returns number of normal vertices

    //vec3 vert(int i);	//returns i-th vertices
    //std::vector<vec3i> face(int idx);	//returns idx-th face
    //vec3 norm(int idx, int i);		//returns normal of idx-th face's i-th vertex 

    ModelParse(std::string filename) : vertices(), faces()
    {
        std::cout << "Filename : " << filename;
        //Reads in file
        std::ifstream in;
        in.open(filename, std::ifstream::in);
        //if file cannot be opened
        if (in.fail())
        {
            std::cout << "Cannot Read" << std::endl;
            return;
        }
        std::string line;

        while (!in.eof())
        {
            //get one line at a time
            std::getline(in, line);
            //string object
            std::istringstream iss(line.c_str());

            char trash;
            if (!line.compare(0, 2, "v "))  //starts with v<space>
            {
                iss >> trash; // first character is v
                vec3 v;
                // followed by xyz co-ords
                iss >> v.x;
                iss >> v.y;
                iss >> v.z;
                vertices.push_back(v);
            }
            else if (!line.compare(0, 2, "f ")) //starts with f<space>
            {
                std::vector<vec3i> f;
                vec3i temp;
                iss >> trash; //first charecter is f
                while (iss >> temp.x >> trash >> temp.y >> trash >> temp.z)  // in the form vert/vertTex/norm (vert is read, the rest are treated as trash)
                {
                    //in wavefront obj all indices start at 1, not zero
                    temp.x--;   //vert
                    temp.y--;   //texture
                    temp.z--; // normal 
                    f.push_back(temp);
                }
                faces.push_back(f);

            }
            else if (!line.compare(0, 3, "vt "))    //starts with vt<space>
            {
                iss >> trash >> trash;//Ignore vt
                vec2 uv;
                iss >> uv.x;
                iss >> uv.y;
                textures.push_back(uv);
            }

            else if (!line.compare(0, 3, "vn "))    //starts with vn<space>
            {
                iss >> trash >> trash;//Ignore vn
                vec3 n;
                iss >> n.x;
                iss >> n.y;
                iss >> n.z;
                normal.push_back(n);
            }

        }
        std::cout << std::endl;
        std::cout << nnormal();
    }

    ~ModelParse()
    {
    }

    int nverts()
    {
        return (int)vertices.size();
    }

    int nfaces()
    {
        return (int)faces.size();
    }

    int ntextures()
    {
        return(int)textures.size();
    }

    int nnormal()
    {
        return(int)normal.size();
    }

    vec3 vert(int i)
    {
        return vertices.at(i);
    }

    std::vector<vec3i> face(int idx)
    {
        return faces.at(idx);
    }

    vec3 norm(int idx, int i)
    {
        int temp = faces.at(idx).at(i).z;   //z is normal
        return normal.at(temp).normalize();
        //cout << normal.at(temp).normalize();
    }




};

