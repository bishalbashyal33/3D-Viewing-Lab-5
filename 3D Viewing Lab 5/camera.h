#pragma once
#pragma once
#include"graphics.h"
#include"maths.h"
#include"matrices.h"
#include"mathtools.h"
#include"structures.h"
#include"modelparser.h"


float clamp(float x, float upper, float lower)
{
    return min(upper, max(x, lower));
}

//// Defines several possible options for camera movement. Used as abstraction to stay away from window-system specific input methods
//enum Camera_Movement {
//    FORWARD,
//    BACKWARD,
//    LEFT,
//    RIGHT
//};
//
//const float YAW = -90.0f;
//const float PITCH = 0.0f;
//const float SPEED = 2.5f;
//const float SENSITIVITY = 0.1f;
//const float ZOOM = 45.0f;
//
//class Camera
//{
//public:
//    // camera Attributes
//    vertex vCamera;
//    vertex Position;
//    vertex Front;
//    vertex Up;
//    vertex Right;
//    vertex WorldUp;
//    // euler Angles
//    float Yaw;
//    float Pitch;
//    // camera options
//    float MovementSpeed;
//    float MouseSensitivity;
//    float Zoom;
//
//
//
//    Camera(vertex position = { 0.0f, 0.0f, 0.0f }, vertex up = { 0.0f, 1.0f, 0.0f }, float yaw = YAW, float pitch = PITCH) : Front(vertex{ 0.0f, 0.0f, -1.0f }), MovementSpeed(SPEED), MouseSensitivity(SENSITIVITY), Zoom(ZOOM)
//    {
//        Position = position;
//        WorldUp = up;
//        Yaw = yaw;
//        Pitch = pitch;
//        updateCameraVectors();
//    }
//    // constructor with scalar values
//    Camera(float posX, float posY, float posZ, float upX, float upY, float upZ, float yaw, float pitch) : Front(vertex{ 0.0f, 0.0f, -1.0f }), MovementSpeed(SPEED), MouseSensitivity(SENSITIVITY), Zoom(ZOOM)
//    {
//        Position = { posX, posY, posZ };
//        WorldUp = { upX, upY, upZ };
//        Yaw = yaw;
//        Pitch = pitch;
//        updateCameraVectors();
//    }
//
//    void ProcessKeyboard(Camera_Movement direction, float deltaTime)
//    {
//        float velocity = MovementSpeed * deltaTime;
//        if (direction == FORWARD)
//            Position.x += Front.x * velocity;
//            Position.y += Front.y * velocity;
//            Position.z += Front.z * velocity;
//            
//        if (direction == BACKWARD)
//        Position.x -= Front.x * velocity;
//        Position.y -= Front.y * velocity;
//        Position.z -= Front.z * velocity;
//        if (direction == LEFT)
//
//        Position.x -= Right.x * velocity;
//        Position.y -= Right.y * velocity;
//        Position.z -= Right.z * velocity;
//            //Position -= Right * velocity;
//        if (direction == RIGHT)
//           // Position += Right * velocity;
//        Position.x += Right.x * velocity;
//        Position.y += Right.y * velocity;
//        Position.z += Right.z * velocity;
//    }
//
//    void ProcessMouseMovement(float xoffset, float yoffset, GLboolean constrainPitch = true)
//    {
//        xoffset *= MouseSensitivity;
//        yoffset *= MouseSensitivity;
//
//        Yaw += xoffset;
//        Pitch += yoffset;
//
//        // make sure that when pitch is out of bounds, screen doesn't get flipped
//        if (constrainPitch)
//        {
//            if (Pitch > 89.0f)
//                Pitch = 89.0f;
//            if (Pitch < -89.0f)
//                Pitch = -89.0f;
//        }
//
//        // update Front, Right and Up Vectors using the updated Euler angles
//        updateCameraVectors();
//    }
//    void ProcessMouseScroll(float yoffset)
//    {
//        Zoom -= (float)yoffset;
//        if (Zoom < 1.0f)
//            Zoom = 1.0f;
//        if (Zoom > 45.0f)
//            Zoom = 45.0f;
//    }
//
//private:
//    // calculates the front vector from the Camera's (updated) Euler Angles
//    void updateCameraVectors()
//    {
//        // calculate the new Front vector
//        vertex front;
//        front.x = cosf((Yaw)) * cos((Pitch));
//        front.y = sinf((Pitch));
//        front.z = sinf((Yaw)) * cosf(Pitch);
//        Front = normalize(front);
//        // also re-calculate the Right and Up vector
//        Right = normalize(crossProduct(Front, WorldUp));  // normalize the vectors, because their length gets closer to 0 the more you look up or down which results in slower movement.
//        Up = normalize(crossProduct(Right, Front));
//    }
//};

struct cam {


    vertex vtarget = { 0,0,1 };;//z  axis
    vertex upvector = { 0,1,0 };
    float yaw = 0.0f;
    vertex vCamera = { 0,0,-5.0f };

    vertex target() {
        vertex temp = lookdir();
        return vCamera + temp;
    }

    vertex lookdir() {
        mat4d temp = matroty();
        return matmultvertex(temp, vtarget);
      
    }

    mat4d matcam() {
        vertex targ = target();
        return  pointAtMatrix(vCamera, targ, upvector);;
    }

    mat4d matroty() {
        mat4d temp = RotationMatrixY(yaw);
        /*  displaymat(temp);*/
        return temp;
    }

    mat4d matview() {
     //   cout << "Eye Vector:" << vCamera;
        mat4d temp = matcam();
        return Matrix_Inverse(temp);
    }


};


struct ptlight {
    float ambientLight = 0.4f;
    float ambientcoefficient = 1.0f;
    float diffusecoefficient = 0.6f;
    vertex pointlight = { 0.0f,0.0f,-10.0f };
    vertex pointToLightvector[3];
    vertex pointtoCameraVector[3];
    float ambientoffset = 0.9f;
    float intensity[3];
    float specularcoefficient = 20.0f;;
    float specularalpha = 120.0f;


    //our  code does not rotate the normal, so, the lgiht is rendered based on the  initial normal positioning, this  needs to be reworked.
    void caclulateintensity(triangle tri, triangle triTranslated) {
        //  Normal Rotation is Remaining
       // cout << "Lighting Vector:" << pointlight;
        for (int i = 0; i < 3; i++) {
            pointToLightvector[i] = normalize(pointlight - triTranslated.p[i]); //vector from light to point

            intensity[i] = dotProduct(vec3tovertex(tri.p[i].normal), pointToLightvector[i]);

            if (intensity[i] < 0.0f) {
                intensity[i] = 0.0f; //for ambitnoffset setting intsity to 0 when intesnity is negative

            }
            //intensity[j] = intensity[j]*0.8f;
            intensity[i] += 0.2f;
            //intensity[i] = 0.2f;
        }
        /*  intensity[1] = 0.2f;
          intensity[2] = 0.2f;*/
    }

  

    // bug found: it is copying the intensity value of last triangle for every triangle face.


    void calculatephongintensity(triangle tri, triangle tritranslated,cam vcam,mat4d modalmatrix) {
       
      
       

        for (int i = 0; i < 3; i++) {
            vertex temp = vec3tovertex(tri.p[i].normal);
          //  tri.p[i].normal = vertextovec3(matmultvertex(modalmatrix, temp));

            pointToLightvector[i] = normalize(pointlight - tritranslated.p[i]);
            float diffuselighting = clamp(diffusecoefficient * (max(dotProduct(vec3tovertex(tri.p[i].normal), pointToLightvector[i]),0)),1.0f,0.0f);
            vertex reflectionvector = normalize(reflect(pointToLightvector[i]*-1, vec3tovertex(tritranslated.p[i].normal)));
            pointtoCameraVector[i]= normalize(vcam.vCamera - tritranslated.p[i]);
           
            float ambientlighting = ambientLight * ambientcoefficient;
            float specularlight = clamp(specularcoefficient * pow(max(dotProduct(reflectionvector, pointtoCameraVector[i]),0.0f),specularalpha),1.0f,0.0f);
           // cout << specularlight<<endl;
            intensity[i] =clamp(specularlight+ambientlighting+diffuselighting,1,0);
           cout << "Intesnity "<<intensity[i]<<endl;

        }

    }

    vertex reflect(vertex incidentlight, vertex normal) {
      
        return subvertex(incidentlight,normal*(2 * dotProduct(incidentlight, normal)));
    }

    //this code segment need debugging
    void calculatefaceintensity(triangle& in, triangle& triTranslated,mat4d& matrotz) {
        std::vector<vec3i> face;
        vertex vert[3];
        triangle t1;
        int i;
        for (int k = 0; k < 3; k++)
                    {
            vertex t1 = vec3tovertex(in.p[k].normal);
           // vertex t1 = matmultvertex(matrotz,temp);
           // cout << v;
                            //  triTranslated.p[j] = vec3tovertex(v);
                            pointToLightvector[k] = normalize(pointlight-triTranslated.p[k]); //vector from light to point

                          //  vertex t1 = vec3tovertex(v);
                            /*   vertex temp = normalize(matmultvertex(rotz, t1));*/
                            intensity[k] = -dotProduct(t1, pointToLightvector[k]);
                            if (intensity[k] < 0.0f) {
                                intensity[k] = 0.0f; //for ambitnoffset setting intsity to 0 when intesnity is negative

                            }
                            //intensity[j] = intensity[j]*0.8f;
                            intensity[k] += 0.2f;
                            // cout << "intensity" << intensity[j]<<endl;
                          

                    }
       

    }

    void getintensity(triangle& triProjected) {
        triProjected.p[0].intensity = intensity[0];

        triProjected.p[1].intensity = intensity[1];

        triProjected.p[2].intensity = intensity[2];

    }



};