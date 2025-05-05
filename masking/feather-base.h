#pragma once

#include "math-array.h"
using Vec3 = std::array<float,3>;
using Mat3x3 = std::array<Vec3,3>; //Vector of rows

class FeatherBase {
protected:
//    Barbules barbules;
    float barb_excentricity_, barbule_inclination_, barbule_length_;
    Mat3x3 global_to_barbs;
    Mat3x3 global_to_left;
    Mat3x3 global_to_right;

public:
    FeatherBase(const Vec3& barb_direction, const Vec3& normal, 
            float barb_excentricity, float barbule_excentricity,
            float barbule_rotation, float barbule_inclination, 
            float barbule_length, //Measured in number of barbs (3 means 3 times the barbs' thickness)
            float barbule_separation //Measured in number of barbules (3 means 3 times the barbules' thickness)
        ) : // barbules(barbule_separation,barbule_excentricity),
        barb_excentricity_(barb_excentricity), barbule_inclination_(barbule_inclination),
        barbule_length_(barbule_length) {

            Mat3x3 barbs_to_global = matrix(
                normalized(cross(normal,barb_direction)),
                normalized(normal),
                normalized(barb_direction)
            );

            Vec3 local_left_z{
                std::cos(barbule_inclination)*std::sin(-barbule_rotation),
                std::sin(barbule_inclination),
                std::cos(barbule_inclination)*std::cos(-barbule_rotation)};
            Vec3 local_left_x = normalized(cross(Vec3{0.0f,1.0f,0.0f},local_left_z));

            Vec3 local_right_z{
                std::cos(barbule_inclination)*std::sin(barbule_rotation),
                std::sin(barbule_inclination),
                std::cos(barbule_inclination)*std::cos(barbule_rotation)};
            Vec3 local_right_x = normalized(cross(Vec3{0.0f,1.0f,0.0f},local_right_z));

            Mat3x3 left_to_barbs = matrix(local_left_x,cross(local_left_z,local_left_x),local_left_z);
            Mat3x3 right_to_barbs = matrix(local_right_x,cross(local_right_z,local_right_x),local_right_z);

            global_to_barbs = inverse(barbs_to_global);
            global_to_left = inverse(matrix_product(barbs_to_global,left_to_barbs));
            global_to_right = inverse(matrix_product(barbs_to_global,right_to_barbs));
/*        std::cerr<<"  BARBS TO GLOBAL "<<barbs_to_global<<std::endl;    
        std::cerr<<"  LEFT TO BARBS "<<left_to_barbs<<std::endl;    
        std::cerr<<"  RIGHT TO BARBS "<<right_to_barbs<<std::endl;    
        std::cerr<<"  GLOBAL TO BARBS "<<global_to_barbs<<std::endl;
        std::cerr<<"  GLOBAL TO LEFT "<<global_to_left<<std::endl;
        std::cerr<<"  GLOBAL TO RIGHT "<<global_to_right<<std::endl; */
    }
};