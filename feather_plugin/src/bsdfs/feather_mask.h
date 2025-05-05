#pragma once
#include "math-array.h"
#include "barbs.h"
#include "barbules.h"

using Vec3 = std::array<float,3>;
using Mat3x3 = std::array<Vec3,3>; //Vector of rows

class Feather_masking {
public:
    class Masking {
        float transmittance_rate_;
        BarbRange left_range_, right_range_;
        std::list<BarbRange> barb_ranges_;
    public:
        Masking(float transmittance_rate, const BarbRange& left_range, const BarbRange& right_range, 
            const std::list<BarbRange>& barb_ranges) :
            transmittance_rate_(transmittance_rate), left_range_(left_range),
            right_range_(right_range),
            barb_ranges_(barb_ranges) {}
        
        constexpr float transmittance_rate() const noexcept { return transmittance_rate_; }
        constexpr const BarbRange& left_range() const noexcept { return left_range_; }
        constexpr const BarbRange& right_range() const noexcept { return right_range_; }
        const std::list<BarbRange>& barb_ranges() const noexcept { return barb_ranges_; }
    };
protected:
    Barbules barbules;
    float barb_excentricity_, barbule_inclination_, barbule_length_;
    Mat3x3 global_to_barbs;
    Mat3x3 global_to_left;
    Mat3x3 global_to_right;

public:
    Feather_masking(const Vec3& barb_direction, const Vec3& normal, 
            float barb_excentricity, float barbule_excentricity,
            float barbule_rotation, float barbule_inclination, 
            float barbule_length, //Measured in number of barbs (3 means 3 times the barbs' thickness)
            float barbule_separation //Measured in number of barbules (3 means 3 times the barbules' thickness)
        ) : barbules(barbule_separation,barbule_excentricity),
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


    Vec3 getLeftDir(const Vec3& dir){
        return product(global_to_left,dir);
    }
    Vec3 getRightDir(const Vec3& dir){
        return product(global_to_right,dir);
    }
    Vec3 getBarbDir(const Vec3& dir){
        return product(global_to_barbs,dir); 
    }

    Masking setDirs(const Vec3& barbdir,const Vec3& leftdir,const Vec3& rightdir){
        auto left_masking = barbules.masking(leftdir);
        auto right_masking = barbules.masking(rightdir);
        auto barb_masking =
            Barbs(barb_excentricity_,barbule_inclination_,barbule_length_,
                left_masking.transmittance_rate(), right_masking.transmittance_rate()).
                    masking(barbdir);
        return Masking(
            barb_masking.transmittance_rate(),
            BarbRange(left_masking.hmin(),left_masking.hmax(),barb_masking.left_rate()),
            BarbRange(right_masking.hmin(),right_masking.hmax(),barb_masking.right_rate()),
            barb_masking.barb_ranges()
        );
         
    }

    Masking masking(const Vec3& dir) {
        auto leftdir = product(global_to_left,dir);
        auto rightdir = product(global_to_right,dir);
        auto barbdir = product(global_to_barbs,dir); 

        auto left_masking = barbules.masking(leftdir);
        auto right_masking = barbules.masking(rightdir);
        auto barb_masking =
            Barbs(barb_excentricity_,barbule_inclination_,barbule_length_,
                left_masking.transmittance_rate(), right_masking.transmittance_rate()).
                    masking(barbdir);
        return Masking(
            barb_masking.transmittance_rate(),
            BarbRange(left_masking.hmin(),left_masking.hmax(),barb_masking.left_rate()),
            BarbRange(right_masking.hmin(),right_masking.hmax(),barb_masking.right_rate()),
            barb_masking.barb_ranges()
        );
    }
};