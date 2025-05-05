#pragma once

#include "barbs.h"
#include "barbules.h"
#include "feather-masking.h"
#include "feather-base.h"

class Feather : public FeatherBase {
public:
    using Masking = FeatherMasking;

protected:
    Barbules barbules;

public:
    Feather(const Vec3& barb_direction, const Vec3& normal, 
            float barb_excentricity, float barbule_excentricity,
            float barbule_rotation, float barbule_inclination, 
            float barbule_length, //Measured in number of barbs (3 means 3 times the barbs' thickness)
            float barbule_separation //Measured in number of barbules (3 means 3 times the barbules' thickness)
        ) : FeatherBase(barb_direction,normal,barb_excentricity,barbule_excentricity,barbule_rotation,
                                                    barbule_inclination,barbule_length,barbule_separation),
            barbules(barbule_separation,barbule_excentricity)  {}

    Masking masking(const Vec3& dir) const {
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