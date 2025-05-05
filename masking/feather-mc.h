#pragma once

#include "barbs-mc.h"
#include "barbules-mc.h"
#include "feather-masking.h"
#include "feather-base.h"

class FeatherMC : public FeatherBase {
public:
    using Masking = FeatherMasking;

protected:
    BarbulesMC barbules;
    unsigned long nsamples_;

public:
    FeatherMC(const Vec3& barb_direction, const Vec3& normal, 
            float barb_excentricity, float barbule_excentricity,
            float barbule_rotation, float barbule_inclination, 
            float barbule_length, //Measured in number of barbs (3 means 3 times the barbs' thickness)
            float barbule_separation, //Measured in number of barbules (3 means 3 times the barbules' thickness)
            unsigned long nsamples
        ) : FeatherBase(barb_direction,normal,barb_excentricity,barbule_excentricity,barbule_rotation,
                                                    barbule_inclination,barbule_length,barbule_separation),
            barbules(barbule_separation,barbule_excentricity,1000,1), nsamples_(nsamples)  {}

    Masking masking(const Vec3& dir) const {
        auto leftdir = product(global_to_left,dir);
        auto rightdir = product(global_to_right,dir);
        auto barbdir = product(global_to_barbs,dir); 

        float transmittance_rate = 0.0f;
        float left_rate = 0.0f;
        float right_rate=0.0f;
        float barb_rate=0.0f;
        for (unsigned long i = 0; i<nsamples_; ++i) { 
            auto left_masking = barbules.masking(leftdir);
            auto right_masking = barbules.masking(rightdir);
            auto barb_masking =
                BarbsMC(barb_excentricity_,barbule_inclination_,barbule_length_,
                    left_masking.transmittance_rate(), right_masking.transmittance_rate(), 1000,1).
                        masking(barbdir);
            transmittance_rate += barb_masking.transmittance_rate();
            left_rate += barb_masking.left_rate();
            right_rate += barb_masking.right_rate();
            for (const BarbRange& range : barb_masking.barb_ranges())
                barb_rate += range.rate();
        } 
        return Masking(
            transmittance_rate/float(nsamples_),
            BarbRange(-1,1,left_rate/float(nsamples_)),
            BarbRange(-1,1,right_rate/float(nsamples_)),
            std::list<BarbRange>{BarbRange(-1,1,barb_rate/float(nsamples_))}
        );
    }
};