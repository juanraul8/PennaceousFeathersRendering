#pragma once
#include "flatland.h"
#include "barbule-masking.h"
#include <random>

class BarbulesMC {
public:
    using Masking = BarbuleMasking;

private:
    float separation_; //Separation in number of barbules
    float excentricity_; //Excentricity (vertical axis)
    float repetitions_; //Repetitions of the barbules for Monte-Carloing everything
    unsigned long samples_; //Number of montecarlo samples
public:
    constexpr float separation() const noexcept { return separation_; }
    constexpr float excentricity() const noexcept { return excentricity_; }
    constexpr float repetitions() const noexcept { return repetitions_; }
    constexpr unsigned long samples() const noexcept { return samples_; }
    BarbulesMC(float separation = 1.0f, float excentricity = 1.0f, unsigned long repetitions = 1000, unsigned long samples = 10000) :
        separation_(separation), excentricity_(excentricity), repetitions_(repetitions), samples_(samples) {}

    //Asumes normal at (0,1), projects towards direction dir
    //Templatized because the vector can have a z component that we ignore
    template<typename V2>
    Masking masking(const V2& d) const {
        flatland::Vec2 dir{d[0],d[1]};
        float factor = 1.0f/(2.0f+2.0f*separation());
        std::list<flatland::ellipse> barbules;
        for (unsigned long n = 0; n<repetitions(); ++n) 
           barbules.push_back(flatland::ellipse(factor+float(n),0.0f,factor,factor*excentricity())); 
        
        flatland::Vec2 o0,o1;
        float distance = float(repetitions())*std::abs(dir[1]);
        if (dir[0]<0.0) {
            o0 = flatland::Vec2{float(repetitions()),0.0f}-10.0f*dir;
            if (dir[1]<0.0) o1 = o0 + distance*flatland::Vec2{dir[1],-dir[0]};
            else            o1 = o0 - distance*flatland::Vec2{dir[1],-dir[0]};     
        } else {
            o1 = -10.0f*dir;
            if (dir[1]<0.0) o0 = o1 - distance*flatland::Vec2{dir[1],-dir[0]};
            else            o0 = o1 + distance*flatland::Vec2{dir[1],-dir[0]};
        } 
        flatland::segment origin(o0,o1);
        std::random_device rd;
        std::mt19937 mt(rd());
        std::uniform_real_distribution<float> dist(-1.0f, 1.0f);

        float hmin = 1.0f;
        float hmax = -1.0f;
        unsigned long barbule_hits = 0;
        for (unsigned long s = 0; s < samples(); ++s) {
            flatland::ray r(origin.parametric(dist(mt)),dir);
            auto [t,barbule] = flatland::intersect(barbules,r);
            if (barbule) {
                ++barbule_hits;
                float h = barbule->h_at(r.at(t),dir);
//                float theta = barbule->parameter_at(r.at(t));
//                float theta_min = barbule->parameter_from_tangent(dir) - M_PI;
//                float theta_max = theta_min + M_PI;
//                float total_theta = barbule->projected_area(dir,theta_min,theta_max);
//                std::cerr<<h<<" vs. "<<(2.0f*(barbule->projected_area(dir,theta_min, theta)/total_theta)-1.0f)<<std::endl;
                if (h<hmin) hmin = std::max(h,-1.0f);
                if (h>hmax) hmax = std::min(h,1.0f);
            }
        }
//        std::cerr<<"["<<hmin<<","<<hmax<<"]"<<std::endl;
        return Masking(1.0f - float(barbule_hits)/float(samples()), hmin, hmax);
    }        
};
