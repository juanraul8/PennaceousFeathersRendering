#pragma once
#include "flatland.h"

class Barbules {
public:
    class Masking {
        float transmittance_rate_;
        float hmin_;
    public:
        Masking(float transmittance_rate, float hmin = -1.0f) :
            transmittance_rate_(transmittance_rate), hmin_(hmin) {}
        
        constexpr float transmittance_rate() const noexcept { return transmittance_rate_; }
        constexpr float barbule_rate() const noexcept { return 1.0f - transmittance_rate(); }
        constexpr float hmin() const noexcept { return hmin_; }
        constexpr float hmax() const noexcept { return 1.0f; }
    };

private:
    float separation_; //Separation in number of barbules
    float excentricity_; //Excentricity (vertical axis)
public:
    constexpr float separation() const noexcept { return separation_; }
    constexpr float excentricity() const noexcept { return excentricity_; }
    Barbules(float separation = 1.0f, float excentricity = 1.0f) :
        separation_(separation), excentricity_(excentricity) {}

    //Asumes normal at (0,1), projects towards direction dir
    //Templatized because the vector can have a z component that we ignore
    template<typename V2>
    Masking masking(const V2& dir) const {
        //We only need from "upper-left", everything else is symmetric
        flatland::Vec2 d{std::abs(dir[0]),std::abs(dir[1])}; 

        float factor = 1.0f/(2.0f+2.0f*separation());
        flatland::segment segment_separation(factor,0.0f,factor*(1.0f + 2.0f*separation()),0.0f);
        flatland::ellipse barbule1(0.0f,0.0f,factor,factor*excentricity());
        flatland::ellipse barbule2(1.0f,0.0f,factor,factor*excentricity());

        float theta_ini = barbule2.parameter_from_tangent(d) - M_PI;
        float theta_end = theta_ini + M_PI;
        flatland::ray r(barbule1.parametric(theta_end),d);
        std::list<float> ts = barbule2.intersect(r);
        float total_barbule = barbule2.projected_area(d,theta_ini,theta_end);
        if (!ts.empty()) {
            float theta_shadow = barbule2.parameter_at(r.at(ts.front()));
            if (theta_shadow > theta_end) theta_shadow -= 2*M_PI;
            float shadow_h = barbule2.projected_area(d,theta_ini,theta_shadow);          
            float hmin = 2.0f*(shadow_h/total_barbule) - 1.0f;
            //float h = barbule2.h_at(r.at(ts.front()),d);
            //if (std::abs(hmin-h)>1.e-3) std::cerr<<hmin<<" vs. "<<h<<std::endl;
            return Masking(0.0f,hmin);
        } else {
            std::list<float> ts1 = segment_separation.intersect(r);
            flatland::ray r2(barbule2.parametric(barbule2.parameter_from_tangent(d)+M_PI),-d);
            std::list<float> ts2 = segment_separation.intersect(r2);
            float total_separation= (length(r.at(ts1.front())-r2.at(ts2.front())))*d[1];
            //dir[1] = cosine with respect to the normal
            return Masking(total_separation/(total_separation+total_barbule));
        }
    }        
};