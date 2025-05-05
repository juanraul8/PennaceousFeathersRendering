#pragma once
#include "flatland.h"
#include "barb-range.h"
#include "barb-masking.h"
#include <random>

class BarbsMC {
public:
    using Masking = BarbMasking;

private:
    float excentricity_; //Excentricity of barb (vertical axis)
    float barbule_inclination_; //Inclination angle of barbules
    float barbule_length_; // Barbule length in number of barbs (horizontally)
    float left_transparency_; //Transparency rate of left barbule
    float right_transparency_; //Transparency rate of right barbule
    unsigned long repetitions_; //Repetitions of the barbules for Monte-Carloing everything
    unsigned long samples_; //Number of montecarlo samples
public:
    constexpr float barbule_inclination() const noexcept { return barbule_inclination_; }
    constexpr float barbule_length() const noexcept { return barbule_length_; }
    constexpr float excentricity() const noexcept { return excentricity_; }
    constexpr float left_transparency() const noexcept { return left_transparency_; }
    constexpr float right_transparency() const noexcept { return right_transparency_; }
    constexpr unsigned long repetitions() const noexcept { return repetitions_; }
    constexpr unsigned long samples() const noexcept { return samples_; }

    BarbsMC(float excentricity, float barbule_inclination, float barbule_length, 
          float left_transparency = 0.0f, float right_transparency = 0.0f, 
          unsigned long repetitions = 1000, unsigned long samples = 10000) :
        excentricity_(excentricity), barbule_inclination_(barbule_inclination),
        barbule_length_(barbule_length), left_transparency_(left_transparency),
        right_transparency_(right_transparency), 
        repetitions_(repetitions), samples_(samples) {}
    

    //Asumes normal at (0,1), projects towards direction dir
    //Templatized because the vector can have a z component that we ignore
    template<typename V2>
    Masking masking(const V2& d) const {
        flatland::Vec2 dir{d[0],d[1]};
        float factor = 1.0/(2.0f*(1.0 + barbule_length()*std::cos(barbule_inclination())));
        std::list<flatland::ellipse> barbs;
        std::list<flatland::segment> lefts, rights;
        for (unsigned long n = 0; n<repetitions(); ++n) {
            barbs.push_back(flatland::ellipse(factor+float(n),0.0f,factor,factor*excentricity()));
            if (n<(repetitions()-1)) {
                flatland::Vec2 corner{float(n) + factor*(2.0f+barbule_length()*std::cos(barbule_inclination())),factor*barbule_length()*std::sin(barbule_inclination())};
                lefts.push_back(flatland::segment(flatland::Vec2{float(n+1),0.0f},corner));
                rights.push_back(flatland::segment(flatland::Vec2{float(n) + 2.0f*factor,0.0f},corner));
            }
        } 

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

        float hmin = 1.0f, hmax = -1.0f;
        float barb_hits = 0.0f, left_hits = 0.0f, right_hits = 0.0f, transmittance_hits = 0.0f;
        for (unsigned long s = 0; s < samples(); ++s) {
            flatland::ray r(origin.parametric(dist(mt)),dir);
            float throughput = 1.0f; unsigned int bounces = 0;
            while ((throughput > 0.0f) && (bounces < 3)) {
                auto [tbarb,barb] = flatland::intersect(barbs,r);
                auto [tleft,left] = flatland::intersect(lefts,r);
                auto [tright,right] = flatland::intersect(rights,r);
                if ((barb) && ((!left) || (tbarb < tleft)) && ((!right) || (tbarb < tright))) {
                    //Hits the barb first
                    barb_hits += 1.0f*throughput;
                    throughput = 0.0f;
                    float h = barb->h_at(r.at(tbarb),dir);
                    if (h<hmin) hmin = std::max(h,-1.0f);
                    if (h>hmax) hmax = std::min(h,1.0f);
                } else if ((left) && ((!barb) || (tleft < tbarb)) && ((!right) || (tleft < tright)) ) {
                    left_hits += (1.0f - left_transparency())*throughput;
                    throughput *= left_transparency();
                    ++bounces;
                    r = flatland::ray(r.at(tleft + 0.001),dir);
                } else if ((right) && ((!barb) || (tright < tbarb)) && ((!left) || (tright < tleft)) ) {
                    right_hits += (1.0f - right_transparency())*throughput;
                    throughput *= right_transparency();
                    ++bounces;
                    r = flatland::ray(r.at(tright + 0.001),dir);
                } else {  //It hits nothing
                    transmittance_hits += 1.0f*throughput;
                    throughput = 0.0f;
                }
            }
        }
        std::list<BarbRange> barb_ranges;
        float barb_rate = float(barb_hits)/float(samples());
        float left_rate = float(left_hits)/float(samples());
        float right_rate = float(right_hits)/float(samples());
        float transmittance_rate = float(transmittance_hits)/float(samples());
        barb_ranges.push_back(BarbRange(hmin,hmax,barb_rate));
        return Masking(transmittance_rate,left_rate,right_rate,barb_ranges);
    }        
};
