#pragma once
#include "flatland.h"


class BarbRange {
    float hmin_, hmax_, rate_;
public:
    BarbRange(float hmin, float hmax, float rate) :
        hmin_(hmin), hmax_(hmax), rate_(rate) {}

    constexpr float hmin() const noexcept { return hmin_; }
    constexpr float hmax() const noexcept { return hmax_; }
    
    //This is a temporary patch, we should fix the H order issue directly on our masking logic
    /*
    constexpr float hmin() const noexcept { 
        
        if (hmax_< hmin_){

            printf("H order is wrong: [%f, %f]\n", hmin_, hmax_);

            return hmax_; 
        }

        return hmin_; 
    }

    constexpr float hmax() const noexcept { 
        
        if (hmax_ < hmin_){

            printf("H order is wrong: [%f, %f]\n", hmin_, hmax_);

            return hmin_; 
        }

        return hmax_; 
    }
    */
    
    constexpr float rate() const noexcept { return rate_; }
    BarbRange& operator*=(float f) {
        rate_*=f;
        return (*this);
    }
    BarbRange& operator/=(float f) {
        rate_/=f;
        return (*this);
    }
};

class Barbs {
public:
    class Masking {
        float transmittance_rate_, left_rate_, right_rate_;
        std::list<BarbRange> barb_ranges_;
    public:
        Masking(float transmittance_rate, float left_rate, float right_rate, 
            const std::list<BarbRange>& barb_ranges) :
            transmittance_rate_(transmittance_rate), left_rate_(left_rate),
            right_rate_(right_rate),
            barb_ranges_(barb_ranges) {}
        
        constexpr float transmittance_rate() const noexcept { return transmittance_rate_; }
        constexpr float left_rate() const noexcept { return left_rate_; }
        constexpr float right_rate() const noexcept { return right_rate_; }
        const std::list<BarbRange>& barb_ranges() const noexcept { return barb_ranges_; }
    };

private:
    float excentricity_; //Excentricity of barb (vertical axis)
    float barbule_inclination_; //Inclination angle of barbules
    float barbule_length_; // Barbule length in number of barbs (horizontally)
    float left_transparency_; //Transparency rate of left barbule
    float right_transparency_; //Transparency rate of right barbule
public:
    constexpr float barbule_inclination() const noexcept { return barbule_inclination_; }
    constexpr float barbule_length() const noexcept { return barbule_length_; }
    constexpr float excentricity() const noexcept { return excentricity_; }
    constexpr float left_transparency() const noexcept { return left_transparency_; }
    constexpr float right_transparency() const noexcept { return right_transparency_; }

    Barbs(float excentricity, float barbule_inclination, float barbule_length, 
          float left_transparency = 0.0f, float right_transparency = 0.0f) :
        excentricity_(excentricity), barbule_inclination_(barbule_inclination),
        barbule_length_(barbule_length), left_transparency_(left_transparency),
        right_transparency_(right_transparency) {}
    

private:
    //Factor should include transparency
    void push_barb_range(std::list<BarbRange>& barb_ranges, const flatland::ellipse& barb, const flatland::Vec2& d, float theta_min, float theta_max, float theta_ini, float theta_end, float factor) const {
        if ((factor <= 0.0f) || (theta_min >= theta_max)) return;
        float total_barb = barb.projected_area(d,theta_ini,theta_end);
        float h_min = 1.0f - 2.0f*barb.projected_area(d,theta_min,theta_end)/total_barb;
        float h_max = 1.0f - 2.0f*barb.projected_area(d,theta_max,theta_end)/total_barb;        
        float rate = factor*barb.projected_area(d,theta_min, theta_max);
        barb_ranges.push_back(BarbRange(h_min,h_max,rate));
    } 

public:
    //Asumes normal at (0,1), projects towards direction dir
    //Templatized because the vector can have a z component that we ignore
    template<typename V2>
    Masking masking(const V2& dir) const {
        //We only need from "upper-left", everything else is symmetric
        flatland::Vec2 d{dir[0],dir[1]};
        float angle = barbule_inclination();
        if (d[1]<0.0f) {
            d[1] *= -1.0f;
            angle *= -1.0f;
        }

        //We are going to need to flip them
        float left_transp = left_transparency(); float right_transp = right_transparency();
        bool flip = false; //Flip left from right because you've flipped the input direction
        if (d[0]<0.0f) {
            d[0] *= -1.0f;
            flip = true;
            std::swap(left_transp,right_transp);
        }

        float factor = 1.0/(2.0f*(1.0 + barbule_length()*std::cos(angle)));
        flatland::ellipse barb1(factor,0.0f,factor,factor*excentricity());
        flatland::ellipse barb2(1.0f+factor,0.0f,factor,factor*excentricity());
        flatland::Vec2 corner{factor*(2.0f+barbule_length()*std::cos(angle)),factor*barbule_length()*std::sin(angle)};
        flatland::segment right1(flatland::Vec2{2.0f*factor,0.0f},corner);
        flatland::segment left1(flatland::Vec2{1.0f,0.0f},corner);
        flatland::segment right2(flatland::Vec2{1.0f + 2.0f*factor,0.0f},corner + flatland::Vec2{1.0f,0.0f});
        flatland::segment left2(flatland::Vec2{2.0f,0.0f},corner + flatland::Vec2{1.0f,0.0f});
 
        float theta_ini = barb2.parameter_from_tangent(d) - M_PI;
        float theta_end = theta_ini + M_PI;
        
        flatland::ray rbarb(barb1.parametric(theta_end),d);
        flatland::ray rbarbule(corner,d);
        
        std::list<BarbRange> barb_ranges;
        float theta_hard_shadow = theta_ini; // We will account for barbule transmittance
        float theta_left_shadow = -M_PI; // It is always here
        float theta_right_shadow = theta_ini; // We will find it if there is one

        std::list<float> ts = barb2.intersect(rbarb); //Barbs cast full shadow
        if (!ts.empty()) {
            float theta_hit = barb2.parameter_at(rbarb.at(ts.front()));
            if (theta_hit > theta_end) theta_hit -=2.0*M_PI;
            theta_hard_shadow = std::max(theta_hard_shadow,theta_hit);
        }
        ts = barb2.intersect(rbarbule);
        if (!ts.empty()) {
            float theta_hit = barb2.parameter_at(rbarbule.at(ts.front()));
            if (theta_hit > theta_end) theta_hit -=2.0*M_PI;
            theta_right_shadow = std::max(theta_right_shadow,theta_hit);
        } else {
            auto top = barb2.parametric(-0.5*M_PI);
            if (corner[1]<top[1]) { 
                //Origin of ray highher than barb top, can be fully ocludded
                auto other_point = top;
                other_point[1] += 2.0f*(corner[1]-top[1]);
                //Straight segment checks for full oclussion
                if (!flatland::segment(top,other_point).intersect(rbarbule).empty()) {
                    theta_right_shadow = theta_end;
                }
            }
        }

        if ((theta_hard_shadow < theta_left_shadow) && (theta_hard_shadow < theta_right_shadow)) {
            if (theta_left_shadow < theta_right_shadow) {
                push_barb_range(barb_ranges,barb2,d,theta_right_shadow,theta_end,theta_ini,theta_end,1.0f);
                push_barb_range(barb_ranges,barb2,d,theta_left_shadow,theta_right_shadow,theta_ini,theta_end,left_transp*right_transp);
                push_barb_range(barb_ranges,barb2,d,theta_hard_shadow,theta_left_shadow,theta_ini,theta_end,right_transp);
            } else {
                push_barb_range(barb_ranges,barb2,d,theta_left_shadow,theta_end,theta_ini,theta_end,1.0f);
                push_barb_range(barb_ranges,barb2,d,theta_right_shadow,theta_left_shadow,theta_ini,theta_end,left_transp);
                push_barb_range(barb_ranges,barb2,d,theta_hard_shadow,theta_right_shadow,theta_ini,theta_end,right_transp);
            }
        } else if (theta_hard_shadow < theta_left_shadow) {
                push_barb_range(barb_ranges,barb2,d,theta_left_shadow,theta_end,theta_ini,theta_end,1.0f);
                push_barb_range(barb_ranges,barb2,d,theta_hard_shadow,theta_left_shadow,theta_ini,theta_end,left_transp);
        } else if (theta_hard_shadow < theta_right_shadow) {
                push_barb_range(barb_ranges,barb2,d,theta_right_shadow,theta_end,theta_ini,theta_end,1.0f);
                push_barb_range(barb_ranges,barb2,d,theta_hard_shadow,theta_right_shadow,theta_ini,theta_end,left_transp*right_transp);
        } else {
                push_barb_range(barb_ranges,barb2,d,theta_hard_shadow,theta_end,theta_ini,theta_end,1.0f);
        }

        float barb_rate = 0.0f;
        for (const BarbRange& range : barb_ranges) barb_rate += range.rate();

//        std::cerr<<"BARB RATE = "<<barb_rate<<" - TOTAL BARB = "<<total_barb<<std::endl;
 
        float t_right_min = -1.0f; float t_right_selfshadow_min = -1.0f;
        ts = right1.intersect_unbounded(rbarb);
        if ((ts.empty()) || (dot(d,right1.normal())<0)) t_right_min = 1.0f;  //Fully occludded by itself 
        else t_right_min = std::min(1.0f,std::max(t_right_min,right1.parameter_at(rbarb.at(ts.front()))));
        ts = right2.intersect_unbounded(rbarbule);
        if (ts.empty()) t_right_selfshadow_min = 1.0f; //Fully occluded by itself
        else t_right_selfshadow_min  = std::min(1.0f,std::max(t_right_selfshadow_min ,right2.parameter_at(rbarbule.at(ts.front()))));
    
        if (t_right_selfshadow_min < t_right_min) t_right_selfshadow_min = t_right_min; 

        //We add the direct light (acknowledging potential transparencies)
        float right_rate = (1.0f - right_transp)*right1.projected_area(d,t_right_selfshadow_min,1.0f);
        //We account for self_shadowing
        right_rate += left_transp*right_transp*right1.projected_area(d, t_right_min, t_right_selfshadow_min);

    
        float left_rate = 0.0f; float t_left_transparency_max = 1.0f;
        //Two options for left: illuminated from below or from above.
        if (dot(d,left1.normal())<0) { //From above
            float t_left_max = 1.0f; 
            ts = left1.intersect_unbounded(rbarb);
            if (ts.empty()) t_left_max = -1.0f;
            else t_left_max = std::max(-1.0f,std::min(t_left_max,left1.parameter_at(rbarb.at(ts.front()))));

            ts = left2.intersect_unbounded(rbarbule);
            if (ts.empty()) t_left_transparency_max = -1.0f;
            else t_left_transparency_max = std::max(-1.0f,std::min(t_left_transparency_max,left2.parameter_at(rbarbule.at(ts.front()))));
            
            if (t_left_transparency_max > t_left_max) t_left_transparency_max = t_left_max;            
            //We add whatever is visible directly
            left_rate += (1.0f-left_transp)*left1.projected_area(d,-1.0f,t_left_transparency_max);
            //We add after a two barbules transparency
            left_rate += left_transp*right_transp*left1.projected_area(d,t_left_transparency_max,t_left_max);
        } else {
            t_left_transparency_max = -1.0f; //No transparency
            //From below
            float t_left_min = -1.0f; float t_left_transparency_min = -1.0f;
            ts = left1.intersect_unbounded(rbarb);
            if (ts.empty()) t_left_min = 1.0f;
            else t_left_min = std::min(1.0f,std::max(t_left_min,left1.parameter_at(rbarb.at(ts.front()))));
            //The same for barbules, single transparency or triple transparency

            ts = left2.intersect_unbounded(rbarbule);
            if (ts.empty()) t_left_transparency_min = 1.0f;
            else t_left_transparency_min = std::min(1.0f,std::max(t_left_transparency_min,left2.parameter_at(rbarbule.at(ts.front()))));
            
            if (t_left_transparency_min < t_left_min) t_left_transparency_min = t_left_min;            
            //We add whatever is visible directly through the right barbule
            left_rate += right_transp*(1.0f-left_transp)*left1.projected_area(d,t_left_transparency_min,1.0f);
            //We add after a two barbules transparency
            left_rate += right_transp*left_transp*right_transp*left1.projected_area(d,t_left_min,t_left_transparency_min);
        }

        flatland::ray rback(barb2.parametric(theta_ini),-d);
        float transparency_rate = 0.0f;
        //t_right_transparency_min = t_right_selfshadow_min
        float t_right_transparency_max = std::min(1.0f,right1.parameter_at(rback.at(right1.intersect_unbounded(rback).front())));
        if (t_right_selfshadow_min < t_right_transparency_max)
            transparency_rate += right_transp*right1.projected_area(d,t_right_selfshadow_min,t_right_transparency_max);
        //t_left_transparency_max calculated above
        float t_left_transparency_min = std::max(-1.0f,left1.parameter_at(rback.at(left1.intersect_unbounded(rback).front())));
        if (t_left_transparency_min < t_left_transparency_max)
            transparency_rate += left_transp*left1.projected_area(d,t_left_transparency_min,t_left_transparency_max);


        if (flip) std::swap(left_rate,right_rate);
        //Normalize rates, although they might cover the cosine we don't want cosine but total values
        float total_rate = barb_rate + left_rate + right_rate + transparency_rate;
        if (std::abs(total_rate)<1.e-5) total_rate=1.0f;
//        std::cerr<<"BARB RATE = "<<barb_rate<<" - LEFT RATE = "<<left_rate<<" - RIGHT RATE = "<<right_rate<<" - TOTAL = "<<total_rate<<std::endl;
        //This should change as we develop this masking term. Basiscally: transmittance
        for (BarbRange& range : barb_ranges) range/=total_rate;
        return Masking(transparency_rate/total_rate,left_rate/total_rate,right_rate/total_rate,barb_ranges);
    }        
};
