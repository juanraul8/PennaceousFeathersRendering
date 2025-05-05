#pragma once

/*
#ifndef __MITSUBA_CORE_VECTOR_H_
#include "math-array.h"
#endif
*/

#include "math-array.h"

namespace flatland {

#ifdef _MATH_ARRAY_H_
using Vec2 = std::array<float,2>;
#endif

/*
#ifdef __MITSUBA_CORE_VECTOR_H_
using Vec2 = mitsuba::Vector2f;
#endif
*/

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

class ray {
    Vec2 o;
    Vec2 d;
public:
    ray(float ox, float oy, float dx, float dy) : o{ox,oy}, d{dx,dy} {}
    template<typename vec>
    ray(const vec& o, const vec& d) : ray(o[0],o[1],d[0],d[1]) {}
    const Vec2& origin() const { return o; }
    const Vec2& direction() const { return d; }
    float origin(std::size_t i) const { return o[i]; }
    float direction(std::size_t i) const { return d[i]; }

    Vec2 at(float t) const { return o+t*d; }
};

class segment {
    Vec2 p0_,p1_;
public:
    constexpr const Vec2& p0() const noexcept { return p0_; }
    constexpr const Vec2& p1() const noexcept { return p1_; }
    constexpr float p0(std::size_t i) const noexcept { return p0()[i]; } 
    constexpr float p1(std::size_t i) const noexcept { return p1()[i]; }

    segment(float x0, float y0, float x1, float y1) : p0_{x0,y0}, p1_{x1,y1} {}
    template<typename vec>
    segment(const vec& v0, const vec& v1) : segment(v0[0],v0[1],v1[0],v1[1]) {}

    //Parametrized with h between 0 and 1
    Vec2 parametric(float h) const {
        return p0()*(0.5f*(1.0f-h)) + p1()*(0.5f*(1.0f+h));
    }

    //Normal is not normalized so it also represents the area differential
    Vec2 normal(float h = 0) const {
        return normalized(Vec2{p0(1)-p1(1),p1(0)-p0(0)});
    }

    //No need for implicit function, I am omitting it for now

    float parameter_at(const Vec2& p) const {
        std::size_t i = 0;
        if (std::abs(p0(1)-p1(1))>std::abs(p0(0)-p1(0))) i = 1;
        return (2.0f*p[i] - p0(i) - p1(i))/(p1(i)-p0(i));
    }

    /*
    float projected_distance(float t0, float t1, const std::array<float,2>& dir) const {
        return (std::min(t1,1.0f)-std::max(t0,0.0f))*dot(dir,normal());
    }
    */

    std::list<float> intersect_unbounded(const ray& ray) const {
        auto n = normal();
        float den = dot(n,ray.direction());
        if (std::abs(den) > 0) {
            return std::list<float>{(dot(n,p0()-ray.origin())/den)};
        }
        return {};
    }

    std::list<float> intersect(const ray& ray) const {
        std::list<float> hits = intersect_unbounded(ray);
        if (!hits.empty()) {
            float t = hits.front();
            if (t>0.0f) {
                float h = parameter_at(ray.at(t));
                if ((h>=-1.0f) && (h<=1.0f)) return std::list<float>{t};
            }
        }
        return {};
    }

    //Returns the area projected towards direction d between two parameters
    //Also works from below
    float projected_area(const Vec2& d, float t0, float t1) const {
        return length(parametric(t1)-parametric(t0))*std::abs(dot(d,normal()));
    }
};

class ellipse {
    Vec2 center_;
    Vec2 axes_;

public:
    constexpr const Vec2& center() const noexcept { return center_; }
    constexpr const Vec2& axes() const noexcept { return axes_; }
    constexpr float center(std::size_t i) const noexcept { return center()[i]; }
    constexpr float axes(std::size_t i) const noexcept { return axes()[i];  }
    ellipse(float c0, float c1, float a0, float a1) : center_{c0,c1}, axes_{a0,a1} {}
    template<typename vec>
    ellipse(const vec& c, const vec& a) : ellipse(c[0],c[1],a[0],a[1]) {}

    Vec2 parametric(float t) const {
        return center() + axes()*Vec2{std::cos(t),std::sin(t)};
    }

    Vec2 normal(float t) const {
        return normalized(Vec2{-axes(1)*std::cos(t),-axes(0)*std::sin(t)});
    }

    float implicit(const Vec2& p) const {
        return sqr((p[0]-center(0))/axes(0)) + sqr((p[1]-center(1))/axes(1)) - 1.0f;
    } 

    Vec2 normal(const Vec2& p) const { 
        return normalized(Vec2{
            2.0f*(p[0]-center(0))/sqr(axes(0)),
            2.0f*(p[1]-center(1))/sqr(axes(1))});
    }

    float parameter_at(const Vec2& p) const {
        return std::atan2((p[1]-center(1))/axes(1),(p[0]-center(0))/axes(0));
    }

    //Returns a list of ts along the ray
    std::list<float> intersect(const ray& ray) const {
        float a = sum(sqr(ray.direction()/axes()));
        float b = 2.0f*sum(ray.direction()*(ray.origin()-center())/sqr(axes()));
        float c = sum(sqr((ray.origin()-center())/axes())) -1.0f;
        float disc = b*b - 4.0f*a*c;
        if (disc < 0) return {};
        else {
            std::list<float> hits;
            float sqrtdisc = std::sqrt(disc);
            float inv2a = 1.0f/(2.0f*a);
            float t;

            t = (-b -sqrtdisc)*inv2a;
            if (t>0) hits.push_back(t);
            t = (-b +sqrtdisc)*inv2a;
            if (t>0) hits.push_back(t);
            return hits;
        }
    }

    //Tangent theta (this is the first solution, while the second is the result of this plus pi)
    float parameter_from_tangent(const Vec2& d) const {
        return std::atan((-axes(1)*d[0])/(axes(0)*d[1]));
    }

    //Returns the segment for the origin points of a ray that follows the direction -d, parametrized by h between -1 and 1
    segment projected_diameter(const Vec2& d, float distance = 10.0f) const {
        float angle = parameter_from_tangent(d);
        Vec2 p0 = parametric(angle); 
        Vec2 p1 = parametric(angle + M_PI);
        //This below helps for extreme symmetry and identifying where is the -1 and 1
        if ((p1[0]!=p0[0]) && ((sgn(d[1])*std::atan((p1[1]-p0[1])/(p1[0]-p0[0])))<0)) std::swap(p0,p1);
        Vec2 displacement = -d*(distance*std::max(axes(0),axes(1))/length(d)); 
        return segment(p0 + displacement, p1 + displacement);
    }

    //Returns the area projected towards direction d between two parameters
    float projected_area(const Vec2& d, float t0, float t1) const {
//        while (t1<=t0)  t1 += 2.0*M_PI
//        std::cerr<<"  "<<d<<"  ["<<(t0*180.0/M_PI)<<","<<(t1*180.0/M_PI)<<"] ";
//        std::cerr<<d[0]<<"*("<<std::sin(t0)<<" - "<<std::sin(t1)<<") + "
//                <<d[1]<<"*("<<std::cos(t1)<<" - "<<std::cos(t0)<<")\n";
        return std::abs(d[0]*axes(1)*(std::sin(t0)-std::sin(t1))+d[1]*axes(0)*(std::cos(t1)-std::cos(t0)));
    }

    float h_at(const Vec2& p, const Vec2& d) const {
        auto seg = projected_diameter(d);
        ray r(p,d);
        return seg.parameter_at(r.at(seg.intersect_unbounded(r).front()));
    }
};

template<typename T>
std::tuple<float,const T*> intersect(const std::list<T>& objects, const ray& r) {
    std::tuple<float,const T*> sol(std::numeric_limits<float>::infinity(),nullptr);
    for (const T& o : objects) {
        std::list<float> ts = o.intersect(r);
        if ( (!ts.empty()) && (ts.front() < std::get<0>(sol)) ) {
            sol = {ts.front(),&o};
        }
    }
    return sol;
}

}
