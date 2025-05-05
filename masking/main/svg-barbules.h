#pragma once
#include "../barbules.h"
#include "colors.h"
//Only include what I use to reduce compilation times
#include <svg.cc/src/circle.h>
#include <svg.cc/src/ellipse.h>
#include <svg.cc/src/line.h>
#include <svg.cc/src/poly.h>
#include <svg.cc/src/group.h>
#include <svg.cc/src/svg.h>

svg::Ellipse svg_ellipse(const flatland::ellipse& e) {
    return svg::Ellipse(e.center(0),e.center(1),e.axes(0),e.axes(1)).stroke_width(0).fill(svg::blue); 
}

svg::Polyline svg_ellipse_range(const flatland::ellipse& e, float theta_min, float theta_max) {
    svg::Polyline poly;
//    std::cerr<<(theta_min*180.0/M_PI)<<" - "<<(theta_max*180.0/M_PI)<< " -> ";
    if (theta_max < theta_min) std::swap(theta_min,theta_max);
    if ((theta_max - theta_min) > M_PI) {
        theta_min += 2.0*M_PI; std::swap(theta_min,theta_max);
    }
    if ((std::abs(theta_max - theta_min - M_PI) < 1.e-3) && (theta_min < 0.0)) {
        theta_min += 2.0*M_PI; std::swap(theta_min,theta_max);  
    }
//    std::cerr<<(theta_min*180.0/M_PI)<<" - "<<(theta_max*180.0/M_PI)<<std::endl;
    for (float d = theta_min; d<theta_max; d+=(M_PI/100.0f)) {
        poly.add_point(e.parametric(d));
    }
    poly.add_point(e.parametric(theta_max));
    return poly.stroke_width(1.0f/25.0f).stroke(svg::yellow).fill(svg::none); 
}

svg::Line svg_segment(const flatland::segment& e) {
    return svg::Line(e.p0(0),e.p0(1),e.p1(0),e.p1(1)).stroke_width(1.0f/25.0f).stroke(svg::yellow).fill(svg::none);
}

svg::Line svg_segment(const flatland::Vec2& p0, const flatland::Vec2& p1) {
    return svg::Line(p0[0],p0[1],p1[0],p1[1]).stroke_width(1.0f/25.0f).stroke(svg::yellow).fill(svg::none);
}

svg::Polygon svg_shadow(const flatland::Vec2& light, const flatland::ellipse& e) {
    svg::Polygon polygon;
    float a1 = e.parameter_from_tangent(light);
    for (float d = M_PI; d>0.0f; d-=(M_PI/100.0f)) {
        polygon.add_point(e.parametric(a1+d));
    }
    polygon.add_point(e.parametric(a1));
    polygon.add_point(e.parametric(a1)+5.0f*light/(std::abs(light[1])+1.e-4));
    polygon.add_point(e.parametric(a1+M_PI)+5.0f*light/(std::abs(light[1])+1.e-4));

    return polygon.stroke_width(0).fill(svg::black);
}

svg::Circle svg_point(const flatland::Vec2& point, float size = 1.0f/10.0f, const svg::Color& color = svg::black) {
    return svg::Circle(point[0],point[1],size).stroke_width(0).fill(color);
}

svg::Line svg_ray(const flatland::ray& r, float width = 1.0f/25.0f, const svg::Color& color = svg::black) {
    return svg_segment(r.origin(),r.at(100.0f)).stroke_width(width).stroke(color).fill(svg::none);
}


svg::Group svg_barbules(const Barbules& barbules, const flatland::Vec2& dir, const svg::Color& color = svg::yellow) {
    svg::Color barbule_color = barbule_proximal_color;
    svg::Color diameter_color = 0.5*barbule_color;
    float separation_width = 1.0/20.0;
    float width_projected_area = 1.0/20.0;
    svg::Color separation_color = transmittance_color;
    float barbule_border = 1.0f/50.0f;
    float reference_line_width = 1/100.0f;
    float shadow_opacity = 0.5;



    svg::Group file;

    flatland::segment separation(1.0f,0.0f,1.0f+2.0f*barbules.separation(),0.0f);
    file.add(svg_segment(separation)).stroke(separation_color).stroke_width(separation_width);

    flatland::ellipse ellipse2(2.0f+2.0f*barbules.separation(),0.0f,1.0f,barbules.excentricity());    
    flatland::ellipse ellipse(0.0f,0.0f,1.0f,barbules.excentricity());
    
    float theta_ini = ellipse2.parameter_from_tangent(dir);
    float theta_end = theta_ini + M_PI;

    file.add(svg_ellipse(ellipse2)).fill(barbule_color).stroke_width(barbule_border).stroke(svg::black);
    file.add(svg_point(ellipse2.parametric(theta_ini)));
    file.add(svg_point(ellipse2.parametric(theta_end)));
    file.add(svg_shadow(-dir,ellipse2)).fill_opacity(shadow_opacity);
    file.add(svg_ellipse(ellipse)).fill(barbule_color).stroke_width(barbule_border).stroke(svg::black);
    file.add(svg_point(ellipse.parametric(theta_ini)));
    file.add(svg_point(ellipse.parametric(theta_end)));
    file.add(svg_shadow(-dir,ellipse)).fill_opacity(shadow_opacity);;

    auto ppro = ellipse2.center() + dir*(barbules.excentricity()*1.2f+3.0f*std::abs(dir[0]));
    auto dpro = flatland::Vec2{dir[1],-dir[0]};
    auto projection = flatland::segment(ppro,ppro+dpro);

    flatland::ray r(ellipse.parametric(theta_ini),-dir);
    file.add(svg_ray(r));
    std::list<float> ts = ellipse2.intersect(r);
    file.add(svg_segment(ellipse2.parametric(theta_ini),ellipse2.parametric(theta_ini)+dir*100.0f).stroke(svg::yellow).stroke_width(1.0/100.0f));
    if (!ts.empty()) {
        file.add(svg_point(r.at(ts.front())));
        flatland::segment diameter(ellipse2.parametric(theta_ini),ellipse2.parametric(theta_end));
        file.add(svg_segment(diameter)).stroke(diameter_color).stroke_width(separation_width);
        file.add(svg_point(r.at(diameter.intersect(r).front())));
    
        flatland::ray rini(ellipse.parametric(theta_ini),dir);
        flatland::ray rend(ellipse2.parametric(theta_ini),dir);
        file.add(svg_ray(rini,reference_line_width,barbule_color));
        file.add(svg_ray(rend,reference_line_width,barbule_color));
        auto pini = rini.at(projection.intersect_unbounded(rini).front());
        auto pend = rend.at(projection.intersect_unbounded(rend).front());
        file.add(svg_segment(pini,pend)).stroke(barbule_color).stroke_width(width_projected_area);
    }
    else {
//           file.add(svg_ellipse_range(ellipse2,theta_end, theta_ini).stroke(svg::yellow));
        std::list<float> ts1 = separation.intersect(r);
        flatland::ray r2(ellipse2.parametric(ellipse2.parameter_from_tangent(dir)+M_PI),dir);
        std::list<float> ts2 = separation.intersect(r2);
        file.add(svg_ray(r2));
        file.add(svg_point(r.at(ts1.front())));
        file.add(svg_point(r2.at(ts2.front())));
//        file.add(svg_segment(r.at(ts1.front()),r2.at(ts2.front())).stroke(svg::ColorRGB(0.0f,200.0f,200.0f)));
//        file.add(svg_segment(r.at(ts1.front()),r.at(ts1.front())+dir*100.0f).stroke(svg::ColorRGB(0.0f,200.0f,200.0f)).stroke_width(1.0/100.0f));
//        file.add(svg_segment(ellipse2.parametric(theta_end),ellipse2.parametric(theta_end)+dir*100.0f).stroke(svg::ColorRGB(0.0f,200.0f,200.0f)).stroke_width(1.0/100.0f));
        flatland::ray rini(ellipse.parametric(theta_ini),dir);
        flatland::ray rsplit(r2.at(ts2.front()),dir);
        flatland::ray rend(ellipse2.parametric(theta_ini),dir);
        auto pini = rini.at(projection.intersect_unbounded(rini).front());
        auto psplit = rsplit.at(projection.intersect_unbounded(rsplit).front());
        auto pend = rend.at(projection.intersect_unbounded(rend).front());
        file.add(svg_ray(rini,reference_line_width,separation_color));
        file.add(svg_ray(rend,reference_line_width,barbule_color));
        file.add(svg_segment(psplit,pend)).stroke(barbule_color).stroke_width(width_projected_area);
        file.add(svg_segment(pini,psplit)).stroke(separation_color).stroke_width(width_projected_area);
    }


    return file;
} 
