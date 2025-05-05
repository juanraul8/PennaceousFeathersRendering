#pragma once
#include "svg-barbules.h"
#include "../barbs.h"

svg::Polygon svg_shadow(const flatland::Vec2& light, const flatland::segment& s) {
    svg::Polygon polygon;
    polygon.add_point(s.p0());
    polygon.add_point(s.p1());
    polygon.add_point(s.p1()+5.0f*light/(std::abs(light[1])+1.e-4));
    polygon.add_point(s.p0()+5.0f*light/(std::abs(light[1])+1.e-4));

    return polygon.stroke_width(0).fill(svg::black);
}



svg::Group svg_barbs(const Barbs& barbs, const flatland::Vec2& dir) {
    float barbule_width = 1.0/10.0;
    float width_projected_area = 1.0/20.0;
    float barb_border = 1.0f/50.0f;
    float reference_line_width = 1/100.0f;
    float shadow_opacity = 0.5;


    svg::Group file;

    flatland::ellipse barb1(0.0,0.0,1.0f,barbs.excentricity());
    float total_distance = 2.0f*(1.0f + barbs.barbule_length()*std::cos(barbs.barbule_inclination()));
    flatland::ellipse barb2(total_distance,0.0,1.0,barbs.excentricity());
    flatland::Vec2 corner{1.0f+barbs.barbule_length()*std::cos(barbs.barbule_inclination()),
                -barbs.barbule_length()*std::sin(barbs.barbule_inclination())};

    flatland::segment proximal1(barb1.parametric(0),corner);
    flatland::segment distal1(barb2.parametric(M_PI),corner);
    corner[0]+=total_distance;
    flatland::segment proximal2(barb2.parametric(0),corner);
    corner[0]-=2*total_distance;
    flatland::segment distal0(barb1.parametric(M_PI),corner);
    corner[0]+=total_distance;

    auto projection = [&] (const flatland::Vec2& p0, const flatland::Vec2& p1, float displacement = 0.0f) {
        auto ppro = barb2.center() + dir*(barbs.excentricity()*1.2f+3.0f*std::abs(dir[0])+displacement*width_projected_area);
        auto dpro = flatland::Vec2{dir[1],-dir[0]};
        flatland::segment pro(ppro,ppro+dpro);
        flatland::ray r0(p0,dir);
        flatland::ray r1(p1,dir);
        auto ts0 = pro.intersect_unbounded(r0);
        auto ts1 = pro.intersect_unbounded(r1);
        return svg_segment(r0.at(ts0.front()),r1.at(ts1.front()));
    };

    file.add(svg_shadow(-dir,proximal2)).fill_opacity(shadow_opacity*(1.0-barbs.left_transparency())); 
    file.add(svg_segment(proximal2)).stroke(barbule_proximal_color).stroke_width(barbule_width);
    file.add(svg_ellipse(barb2)).fill(barb_color).stroke(svg::black).stroke_width(barb_border);
    file.add(svg_shadow(-dir,barb2)).fill_opacity(shadow_opacity);
    file.add(svg_shadow(-dir,distal1)).fill_opacity(shadow_opacity*(1.0f-barbs.left_transparency()));
    file.add(svg_segment(distal1)).stroke(barbule_distal_color).stroke_width(barbule_width);
    file.add(svg_shadow(-dir,proximal1)).fill_opacity(shadow_opacity*(1.0f-barbs.right_transparency()));
    file.add(svg_segment(proximal1)).stroke(barbule_proximal_color).stroke_width(barbule_width);
    file.add(svg_ellipse(barb1)).fill(barb_color).stroke(svg::black).stroke_width(barb_border);
    file.add(svg_shadow(-dir,barb1)).fill_opacity(shadow_opacity);
//    file.add(svg_shadow(-dir,distal0)).fill_opacity(shadow_opacity*(1.0f-barbs.left_transparency()));
    file.add(svg_segment(distal0)).stroke(barbule_distal_color).stroke_width(barbule_width);



    float theta_ini = barb2.parameter_from_tangent(dir) - M_PI;
    float theta_end = theta_ini + M_PI;
    file.add(svg_point(barb1.parametric(theta_ini)));
    file.add(svg_point(barb1.parametric(theta_end)));
    file.add(svg_point(barb2.parametric(theta_ini)));
    file.add(svg_point(barb2.parametric(theta_end)));
    file.add(svg_point(barb1.parametric(0)));
    file.add(svg_point(barb2.parametric(M_PI)));
    file.add(svg_point(corner));

    flatland::ray r1(barb1.parametric(theta_end),-dir);
    file.add(svg_ray(r1));
    std::list<float> ts = barb2.intersect(r1);

    if (ts.empty()) {
        flatland::ray r2(barb2.parametric(theta_ini),dir);
        file.add(svg_ray(r2));
        std::list<float> ts2 = distal1.intersect(r2);
        if (ts2.empty()) proximal1.intersect(r2);
        file.add(svg_point(r2.at(ts2.front())));
    } else {
        file.add(svg_point(r1.at(ts.front())));
    }

    ts = proximal1.intersect(r1);
    if (!ts.empty()) file.add(svg_point(r1.at(ts.front())));
    ts = distal1.intersect(r1);
    if (!ts.empty()) file.add(svg_point(r1.at(ts.front())));
    
    flatland::ray r3(corner,-dir);
    file.add(svg_ray(r3));
    ts = barb2.intersect(r3);
    if (!ts.empty()) file.add(svg_point(r3.at(ts.front())));
/*
    file.add(svg_segment(barb1.parametric(theta_end),barb1.parametric(theta_end)+dir*100.0f).stroke(svg::yellow).stroke_width(1.0/100.0f));
    file.add(svg_segment(barb2.parametric(theta_end),barb2.parametric(theta_end)+dir*100.0f).stroke(svg::yellow).stroke_width(1.0/100.0f));
    file.add(svg_segment(barb2.parametric(-M_PI),barb2.parametric(-M_PI)+dir*100.0f).stroke(svg::yellow).stroke_width(1.0/100.0f));
    file.add(svg_segment(barb3.parametric(-M_PI),barb3.parametric(-M_PI)+dir*100.0f).stroke(svg::yellow).stroke_width(1.0/100.0f));

    flatland::ray r1(barb1.parametric(theta_ini),-dir);
    std::list<float> ts = barb2.intersect(r1);
    float theta_shadow = -M_PI; //Shadow of the barb
    if (!ts.empty()) {
        float theta_hit = barb2.parameter_at(r1.at(ts.front()));
//        std::cerr<<"HIT at "<<(theta_hit*180.0/M_PI)<<std::endl;
        if (theta_hit > theta_end) theta_hit -= 2.0*M_PI;
        theta_shadow = std::max(theta_shadow,theta_hit);
    }
    file.add(svg_segment(left1.p1(),left1.p1()+dir*100.0f).stroke(svg::yellow).stroke_width(1.0/100.0f));
    file.add(svg_segment(left2.p1(),left2.p1()+dir*100.0f).stroke(svg::yellow).stroke_width(1.0/100.0f));
    flatland::ray r2(left1.p1(),-dir);
    ts = barb2.intersect(r2);
    if (!ts.empty()) {
        float theta_hit = barb2.parameter_at(r2.at(ts.front()));
 //       std::cerr<<"HIT at "<<(theta_hit*180.0/M_PI)<<std::endl;
        if (theta_hit > theta_end) theta_hit -= 2.0*M_PI;
        theta_shadow = std::max(theta_shadow,theta_hit);
    }
//    std::cerr<<(theta_end*180.0/M_PI)<<" - "<<(theta_shadow*180.0/M_PI)<<" - "<<(theta_ini*180.0/M_PI)<<std::endl;
//    file.add(svg_ellipse_range(barb2,theta_shadow, theta_end).stroke(svg::yellow));
//    if (theta_shadow > theta_ini)
//        file.add(svg_ellipse_range(barb2,theta_ini, theta_shadow).stroke(svg::gray));

//    auto projection = barb2.projected_diameter(-dir,3.0f);
    auto masking = barbs.masking(dir);

    auto ppro = barb2.center() + dir*(std::max(barbs.excentricity(),std::abs(barbs.barbule_length()))*2.0f+4.0f*std::abs(dir[0]));
    auto dpro = flatland::Vec2{dir[1],-dir[0]};
    auto projection = flatland::segment(ppro,ppro+dpro);

    flatland::ray rini(barb1.parametric(theta_end),-dir);
    flatland::ray rend(barb2.parametric(theta_end),-dir);

    auto pini = rini.at(projection.intersect_unbounded(rini).front());
    auto pend = rend.at(projection.intersect_unbounded(rend).front());
    file.add(svg_segment(pini,pend)).stroke(svg::yellow);

*/
/*
    file.add(svg_segment(projection)).stroke(svg::gray);
    auto range = masking.barb_ranges().front();
    file.add(svg_segment(projection.parametric(range.hmin()),projection.parametric(range.hmax()))).stroke(svg::yellow);

    float barb_projection_rate = 2.0*range.rate()/(range.hmax()-range.hmin());
    auto pi = projection.p1();
    auto di = projection.p1()-projection.p0();
    auto po = pi + di*(masking.left_rate()/(barb_projection_rate*(1.0f - barbs.left_transparency())));
    file.add(svg_segment(pi,po)).stroke(svg::red);
    pi=po;
    po = pi + di*(masking.right_rate()/(barb_projection_rate*(1.0f - barbs.right_transparency())));
    file.add(svg_segment(pi,po)).stroke(svg::green);
*/
    return file;
} 
