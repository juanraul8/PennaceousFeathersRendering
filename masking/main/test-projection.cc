#include <iostream>
#include "../flatland.h"
#include <patterns/patterns.h>
#include "svg-barbules.h"

struct Config : public pattern::Reflectable<Config> {
    float excentricity = 1.0f, angle_degrees = 45.0f, scale = 500.0f;
    std::string output = "";
    auto reflect() { return std::tie(excentricity,angle_degrees,scale,output);}
    auto reflect_names() const { return std::tuple("excentricity","angle","scale","output"); }
};


int main(int argc, char** argv) {
    Config config;
    pattern::load_commandline(config,argc,argv);
    float angle = config.angle_degrees*M_PI/180.0;
    std::list<flatland::Vec2> dirs;
    dirs.push_back(flatland::Vec2{-std::sin(angle),-std::cos(angle)});
    dirs.push_back(flatland::Vec2{std::sin(angle),-std::cos(angle)});
    dirs.push_back(flatland::Vec2{-std::sin(angle),std::cos(angle)});
    dirs.push_back(flatland::Vec2{std::sin(angle),std::cos(angle)});

    svg::SVG file;
    flatland::ellipse ellipse(0.0f,0.0f,1.0f,config.excentricity);
    auto& g = file.add(svg::Group());
    g.add(svg_ellipse(ellipse)).fill(svg::black);
    for (auto dir : dirs) {
        auto seg = ellipse.projected_diameter(dir,3.0*std::max(1.0f,config.excentricity));
        g.add(svg_segment(seg)).stroke(svg::ColorRGB(255.0f*(0.5f*dir[0]+1.0f),255.0f*(0.5*dir[1]+1.0f),0.0f));
        g.add(svg_point(seg.parametric(1.0f))).fill(svg::blue);
        g.add(svg_point(seg.parametric(-1.0f))).fill(svg::black);
    }
    g.add_transform(svg::Scale(config.scale));
    file.viewBox((-1.0f - 3.0f*config.excentricity)*config.scale,
        (-1.0f - 3.0f*config.excentricity)*config.scale,
        2.0f*(1.0f + 3.0f*config.excentricity)*config.scale,
        2.0f*(1.0f + 3.0f*config.excentricity)*config.scale);
    if (config.output.empty()) file.save_to_stream(std::cout);
    else file.save_to_file(config.output);
}
