#include <iostream>
#include "../flatland.h"
#include "../barbules.h"
#include "svg-barbules.h"

struct Config : public pattern::Reflectable<Config> {
    float excentricity = 1.0f, angle_degrees = 0.0f, separation = 1.0f, scale = 500.0f;
    std::string output = "";
    auto reflect() { return std::tie(excentricity,angle_degrees,separation,scale,output);}
    auto reflect_names() const { return std::tuple("excentricity","angle","separation","scale","output"); }
};

int main(int argc, char** argv) {
    Config config;
    pattern::load_commandline(config,argc,argv);
    float angle = config.angle_degrees*M_PI/180.0;
    flatland::Vec2 dir{-std::sin(angle),-std::cos(angle)};

    svg::SVG file;
    
    Barbules barbules(config.separation,config.excentricity);
    file.add(svg_barbules(barbules,dir)).add_transform(svg::Scale(config.scale));
    file.viewBox(-1.5f*config.scale,
        (-config.excentricity-2.5f)*config.scale,
        (3.0f+config.separation+3.0f)*config.scale,
        (config.excentricity+5.0f)*config.scale);
    if (config.output.empty()) file.save_to_stream(std::cout);
    else file.save_to_file(config.output);
}
