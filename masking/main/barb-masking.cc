#include <iostream>
#include "../flatland.h"
#include "svg-barbs.h"

struct Config : public pattern::Reflectable<Config> {
    float excentricity = 1.0f, 
          angle_degrees = 0.0f, 
          barbule_inclination = 0.0f,
          barbule_length = 0.5f,
          left_transparency = 0.0f,
          right_transparency = 0.0f,
          scale = 500.0f;
    std::string output = "";
    auto reflect() { return std::tie(excentricity,angle_degrees,barbule_inclination,barbule_length,left_transparency,right_transparency,scale, output);}
    auto reflect_names() const { return std::tuple("excentricity","angle","barbule-inclination","barbule-length","left-transparency", "right-transparency","scale", "output"); }
};


int main(int argc, char** argv) {

    Config config;
    pattern::load_commandline(config,argc,argv);
    float angle = config.angle_degrees*M_PI/180.0;
    flatland::Vec2 dir{-std::sin(angle),-std::cos(angle)};

    svg::SVG file;
    
    Barbs barbs(config.excentricity,config.barbule_inclination*M_PI/180.0,config.barbule_length,config.left_transparency,config.right_transparency);
    file.add(svg_barbs(barbs,dir)).add_transform(svg::Scale(config.scale));
    file.viewBox(-1.5f*config.scale,
        (-config.excentricity-2.0f)*config.scale,
        (1.0f+3.0*config.barbule_length+0.5f)*config.scale,
        (config.excentricity+4.0f)*config.scale);
    if (config.output.empty()) file.save_to_stream(std::cout);
    else file.save_to_file(config.output);
}
