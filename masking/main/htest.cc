#include <iostream>
#include <svg.cc/plot/svgplot.h>
#include <patterns/patterns.h>
#include "../flatland.h"
#include <cmath>

struct Config : public pattern::Reflectable<Config> {
    float excentricity = 1.0f, angle = 0;
    unsigned long samples = 1000;
    std::string output = "htest.svg";
    auto reflect() { return std::tie(excentricity,angle, samples, output);}
    auto reflect_names() const { return std::tuple("excentricity","angle","samples","output"); }
};

int main(int argc, char** argv) {
    Config config;
    pattern::load_commandline(config,argc,argv);

    std::list<float> hin, hout;
    flatland::Vec2 dir{std::cos(float(M_PI*config.angle/180.0f)),std::sin(float(M_PI*config.angle/180.0f))};
    flatland::ellipse ellip(0,0,1,config.excentricity);
    flatland::segment casting = ellip.projected_diameter(dir);
    for (float h : svg::plot::linspace(-1.0f,1.0f,config.samples)) {
        flatland::ray r(casting.parametric(h),dir);
        auto ts = ellip.intersect(r);
        if (ts.empty()) std::cerr<<"No hit at "<<h<<std::endl;
        else {
            hin.push_back(h);
            hout.push_back(ellip.h_at(r.at(ts.front()),dir));
        }
    }
    
    svg::plot::SVGPlot plt;
    plt.plot(hin,hout);
    plt.xticks({-1.0,0.0,1.0}).yticks({-1.0,0.0,1.0});
    plt.title("Test for h").ylabel("H Out").xlabel("H In").figsize({500,500});
    plt.savefig(config.output);
}
