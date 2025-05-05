#include <iostream>
#include <svg.cc/plot/svgplot.h>
#include <patterns/patterns.h>
#include "../feather.h"
#include "../feather-mc.h"
#include <cmath>

struct Config : public pattern::Reflectable<Config> {
    float barb_excentricity = 1.0f,
          barbule_excentricity = 1.0f,
          barbule_separation = 1.0f,
          barbule_rotation = 90.0f,
          barbule_inclination = 0.0f,
          barbule_length = 0.5f;
    unsigned long plot_samples = 100,
                  montecarlo_samples = 0;
    std::string output = "feathers.svg";
    auto reflect() { return std::tie(barb_excentricity,barbule_excentricity,barbule_separation,barbule_rotation,barbule_inclination,barbule_length,plot_samples,montecarlo_samples,output);}
    auto reflect_names() const { 
        return std::tuple("barb-excentricity","barbule-excentricity","barbule-separation","barbule-rotation","barbule-inclination","barbule-length","plot-samples","montecarlo-samples","output"); }
};

template<typename F>
void test(const F& feathers, unsigned long plot_samples, const std::string& title, const std::string& output) {
    svg::plot::SVGPlot plt;
    plt.figsize({800,400});
    auto background = [] (float t, float p) {
        int it = int(t/4.5)-(t<0.0f?1:0); int ip = int(p/4.5)-(p<0.0f?1:0);
        if (((it+ip)%2)==0) return std::array<float,3>{1,1,1};
        else return std::array<float,3>{0.7,0.7,0.7};  
    };
    auto f = [&feathers] (float t, float p) {

//        std::cerr<<"THETA = "<<t<<" - PHI = "<<p<<std::endl;
        float theta = t*M_PI/180.0f;
        float phi = p*M_PI/180.0f;
        Vec3 view{std::cos(phi)*std::sin(theta), std::sin(phi), std::cos(theta)*std::cos(phi)};

//        std::cerr<<"DIR = "<<view<<std::endl;
        auto masking = feathers.masking(view);
        float r = masking.left_range().rate();
        float g = masking.right_range().rate();
        float b = 0.0f;
        for (auto range : masking.barb_ranges()) b += range.rate();
        float a = masking.transmittance_rate();
        if (a >= 1.0) return std::array<float,4>{r,g,b,0.0f};
        return std::array<float,4>{std::min(1.0f,r/(1.0f-a)),std::min(1.0f,g/(1.0f-a)),std::min(1.0f,b/(1.0f-a)),(1.0f-a)};
    };

//    plt.imshow(svg::plot::linspace(-180.0,180.0,2*plot_samples-1),
//           svg::plot::linspace(-88,88,plot_samples),background).interpolation("bicubic");
    plt.imshow(svg::plot::linspace(-180.0,180.0,2*plot_samples-1),
           svg::plot::linspace(-88,88,plot_samples),f).interpolation("bicubic");
    plt.title(title).xlabel("Azimuth w.r.t. barb").ylabel("Inclination w.r.t. barb");
    plt.savefig(output); 
}

int main(int argc, char** argv) {
    Config config;
    pattern::load_commandline(config,argc,argv);
    Feather feathers(
        Vec3{0,-1,0},Vec3{0,0,1},config.barb_excentricity,config.barbule_excentricity,
        config.barbule_rotation*M_PI/180.0f,config.barbule_inclination*M_PI/180.0f,
        config.barbule_length, config.barbule_separation
    );
    std::stringstream sstr;
    sstr<<"Feathers ";
    if (config.montecarlo_samples > 0) sstr<<"(MC "<<config.montecarlo_samples<<")";
    else sstr<<"(analytical)";
    sstr<<" | exc="<<config.barb_excentricity<<", "<<config.barbule_excentricity<<" - len="
            <<config.barbule_length<<" - inc="<<config.barbule_inclination<<" - rot="<<config.barbule_rotation
            <<" - sep="<<config.barbule_separation;
            
    if (config.montecarlo_samples>0){
        FeatherMC feathers(
            Vec3{0,-1,0},Vec3{0,0,1},config.barb_excentricity,config.barbule_excentricity,
            config.barbule_rotation*M_PI/180.0f,config.barbule_inclination*M_PI/180.0f,
            config.barbule_length, config.barbule_separation,config.montecarlo_samples
        );
        test(feathers,config.plot_samples,sstr.str(),config.output);  
    } else {
        Feather feathers(
            Vec3{0,-1,0},Vec3{0,0,1},config.barb_excentricity,config.barbule_excentricity,
            config.barbule_rotation*M_PI/180.0f,config.barbule_inclination*M_PI/180.0f,
            config.barbule_length, config.barbule_separation
        );
        test(feathers,config.plot_samples,sstr.str(),config.output);  
    }        
}
