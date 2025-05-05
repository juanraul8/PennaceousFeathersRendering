#include <iostream>
#include <svg.cc/plot/svgplot.h>
#include <patterns/patterns.h>
#include "../barbules.h"
#include "../barbules-mc.h"
#include <cmath>

struct Config : public pattern::Reflectable<Config> {
    float excentricity = 1.0f,
          separation = 0.5f;
    unsigned long montecarlo_samples = 0,
          montecarlo_barbules = 1000,
          plot_samples = 100;
    bool plot_h = false;
    std::string output = "";
    auto reflect() { return std::tie(excentricity,separation,plot_samples, montecarlo_samples, montecarlo_barbules, plot_h, output);}
    auto reflect_names() const { return std::tuple("excentricity","separation","plot-samples", "montecarlo-samples","montecarlo-barbules","plot-h","output"); }
};

template<typename B>
void test(const B& barbules, unsigned long plot_samples, const std::string& title, const std::string& output, bool plot_h) {
    const float lw = 4.0f;

    svg::plot::SVGPlot plt;
    std::list<float> angle, barbule, transmittance, total;
    std::list<float> hmin, hmax;
    for (unsigned long ia = 0; ia<plot_samples; ++ia) {
        float a = (2.0f*M_PI)*float(ia)/float(plot_samples-1) - M_PI; 
        auto masking = barbules.masking(flatland::Vec2{-std::sin(a),-std::cos(a)});
        angle.push_back(a*180.0f/M_PI);
        barbule.push_back(masking.barbule_rate());
        transmittance.push_back(masking.transmittance_rate());
        hmin.push_back(masking.hmin());
        hmax.push_back(masking.hmax());
        total.push_back(masking.barbule_rate() + masking.transmittance_rate());
    }
    if (plot_h) {
        plt.subplot(2,1,0).plot(angle,barbule).linewidth(lw).color("y");
        plt.subplot(2,1,0).plot(angle,transmittance).linewidth(lw).color("c");
        plt.subplot(2,1,0).plot(angle,total).linewidth(lw).color("k");

        plt.subplot(2,1,1).plot(angle,hmin).linewidth(lw).color("y");
        plt.subplot(2,1,1).plot(angle,hmax).linewidth(lw).color("y");
        plt.subplot(2,1,0).ylabel("Rate").xlabel("Angle").figsize({800,300});
        plt.subplot(2,1,1).ylabel("H").xlabel("Angle").figsize({800,300});
    } else {
    plt.plot(angle,barbule).linewidth(lw).color("y");
    plt.plot(angle,transmittance).linewidth(lw).color("c");
    plt.plot(angle,total).linewidth(lw).color("k");
    plt.ylabel("Rate").xlabel("Angle").figsize({800,300});
    }
    plt.title(title);
    plt.savefig(output);
}

int main(int argc, char** argv) {
    Config config;
    pattern::load_commandline(config,argc,argv);
    
    std::stringstream sstr;
    sstr<<"Barbules ";
    if (config.montecarlo_samples > 0) sstr<<"(MC "<<config.montecarlo_samples<<"sp)";
    else sstr<<"(analytical)";
    sstr<<" | excentricity="<<config.excentricity<<" - separation="<<config.separation;
    if (config.montecarlo_samples > 0) {
        test(BarbulesMC(config.separation,config.excentricity,config.montecarlo_barbules,config.montecarlo_samples), config.plot_samples,sstr.str(),config.output, config.plot_h);
    } else {
        test(Barbules(config.separation,config.excentricity), config.plot_samples, sstr.str(),config.output, config.plot_h);
    }
}
