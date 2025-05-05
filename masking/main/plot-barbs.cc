#include <iostream>
#include <svg.cc/plot/svgplot.h>
#include <patterns/patterns.h>
#include "../barbs.h"
#include "../barbs-mc.h"
#include <cmath>

struct Config : public pattern::Reflectable<Config> {
    float excentricity = 1.0f,
          barbule_inclination = 0.0f,
          barbule_length = 0.5f,
          left_transparency = 0.0f,
          right_transparency = 0.0f;
    unsigned long montecarlo_samples = 0,
          montecarlo_barbs = 1000,
          plot_samples = 100;
    bool plot_h = false;
    std::string output = "";
    auto reflect() { return std::tie(excentricity,barbule_inclination,barbule_length,left_transparency,right_transparency,plot_samples,montecarlo_samples, montecarlo_barbs,plot_h,output);}
    auto reflect_names() const { return std::tuple("excentricity","barbule-inclination","barbule-length","left-transparency", "right-transparency","plot-samples","montecarlo-samples","montecarlo-barbs","plot-h","output"); }
};

template<typename B>
void test(const B& barbs, unsigned long plot_samples, const std::string& title, const std::string& output, bool plot_h) {
    const float lw = 4.0f;

    svg::plot::SVGPlot plt;
    std::list<float> angle, barb, left, right, trans, total;
    std::vector<std::list<float>> hmin(3), hmax(3); //This is the maximum number of allowed barb_ranges 
    for (unsigned long ia = 0; ia<plot_samples; ++ia) {
        float a = (2.0f*M_PI)*float(ia)/float(plot_samples-1) - M_PI; 
        auto masking = barbs.masking(flatland::Vec2{-std::sin(a),-std::cos(a)});
        angle.push_back(a*180.0f/M_PI);
        float barb_rate = 0.0f;
        std::size_t i = 0;
        for (auto range : masking.barb_ranges()) { 
            if (i>=hmin.size()) break;
            if (i>=hmax.size()) break;
            hmin[i].push_back(range.hmin());
            hmax[i].push_back(range.hmax());
            barb_rate+=range.rate(); 
//            std::cerr<<hmin[i].back()<<" - "<<hmax[i].back()<<" | ";
            ++i;
        }
        //Aligning
        for (; i<hmin.size();++i) {
            hmin[i].push_back(-1); hmax[i].push_back(-1);
        }
//        std::cerr<<barb_rate<<" | "<<masking.left_rate()<<" | "<<masking.right_rate()<<std::endl;
        barb.push_back(barb_rate);
        left.push_back(masking.left_rate());
        right.push_back(masking.right_rate());
        trans.push_back(masking.transmittance_rate());
        total.push_back(barb_rate+masking.left_rate()+masking.right_rate()+masking.transmittance_rate());
    }
    if (plot_h) {
        plt.subplot(2,1,0).plot(angle,left).linewidth(lw).color("r");
        plt.subplot(2,1,0).plot(angle,right).linewidth(lw).color("g");
        plt.subplot(2,1,0).plot(angle,barb).linewidth(lw).color("b");
        plt.subplot(2,1,0).plot(angle,trans).linewidth(lw).color("c");
        plt.subplot(2,1,0).plot(angle,total).linewidth(lw).color("k");
        for (auto h : hmax) plt.subplot(2,1,1).plot(angle,h).linewidth(lw).color("b");
        for (auto h : hmin) plt.subplot(2,1,1).plot(angle,h).linewidth(0.5*lw).color("k");
        plt.subplot(2,1,0).ylabel("Rate").xlabel("Angle").figsize({800,300});
        plt.subplot(2,1,1).ylabel("H").xlabel("Angle").figsize({800,300});
    } else {
        plt.plot(angle,left).linewidth(lw).color("r");
        plt.plot(angle,right).linewidth(lw).color("g");
        plt.plot(angle,barb).linewidth(lw).color("b");
        plt.plot(angle,trans).linewidth(lw).color("c");
        plt.plot(angle,total).linewidth(lw).color("k");
    }
    
    plt.title(title);
    plt.savefig(output);
}

int main(int argc, char** argv) {
    Config config;
    pattern::load_commandline(config,argc,argv);

    std::stringstream sstr;
    sstr<<"Barbs ";
    if (config.montecarlo_samples > 0) sstr<<"(MC "<<config.montecarlo_samples<<"sp)";
    else sstr<<"(analytical)";
    sstr<<" | exc="<<config.excentricity<<" - blength="
            <<config.barbule_length<<" - binc="<<config.barbule_inclination<<" - lefttr="<<config.left_transparency<<" - righttr="<<config.right_transparency;
    if (config.montecarlo_samples > 0) {
        test(BarbsMC(config.excentricity,config.barbule_inclination*M_PI/180.0,config.barbule_length, config.left_transparency, config.right_transparency,config.montecarlo_barbs,config.montecarlo_samples),config.plot_samples, sstr.str(),config.output,config.plot_h);
    } else {
        test(Barbs(config.excentricity,config.barbule_inclination*M_PI/180.0,config.barbule_length, config.left_transparency, config.right_transparency),config.plot_samples, sstr.str(),config.output, config.plot_h);
    }
}
