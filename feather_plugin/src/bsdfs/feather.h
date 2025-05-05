/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2014 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

/*
A feather BSDF is a mixture of materials based on the particular hierarchical structure of feathers: rachis, barb and barbules. 

References:
[0] Diffuse BSDF: Smooth diffuse material (diffuse.cpp) (debug material)
[1] Hair BSDF: https://www.pbrt.org/hair.pdf --> hair.cpp (PBRT 3)
[2] Fur BSDF: https://sites.cs.ucsb.edu/~lingqi/project_page/fur2/index.html --> fur.cpp 
[3] Mixture of BSDFs: Mixture material (mixturebsdf.cpp) (https://github.com/mitsuba-renderer/mitsuba/blob/cfeb7766e7a1513492451f35dc65b86409655a7b/src/bsdfs/mixturebsdf.cpp)
[4] Opacity Mask: Opacity mask material (mask.cpp)(https://github.com/mitsuba-renderer/mitsuba/blob/master/src/bsdfs/mask.cpp)
[5] Two-sided BRDF: Two-sided BRDF adapter (twosided.cpp)(https://github.com/mitsuba-renderer/mitsuba/blob/master/src/bsdfs/twosided.cpp)
*/ 

#pragma once

#if !defined(__FEATHER_H)
#define __FEATHER_H

#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/texture.h>
#include <mitsuba/render/scene.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/hw/basicshader.h>
#include <mitsuba/core/warp.h>
#include "ior.h"

//Bark/Masking
#include <cmath>
#include <optional>
#include <array>
#include <tuple>
#include "barbs.h"
#include "barbules.h"
#include "feather_mask.h"

#include <mutex> 

MTS_NAMESPACE_BEGIN

//Debugging
Float w_min = 100;
Float w_max = -100;

// ----------------------------------------------------Feather BSDF Declarations---------------------------------------------------------------------
class Feather : public BSDF {
    
    public:

    // Feather Public Methods
    Feather(const Properties &props);
    Feather(Stream *stream, InstanceManager *manager);
    
    void configure();
    void serialize(Stream *stream, InstanceManager *manager) const;

//    Frame getFrame(const Intersection &its) const;
    Frame getFrame(Float barb_azimuthal, Float barbule_azimuthal, Float barbule_longitudinal,BSDFSamplingRecord bRec ) const;
    
    void setLocalFrames(Frame &barb_f, Frame &left_barbule_f,Frame &right_barbule_f, Float barb_angle) const;
 
    Transform getTransform(Float barb_azimuthal, Float azimuthal, Float longitudinal) const;
    
    Frame getBarbuleFrame(Float barbule_angle, const Intersection &its) const;

    bool isFeather(const Intersection &its, Float *texture, Float *color) const;
    
    std::vector<Float> getMaskWeights(Frame const barb_f, Frame const left_barbule_f,Frame const right_barbule_f, Vector wi_dir) const ;

    Spectrum eval(const BSDFSamplingRecord &bRec, EMeasure measure) const;
    Spectrum featherEval(const BSDFSamplingRecord &bRec, EMeasure measure, Float texture_mask) const;
    Float pdf(const BSDFSamplingRecord &bRec, EMeasure measure) const;
    Float featherPdf(const BSDFSamplingRecord &bRec, EMeasure measure, Float texture_mask) const;
    Spectrum sample(BSDFSamplingRecord &bRec, const Point2 &sample) const;
    Spectrum sample(BSDFSamplingRecord &bRec, Float &pdf, const Point2 &sample) const;
    Spectrum featherSample(BSDFSamplingRecord &bRec, Float &pdf, const Point2 &_sample, Float texture_mask) const;
    
    //Masking functions
    Float getBarbuleMasks(BSDFSamplingRecord bRec,Float barb_angle, Float barbule_angle_azimuthal, Float barbule_angle_longitudinal,
    Barbules::Masking &masking_right_barbule, Barbules::Masking &masking_left_barbule,Barbs &barbs ) const;
    Feather_masking::Masking getMasking(Frame const barb_f, Frame const left_barbule_f,Frame const right_barbule_f, Vector wi_dir) const;

    //Virtual measurements
    std::string getExperimentFileName(std::string experiment_name);
    void getVectorMeasurements(Float phi_i, Float theta_i, Float phi_o, Float theta_o, Vector &wi, Vector &wo);
    Spectrum measureBRDF(Float phi_i, Float theta_i, Float phi_o, Float theta_o);
    
    //void latLongQuery(int n_phi_o, int n_theta_o, int step_phi_o, int step_theta_o, int N, const char* file_path);

    void latLongQuery(std::vector<Float> phi_i_values, std::vector<Float> theta_i_values, std::vector<Float> phi_o_values, std::vector<Float> theta_o_values,
                      int N, const char* file_path);
    void latLongQueryParallel(int n_phi_o, int n_theta_o, int step_phi_o, int step_theta_o, int N, char* file_path);

    //Debugging functions
    Spectrum debuggingGeometry(const BSDFSamplingRecord &bRec, EMeasure measure) const;
    void getTransoformedRec(Float barb_angle, Float barbule_angle_azimuthal, Float barbule_angle_longitudinal,BSDFSamplingRecord bRec, BSDFSamplingRecord &bRec_barb,BSDFSamplingRecord &bRec_left_barbule, BSDFSamplingRecord &bRec_right_barbule) const;
    
    flatland::Vec2 getBarbDir(const Frame transformed_frame, const BSDFSamplingRecord &bRec,BSDFSamplingRecord &bRecTransformed ) const;
    Spectrum evalTransmittanceEvent(const BSDFSamplingRecord &bRec, EMeasure measure) const;
    Float pdfTransmittanceEvent(const BSDFSamplingRecord &bRec, EMeasure measure) const;
    Spectrum sampleTransmittanceEvent(BSDFSamplingRecord &bRec, Float &pdf, const Point2 &_sample) const;
    Spectrum getDiffuseReflectance(const Intersection &its) const;
    
    Float getBarbAngle(const Intersection its) const;
    void addChild(const std::string &name, ConfigurableObject *child);
    std::string toString() const;

    MTS_DECLARE_CLASS()

    private:
    
    //BRDR materials
    std::vector<BSDF*> m_bsdfs;//Order: rachis, barb, barbule

    // Masking / shadowing;
    Barbules barbules;
    Vector barb_axis;
    Float barbule_length = 1.0f; 
    Float barbule_sep = 1.0f; 

    //NB: Eccentricity values should be computed from the corresponding pigmentation BSDF (barb/barbule) in order to be consistent
    //Float barbule_exc; 
    //Float barb_exc; 

    Float barbule_angle_azimuthal, barbule_angle_longitudinal;

    //Textures
    ref<Texture> m_texture; //Geometry texture
    ref<Texture> m_masking_texture; //Geometry texture
    ref<Texture> m_color; //Color texture

    //Feather masks --> Background is always 0.0
    Float shaft_mask;
    Float barb_mask;
    //Float barb_angle;
    Float barb_angle2;
    
    // Color variables
    Spectrum m_color_rachis; // Rachis color
    Spectrum barbs_rgb; // when rgb color passed 

    // Debugging variables
    bool use_masking = false;//Enable/disable the masking term. If it is false, we only render the Barb BCSDF.
    bool debugging_geometry = false;
    bool debugging_bsdf = false; //Debugging the contribution of each term using a particular color (Barb: Blue, Red: Left Barbule, Green: Right Barbule)
    Spectrum barb_color, left_barbule_color, right_barbule_color;
    mutable std::mutex mutex_feather;//Careful with this!!!

    // Virtual measurements
    bool ferrero_setup = false;// In our measurements, the outgoing azimuthal angle depends on the incident angle.
    bool measure_eval = false;// If it is enabled, we measure this BSDF eval function.
    bool measure_fiber = false;
    Float virtual_barb_angle = 0.0f; // Arbitrary barb angle we use for the virtual measurements
    bool hit_rachis = false; // In case, we need to consider the rachis radiance too.
    int n_experiments = 1000;
    std::string experiment_name = "test";//We use this variable to identify the free variable we are testing from the masking.
    std::string measurements_file = "test";//Output file for BRDF virtual measurements.
    std::string angles_file = "test";//Input angles file for BRDF virtual measurements.
    int measurements_seed = -1;//Seed for BRDF virtual measurements.-1: denotes original seed.
    bool fitting = false;

    Float eval_max = 1000.0f;//Maximum value to clamp the BCSDF evaluation
    Float sample_max = 2.0f;//Maximum value to clamp the weights of the BCSDF sampling
};

MTS_IMPLEMENT_CLASS_S(Feather, false, BSDF)
MTS_EXPORT_PLUGIN(Feather, "Feather BRDF")

MTS_NAMESPACE_END

#endif /* __FEATHER_H */