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

#include "feather.h"
#include <omp.h>

MTS_NAMESPACE_BEGIN

//----------------------Coordinate system auxiliar transformations-------------------------

Vec3 inline getTransformedVector(Float barb_angle, Vector input_vector)
{
   
    Transform Rn = Transform::rotate(Vector(0,0,1), -barb_angle); 
    
    Vector temp_view = normalize(Rn.inverse()(input_vector));
    Vector t = temp_view;
    
    Vec3 view{t[0],t[1],-t[2]};
    return view;
}

// Barb/Barbule local coordinate system transformation for internal shading computations (new)
Transform Feather::getTransform(Float barb_azimuthal, Float barbule_azimuthal, Float barbule_longitudinal) const
{

    // We would like to compute the inverse transformation (world coordinate to local coordinates): local transformation (elevation and rotation around the resulting plane)
    // NB: Matrices muliplication order (right to left order) --> Order was validated with MATLAB.
    // NB: Euler angles order R = RxRyRz

    // Compute inverse using negative angles (inverse rotation)
    Vector axis = Vector(0.0, 0.0, 1.0);
    //printf("rot along z %.2f",-barb_azimuthal - barbule_azimuthal);
    Transform Rn = Transform::rotate(axis, -barb_azimuthal - barbule_azimuthal);
    axis = Vector(1.0, 0.0, 0.0);
    // rotation along the local x axis
    Transform Rt = Transform::rotate(axis, -barbule_longitudinal);
    // Transform frame = Rt;//Rn * Rt;
    Transform _transform = Rt * Rn; // Change order due to the inverse product identity

    return _transform;
}

// Barbule local coordinate system for internal shading computations (original)
Frame Feather::getBarbuleFrame(Float barbule_angle, const Intersection &its) const
{
    Frame result;
    Vector tex;
    //get the geometry texture rgb values
    //printf("barbule_angle %f\n", barbule_angle);

    m_texture->eval(its, false).toLinearRGB(tex.x, tex.y, tex.z);
    //m_texture->eval(its, false).toSRGB(tex.x, tex.y, tex.z);
    barbule_angle = tex.x;

    // printf("barbule_angle %f\n", barbule_angle);

    Float theta = barbule_angle * M_PI / 180; // degrees to radians

    // recover tangent vector
    // NB: When s is defined first, the transformation is similar to Transform::rotate(axis, barb_angle);
    Vector t;
    t.x = std::cos(theta);
    t.y = std::sin(theta);
    t.z = 0.0; // we don't use z here --> blue channel correspond to the curve mask (background, shaft, barbs)

    Vector s(-t.y, t.x, 0.0); // bitangent vector

    // tangent plane along the curve
    result.t = normalize(t);
    result.s = normalize(s);
    result.n = normalize(cross(result.t, result.s));

    return result;
}

// Barb local coordinate system for internal shading computations --> rotate coordinate system based on feather tangents
// Reference code: bumpmapcpp, normalmap.cpp
Frame Feather::getFrame(Float barb_azimuthal, Float barbule_azimuthal, Float barbule_longitudinal,BSDFSamplingRecord bRec) const
{
    //Transform m_transform = getTransform(barb_azimuthal, barbule_azimuthal, barbule_longitudinal);
    barb_azimuthal = 0.0f;
    
    Vector n = bRec.its.shFrame.n;
    Vector s = bRec.its.shFrame.s;
    Vector t = bRec.its.shFrame.t;
    Transform Rn = Transform::rotate(n, -barb_azimuthal - barbule_azimuthal);
    Transform Rt = Transform::rotate(s, -barbule_longitudinal);
    Transform m_transform = Rt * Rn;
    Vector m_x = normalize(m_transform(s));
    Vector m_y = normalize(m_transform(t));
    Vector m_z = normalize(m_transform(n));
    //printf(" mx %s my %s mz %s n %s s %s t %s \n",m_x.toString().c_str(),m_y.toString().c_str(),m_z.toString().c_str(),n.toString().c_str(),s.toString().c_str(),t.toString().c_str());
    
    return Frame(m_x, m_y, m_z);
}

inline Vector getCorrectFrame(const BSDFSamplingRecord &bRec, Vector input_vector){
    
    //return input_vector;
    //return input_vector;
    //return input_vector;
    //check pigmentation code for original implentation
    const Frame shFrame = bRec.its.shFrame;
    return shFrame.toWorld(input_vector);
    //const Frame geoFrame = bRec.its.geoFrame;
    //return geoFrame.toWorld(input_vector);
    //Vector tmpvector = geoFrame.toWorld(input_vector);
    
    //Vector tmpvector = geoFrame.toWorld(input_vector);
    ////tmpvector = shFrame.toLocal(tmpvector);
    //return tmpvector;
    
}


// Frame frame_barb, Frame frame_left_barbule,Frame frame_right_barbule
// Mitsuba works with angles in degrees? 
inline Frame getNewFrame(Float barbule_angle_longitudinal, Float barbule_angle_azimuthal, Float barb_angle){
    
    Vector axis(0.0, 0.0, 1.0);
    Transform Rn = Transform::rotate(axis, barb_angle);

    Transform _transform = Rn; // Change order due to the inverse product identity
    
    Vector local_x = _transform.inverse()(Vector(1,0,0)); // (Vector(1,0,0));
    Vector local_y = _transform.inverse()(Vector(0,1,0)); // (Vector(0,-1,0));
    Vector local_z = _transform.inverse()(Vector(0,0,1));

    return  Frame(normalize(local_x),normalize(local_y),normalize(local_z));
}

inline Frame getBarbulesFrame(Float barbule_angle_longitudinal,Float barbule_angle_azimuthal ){

    Vector axis(0.0, 0.0, 1.0);// rotation along the local z axis
    Transform Rn = Transform::rotate(axis, -barbule_angle_azimuthal);

    axis = Vector(1.0, 0.0, 0.0);// rotation along the local x axis
    Transform Rt = Transform::rotate(axis, -barbule_angle_longitudinal);

    Transform _transform = Rn * Rt; // Change order due to the inverse product identity
    
    Vector local_x = _transform.inverse()(Vector(1,0,0)); // (Vector(1,0,0));
    Vector local_y = _transform.inverse()(Vector(0,1,0)); // (Vector(0,-1,0));
    Vector local_z = _transform.inverse()(Vector(0,0,1)); // (Vector(0,0,-1));

    return  Frame(normalize(local_x),normalize(local_y),normalize(local_z));
}

void Feather::setLocalFrames(Frame &barb_f, Frame &left_barbule_f,Frame &right_barbule_f, Float barb_angle) const  {
    
    //barb_angle = 0.0f;
    barb_f = getNewFrame(0.0f,0.0f,-barb_angle);

    // Original transformation: latlong reflectance plots are asymmetric with discontinuity artifacts      
    //left_barbule_f = getBarbulesFrame(barbule_angle_longitudinal, -barbule_angle_azimuthal - barb_angle);
    //right_barbule_f = getBarbulesFrame(barbule_angle_longitudinal, barbule_angle_azimuthal - barb_angle);

    //left_barbule_f = getBarbulesFrame(-barbule_angle_longitudinal, -barbule_angle_azimuthal - barb_angle);
    //right_barbule_f = getBarbulesFrame(-barbule_angle_longitudinal, barbule_angle_azimuthal - barb_angle);

    //left_barbule_f = getBarbulesFrame(barbule_angle_longitudinal, -barbule_angle_azimuthal - barb_angle);
    //right_barbule_f = getBarbulesFrame(-barbule_angle_longitudinal, barbule_angle_azimuthal - barb_angle);

    left_barbule_f = getBarbulesFrame(-barbule_angle_longitudinal, -barbule_angle_azimuthal - barb_angle);
    right_barbule_f = getBarbulesFrame(barbule_angle_longitudinal, barbule_angle_azimuthal - barb_angle);
}

flatland::Vec2 Feather::getBarbDir(const Frame transformed_frame, const BSDFSamplingRecord &bRec, BSDFSamplingRecord &bRecTransformed) const
{
    flatland::Vec2 barb_dir;
    
//    getCorrectFrame(const BSDFSamplingRecord &bRec, Vector input_vector)
    
    Vector v1 = getCorrectFrame(bRec,bRec.wi);
    //if (dot(v1, bRec.its.shFrame.n) < 0)
    //        v1[2] *= -1;
    Vector v2 = getCorrectFrame(bRec,bRec.wo);
    //if (dot(v2, bRec.its.shFrame.n) < 0)
    //        v2[2] *= -1;
    bRecTransformed.wi = transformed_frame.toLocal(v1);
    bRecTransformed.wo = transformed_frame.toLocal(v2);

    barb_dir[0] = bRecTransformed.wi[0];
    barb_dir[1] = bRecTransformed.wi[2];
    return barb_dir;
}

//------------------------------Masking auxiliar functions --------------------------------

inline Float computeEccentricity(Float a, Float b){    
    Float e = 0.0f;
    bool debug = false;

    if (a > b){
        e = sqrt(1.0 - (b*b)/(a*a));
    }else{
        e = sqrt(1.0 - (a*a)/(b*b));
    }

    if (debug){
        printf("Ellipse (eccentricity) a = %f, b = %f, e = %f\n", a, b, e);
    }

    return e;
} 

std::vector<Float> Feather::getMaskWeights(Frame const barb_f, Frame const left_barbule_f,Frame const right_barbule_f, Vector wi_dir) const  {
   
    bool debug = false;
   
    Vector leftdir = left_barbule_f.toLocal(wi_dir);
    Vec3 left_barbules_dir{leftdir[0],-leftdir[2],-leftdir[1]};
    
    Vector rightdir = right_barbule_f.toLocal(wi_dir);
    Vec3 right_barbules_dir{rightdir[0],-rightdir[2],rightdir[1]}; //previously the last component was positive
    
    //Vec3 right_barbules_dir{rightdir[0],rightdir[2],rightdir[1]};
    Vector barbdir = barb_f.toLocal(wi_dir);
    Vec3 barb_dir{barbdir[0],-barbdir[2],barbdir[1]};
    //Vec3 barb_dir{barbdir[0],barbdir[1],barbdir[2]};
    
    //Barb and barbule eccentricity should be computed based on the underlying BSDF parameters.
    //Only makes sense for Pigmentation BCSDF
    Float barb_eccentricity = computeEccentricity(m_bsdfs[1]->a_outer, m_bsdfs[1]->b_outer);
    Float barbule_eccentricity = computeEccentricity(m_bsdfs[2]->a_outer, m_bsdfs[2]->b_outer);

    Feather_masking feather_m(Vec3{0,-1,0},Vec3{0,0,1},barb_eccentricity, barbule_eccentricity,
                            barbule_angle_azimuthal * M_PI/180.0f, barbule_angle_longitudinal*M_PI/180.0f,
                            barbule_length, barbule_sep);
    
    auto masking = feather_m.setDirs(barb_dir,left_barbules_dir,right_barbules_dir);
    
    // printf("widir %f %f %f ",barbdir[0],barbdir[1],barbdir[2]);
    // printf("\n");
    
    if(debug){
        printf("v2 %.2f %.2f %.2f",wi_dir[0],wi_dir[1],wi_dir[2]);
        printf(" leftdir2 %f %f %f ",left_barbules_dir[0],left_barbules_dir[1],left_barbules_dir[2]);    
        printf(" rightdir2 %f %f %f ",right_barbules_dir[0],right_barbules_dir[1],right_barbules_dir[2]);
        printf( " barbdir3 %.2f %.2f %.2f",barb_dir[0],barb_dir[1],barb_dir[2]);
    }

    std::vector<Float> m_weights;  
    float total_b = 0.0f;

    for (auto range : masking.barb_ranges()){
        total_b += range.rate();
    }
    
    m_weights.push_back(total_b);
    m_weights.push_back(masking.left_range().rate());
    m_weights.push_back(masking.right_range().rate());
    m_weights.push_back(masking.transmittance_rate());
    
    if(debug){
        printf(" w ");
        printf(" %f %f %f %f ",m_weights[0],m_weights[1],m_weights[2],m_weights[3]);
        printf(" \n");
    }

    return m_weights;        
}

Feather_masking::Masking Feather::getMasking(Frame const barb_f, Frame const left_barbule_f,Frame const right_barbule_f, Vector wi_dir) const  {
   
    bool debug = false;
   
    Vector leftdir = left_barbule_f.toLocal(wi_dir);
    Vec3 left_barbules_dir{leftdir[0],-leftdir[2],-leftdir[1]};
    
    Vector rightdir = right_barbule_f.toLocal(wi_dir);
    Vec3 right_barbules_dir{rightdir[0],-rightdir[2],rightdir[1]}; //previously the last component was positive
    
    //Vec3 right_barbules_dir{rightdir[0],rightdir[2],rightdir[1]};
    Vector barbdir = barb_f.toLocal(wi_dir);
    Vec3 barb_dir{barbdir[0],-barbdir[2],barbdir[1]};
    //Vec3 barb_dir{barbdir[0],barbdir[1],barbdir[2]};

    //printf("getMasking = barb = %s, left_barbule = %s, right_barbule = %s\n", barbdir.toString().c_str(), leftdir.toString().c_str(), rightdir.toString().c_str());
    
    //Barb and barbule eccentricity should be computed based on the underlying BSDF parameters. Only makes sense for Pigmentation BCSDF
    //Float barb_eccentricity = computeEccentricity(m_bsdfs[1]->a_outer, m_bsdfs[1]->b_outer);
    //Float barbule_eccentricity = computeEccentricity(m_bsdfs[2]->a_outer, m_bsdfs[2]->b_outer);

    //Masking code is assumming a = 1.0 (minor axis) and the eccentricity parameter is actually the major axis (b)
    Float barb_eccentricity = m_bsdfs[1]->b_outer;
    Float barbule_eccentricity = m_bsdfs[2]->b_outer;
    
    //printf("Barbule separation = %f, length = %f, azimuthal angle = %f, longitudinal angle = %f, barb eccentricity = %f, barbule eccentricity = %f\n", barbule_sep, barbule_length, barbule_angle_azimuthal, 
    //        barbule_angle_longitudinal, barb_eccentricity, barbule_eccentricity);

    Feather_masking feather_m(Vec3{0,-1,0},Vec3{0,0,1}, barb_eccentricity, barbule_eccentricity,
                            degToRad(barbule_angle_azimuthal), degToRad(barbule_angle_longitudinal),
                            barbule_length, barbule_sep);
    
    Feather_masking::Masking masking = feather_m.setDirs(barb_dir,left_barbules_dir,right_barbules_dir);
    
    return masking;       
}

Feather::Feather(const Properties &props) : BSDF(props)
{
    // props is an associative map where the BSDF configuration is stored --> helpful to obtain the values of the xml configuration file

    printf("Creating Feather BRDF\n");

    // Masking-shadowing parameters
    // barb_w = 0.8f;barbule_left_w = 0.1f; barbule_right_w = 0.1f;
    barbule_angle_azimuthal = props.getFloat("barbule_azimuthal_angle", 45.0f);
    barbule_angle_longitudinal = props.getFloat("barbule_longitudinal_angle", 0.0f);
    barbule_sep = props.getFloat("barbule_separation", 1.0f);
    barbule_length = props.getFloat("barbule_length", 1.0f);
    
    //printf("barbule_longitudinal_angle = %f\n", barbule_angle_longitudinal);
    //printf("barbule separation = %f\n", barbule_sep);

    // Texture masks
    shaft_mask = props.getFloat("shaft_mask", 0.5f);
    barb_mask = props.getFloat("barb_mask", 1.0f);

    // Debugging parameters
    use_masking = props.getBoolean("use_masking", true);
    debugging_geometry = props.getBoolean("debugging_geometry", false);
    debugging_bsdf = props.getBoolean("debugging_bsdf", false);

    barb_color[2] = 1.0f;
    left_barbule_color[0] = 1.0f;
    right_barbule_color[1] = 1.0f;

    eval_max = props.getFloat("eval_max", 1000.0f);
    sample_max = props.getFloat("sample_max", 2.0f);
    
    Float barbule_eccentricity = 0.0f;
    barbules = Barbules(barbule_sep, barbule_eccentricity);
    
    barb_angle2 = props.getFloat("barb_angle", 0.0f);

    measure_eval = props.getBoolean("measure_eval", false);

    //printf("measure_eval = %d\n", measure_eval);

    if (measure_eval){
        virtual_barb_angle = props.getFloat("barb_angle", 30.0f);//Only debugging purpose, barb angle should be obtained from the geometry texture
        hit_rachis = props.getBoolean("hit_rachis", false);
        n_experiments = props.getInteger("n_experiments", 1000);
        experiment_name = props.getString("experiment_name", "");
        measurements_file = props.getString("measurements_file", "");
        angles_file = props.getString("angles_file", "");
        ferrero_setup = props.getBoolean("ferrero_setup", false);

        measurements_seed = props.getInteger("measurements_seed", -1);

        if (measurements_seed != -1){

            printf("Setting seed = %d\n", measurements_seed);

            std::srand(measurements_seed);
        }
    }
}

Feather::Feather(Stream *stream, InstanceManager *manager) : BSDF(stream, manager)
{

    m_texture = static_cast<Texture *>(manager->getInstance(stream));
    m_color = static_cast<Texture *>((manager->getInstance(stream)));
    m_masking_texture = static_cast<Texture *>(manager->getInstance(stream));

    m_bsdfs.push_back(static_cast<BSDF *>(manager->getInstance(stream)));
    m_bsdfs.push_back(static_cast<BSDF *>(manager->getInstance(stream)));
    m_bsdfs.push_back(static_cast<BSDF *>(manager->getInstance(stream)));
    m_bsdfs.push_back(static_cast<BSDF *>(manager->getInstance(stream)));
    m_bsdfs[0]->incRef();//Rachis
    m_bsdfs[1]->incRef();//Barb
    m_bsdfs[2]->incRef();//Left barbule
    m_bsdfs[3]->incRef();//Right barbule

    configure();
}

std::string Feather::getExperimentFileName(std::string experiment_name){

    std::string result = "";
    char variable[100];

    printf("experiment_name = %s\n", experiment_name.c_str());

    if (experiment_name == "barb_angle"){
        sprintf(variable, "barb_angle_%0.2f", virtual_barb_angle);
    }

    if (experiment_name == "barbule_azimuthal_angle"){
        sprintf(variable, "barbule_azimuthal_angle_%0.2f", barbule_angle_azimuthal);
    }

    if (experiment_name == "barbule_longitudinal_angle"){
        sprintf(variable, "barbule_longitudinal_angle_%0.2f", barbule_angle_longitudinal);
    }

    if (experiment_name == "barbule_separation"){
        sprintf(variable, "barbule_separation_%0.2f", barbule_sep);
    }

    if (experiment_name == "barbule_length"){
        sprintf(variable, "barbule_length_%0.2f", barbule_length);
    }

    // This is uggly and prone to bugs!
    if (experiment_name == "a_inner_barb"){
        sprintf(variable, "a_inner_barb_%0.2f", m_bsdfs[1]->a_inner);
    }

    if (experiment_name == "b_outer_barbule"){
        sprintf(variable, "b_outer_barbule_%0.2f", m_bsdfs[2]->b_outer);
    }

    result = variable;

    printf("variable = %s\n", variable);

    return variable;
}

//Query the Feather BSDF eval function at the incident and outgoing angles provided. For a specific configuration, the BSDF is called N times!
void Feather::latLongQuery(std::vector<Float> phi_i_values, std::vector<Float> theta_i_values, std::vector<Float> phi_o_values, 
                             std::vector<Float> theta_o_values, int N, const char* file_path){
    
    // Pointer to file
    FILE* fp = fopen(file_path, "w");

    //printf("Creating file = %s\n", file_path);
    //printf("Virtual measurements: N = %d\n", N);

    Float barb_eccentricity = 1.0f, barbule_eccentricity = 1.0f;

    fprintf(fp, "%d %d %d %d %d\n", theta_i_values.size(), theta_o_values.size(), phi_i_values.size(), phi_o_values.size(), 3);
    fprintf(fp, "%f %f %f %f %f %f %f\n", virtual_barb_angle, barbule_angle_azimuthal, barbule_angle_longitudinal, barbule_sep, barbule_length, barb_eccentricity, barbule_eccentricity);

    /*
    for (auto theta_o: theta_o_values)
    {
        printf("theta_o = %f\n", theta_o);
    }
    */

    for (auto phi_i: phi_i_values)
    {
        for (auto theta_i: theta_i_values)
        {
            for (auto phi_o: phi_o_values)
            {
                for (auto theta_o: theta_o_values)
                {    
                    //printf("phi_i = %f, theta_i = %f, phi_o = %f, theta_o = %f\n", phi_i, phi_o, theta_i, theta_o);

                    Spectrum measured_radiance = Spectrum(0.0f);

                    for (int n = 0; n < N;n++){
                        Spectrum radiance = measureBRDF(phi_i, theta_i, phi_o, theta_o);

                        //printf("radiance %s\n", radiance.toString().c_str());

                        measured_radiance += radiance;
                        //measured_radiance += measureBRDF(phi_i, theta_i, phi_o, theta_o);
                    }

                    measured_radiance /= N;

                    // Write measurements
    
                    //printf("Writing measurements\n");
                    fprintf(fp, "%f %f %f %f %f %f %f\n", phi_i, theta_i, phi_o, theta_o, measured_radiance[0], measured_radiance[1], measured_radiance[2]);
                }   
            }   
        }
    }

    fclose(fp);// Close file

    //printf("Virtual measurements: Done \n", N);
}

void Feather::configure()
{
    // NB: This function is called after the XML parsing is completed

    m_components.clear();

    unsigned int extraFlags = 0;
    extraFlags |= ESpatiallyVarying;

    //m_components.push_back(ENull | EFrontSide | EBackSide | EGlossyReflection | EAnisotropic | ESpatiallyVarying);
    m_components.push_back(m_bsdfs[0]->getType(0) | extraFlags);
    m_components.push_back(m_bsdfs[1]->getType(0) | extraFlags);
    m_components.push_back(m_bsdfs[2]->getType(0) | extraFlags);
    m_components.push_back(m_bsdfs[3]->getType(0) | extraFlags);
    m_components.push_back(ENull | EFrontSide | EBackSide | extraFlags);

    m_texture = ensureEnergyConservation(m_texture, "feather_texture", 1.0f);
    m_color = ensureEnergyConservation(m_color, "feather_color", 1.0f);
    m_masking_texture = ensureEnergyConservation(m_masking_texture, "feather_masking_texture", 1.0f);

    BSDF::configure();

    if (measure_eval){

        printf("Measuring the BRDF function \n");

        // Feather measurements (angles should be defined outside of this file)
        /*
        std::vector<Float> theta_i_values = {5, 20, 35, 50, 65};
        std::vector<Float> theta_o_values = {0, 15, 30, 45, 60};
        std::vector<Float> phi_i_values = {0, 90};
        std::vector<Float> phi_o_values = {0, 180};
        */

        std::vector<Float> theta_i_values, theta_o_values, phi_i_values, phi_o_values;

        // -----------------Read angles from configuration file-----------------
        //std::string angles_file = "./scenes/measurements/angles.txt";

        FILE* fp = fopen(angles_file.c_str(), "r");

        if (fp == NULL) {
            printf("Failed to open file the angles configuration file: %s\n", angles_file.c_str());
            return;
        }

        char line[1000];
        int n_values;

        fscanf(fp, "%*[^\n]\n");//Skip first line
        

        fscanf(fp, "%s %d\n",line, &n_values);
        //printf("%s\n", line);

        theta_i_values.reserve(n_values);
        theta_i_values.resize(n_values);

        for (int i = 0; i < n_values; i++){
            fscanf(fp, "%f\n", &theta_i_values[i]);
        }

        /*
        printf("Theta i: %d\n", theta_i_values.size());

        for (int i = 0; i < theta_i_values.size(); i++){
            printf("%f ", theta_i_values[i]);
        }

        printf("\n");
        */

        fscanf(fp, "%s %d\n",line, &n_values);
        theta_o_values.reserve(n_values);
        theta_o_values.resize(n_values);

        for (int i = 0; i < n_values; i++){
            fscanf(fp, "%f\n", &theta_o_values[i]);
        }

        /*
        printf("Theta o: %d\n", theta_o_values.size());

        for (int i = 0; i < theta_o_values.size(); i++){
            printf("%f ", theta_o_values[i]);
        }

        printf("\n");
        */

        fscanf(fp, "%s %d\n",line, &n_values);
        phi_i_values.reserve(n_values);
        phi_i_values.resize(n_values);

        for (int i = 0; i < n_values; i++){
            fscanf(fp, "%f\n", &phi_i_values[i]);
        }

        fscanf(fp, "%s %d\n",line, &n_values);
        phi_o_values.reserve(n_values);
        phi_o_values.resize(n_values);

        for (int i = 0; i < n_values; i++){
            fscanf(fp, "%f\n", &phi_o_values[i]);
        }

        //printf("Creating measurements: %s\n", measurements_file.c_str());
        //printf("Ferrero measurements setup: %d\n", ferrero_setup);

        latLongQuery(phi_i_values, theta_i_values, phi_o_values, theta_o_values, n_experiments, measurements_file.c_str());
    }
}

void Feather::serialize(Stream *stream, InstanceManager *manager) const
{
    // The managet is responsible for the serialization and unserialization of objects in the scene
    // Reference textures defined in the scene configuration file can be obtained here.

    BSDF::serialize(stream, manager);
    manager->serialize(stream, m_texture.get());
    manager->serialize(stream, m_color.get());;
    manager->serialize(stream, m_masking_texture.get());

    Assert(m_bsdfs.size() == 4);//Rachis, Barb, Left Barbule and Right Barbule
    manager->serialize(stream, m_bsdfs[0]);
    manager->serialize(stream, m_bsdfs[1]);
    manager->serialize(stream, m_bsdfs[2]);
    manager->serialize(stream, m_bsdfs[3]);
}

Spectrum Feather::getDiffuseReflectance(const Intersection &its) const
{
    return m_color->eval(its);
}

inline Float Feather::getBarbAngle(const Intersection its) const
{
    Vector tex;
    m_texture->eval(its, false).toLinearRGB(tex.x, tex.y, tex.z);
    //m_texture->eval(its, false).toSRGB(tex.x, tex.y, tex.z);

    Float barb_angle = 2.0 * M_PI * tex.x - M_PI; //[0-1] --> [-PI, PI]
    barb_angle = radToDeg(barb_angle);
    barb_angle -= 90;
    
    //printf("barb_angle = %f vs barb_angle2 = %f\n", barb_angle, barb_angle2);

    //barb_angle = barb_angle2;
    //barb_angle = 0.0f;
    
    return barb_angle;
}

// Boolean function that indicates whether a point to evaluate is part of the feather or the background
// NB: In addition, it obtain the color from the geometry and color texture maps
bool Feather::isFeather(const Intersection &its, Float *texture, Float *color) const
{

    // Get the value from the geometry texture if it's available
    m_texture->eval(its, false).toLinearRGB(texture[0], texture[1], texture[2]);
    //m_texture->eval(its, false).toSRGB(texture[0], texture[1], texture[2]);

    // Get the value from the color texture if it's available
    m_color->eval(its, false).toLinearRGB(color[0], color[1], color[2]);
    //m_color->eval(its, false).toSRGB(color[0], color[1], color[2]);

    // printf("Feather geometry texture (%f, %f, %f)\n", texture[0], texture[1], texture[2]);//Debugging

    // Detect the background using the geometry feature, which encodes the geometry of the underlying feather
    if (texture[2] <= 0.1)
        return false;

    return true;
}

// Obtains the alpha matting texture from the geometry texture. The blue channel contains the semantic segmentation of the corresponding
inline Spectrum getOpacity(Float *geometry_texture){

    //return Spectrum( (geometry_texture[2] <= 0.1) ? 0.0f : 1.0f );

    return Spectrum( (geometry_texture[2] <= 0.9) ? geometry_texture[2] : 1.0f );
}

// Spectrum to RGB safe transformation
// NB: Taking of numerical issues: (1) NAN --> Black (2) Clipping to [0.0, max_value]. For eval is usually 1000, while for sample_max is usually 2.0
inline Spectrum toSpectrum(Spectrum value, Float max_value)
{
    Spectrum new_value;

    if (std::isnan(value[0]) || std::isnan(value[1]) || std::isnan(value[2])){

        printf("Feather BSDF(NAN cases): %f\n", max_value);

        return Spectrum(0.0f);
    }

    if (std::isinf(value[0]) || std::isinf(value[1]) || std::isinf(value[2])){
        
        printf("Feather BSDF (infinity cases): %f\n", max_value);

        return Spectrum(0.0f);
    }
    
    //Clamp weights [0, max_value]
    for (int i = 0; i < 3; i++){
        new_value[i] = std::min(std::max(Float(0.0f), value[i]), max_value);
    }

    return new_value;
}

inline Spectrum Feather::debuggingGeometry(const BSDFSamplingRecord &bRec, EMeasure measure) const
{

    // printf("Debugging Geometry\n");

    Float texture_rgb[3] = {0.0f, 0.0f, 0.0f};
    Float color_rgb[3] = {0.0f, 0.0f, 0.0f};
    Spectrum rachis_color;

    if (measure != ESolidAngle || Frame::cosTheta(bRec.wi) <= 0 || Frame::cosTheta(bRec.wo) <= 0)
    {
        return Spectrum(0.0);
    }

    // Check if the evaluation point (underlying surface interaction) is inside the feather --> Obtain the values from the corresponding feather textures (geometry and color)
    if (!isFeather(bRec.its, texture_rgb, color_rgb))
    {
        return Spectrum(((bRec.typeMask & ENull) && measure == EDiscrete) ? 1.0f : 0.0f);
    }

    // Debugging geometry texture

    // printf("debugging_geometry\n");
    return Spectrum(texture_rgb); // Geometry texture
    // return Spectrum(color_rgb);//Color texture
}

// Dario Lanza
Spectrum Feather::evalTransmittanceEvent(const BSDFSamplingRecord &bRec, EMeasure measure) const
{
    if (measure == EDiscrete && std::abs(1-dot(bRec.wi, -bRec.wo)) < DeltaEpsilon){
        return Spectrum(1.0f);
    }
    
    return Spectrum(0.0f);
}

// Evaluate the BSDF f(wi, wo), the blue channel of the feather texture (texture_rgb) is used to discriminate the corresponding material to be evaluated
Spectrum Feather::featherEval(const BSDFSamplingRecord &bRec, EMeasure measure, Float texture_mask) const
{

    //int thread_id = omp_get_thread_num();

    //printf("thread_id = %d\n", thread_id);

    if (texture_mask <= shaft_mask)
    { // Rachis

        Spectrum rachis_spectrum = m_bsdfs[0]->eval(bRec, measure);
        //return rachis_spectrum;
        return toSpectrum(rachis_spectrum, eval_max);
    }

    //-------------------------------Coordinate system transformations -----------------------------------
    BSDFSamplingRecord bRec_barb(bRec), bRec_left_barbule(bRec), bRec_right_barbule(bRec);
        
    Float barb_angle = getBarbAngle(bRec.its);
    Frame barb_f,left_barbule_f,right_barbule_f;
    setLocalFrames(barb_f, left_barbule_f,right_barbule_f, barb_angle);
    
    Vector temp_wi = bRec.wi;
    Vector temp_wo = bRec.wo;
    bRec_left_barbule.wi = left_barbule_f.toLocal(temp_wi);
    bRec_left_barbule.wo = left_barbule_f.toLocal(temp_wo);
    
    temp_wi = bRec.wi;
    temp_wo = bRec.wo;
    bRec_right_barbule.wi = right_barbule_f.toLocal(temp_wi);
    bRec_right_barbule.wo = right_barbule_f.toLocal(temp_wo);

    temp_wi = bRec.wi;
    temp_wo = bRec.wo;
    bRec_barb.wi = barb_f.toLocal(temp_wi);
    bRec_barb.wo = barb_f.toLocal(temp_wo);

    Spectrum result(0.0f);

    //printf("Masking enabled? %d\n", use_masking);

    if (use_masking){

        Feather_masking::Masking masking = getMasking(barb_f, left_barbule_f,right_barbule_f,bRec.its.wi);

        //printf("Barb range size: %d\n", masking.barb_ranges().size());

        int count = 0;

        for (auto barb_range : masking.barb_ranges()){
            
            //mutex_feather.lock();
            //m_bsdfs[1]->hmin = barb_range.hmin(); m_bsdfs[1]->hmax = barb_range.hmax();

            /*
            if (barb_range.hmin() > barb_range.hmax()){
                printf("Wrong h order for Barb (eval): [%f, %f]. Barb range size = %d\n", barb_range.hmin(), barb_range.hmax(), masking.barb_ranges().size());
            
                Spectrum test;
                test[2] = 1.0;

                return test;
            }else{
                return Spectrum();
            }
            */

            /*
            if (barb_range.hmin() > barb_range.hmax()){
                printf("Wrong h order for Barb (eval): [%f, %f]. Barb range size = %d (case %d), angles (phi_i = %f, theta_i = %f)\n", barb_range.hmin(), barb_range.hmax(), masking.barb_ranges().size(), count);
            }
            */
            

            //printf("Barb H = [%f, %f]\n", m_bsdfs[1]->hmin, m_bsdfs[1]->hmax);

            if (debugging_bsdf){
                result += barb_color * barb_range.rate();
            }else{
                result += m_bsdfs[1]->eval(bRec_barb, measure) * barb_range.rate();
            }

            //mutex_feather.unlock();

            count++;
        }

        /*
        if (masking.left_range().hmin() > masking.left_range().hmax()){
            printf("Wrong h order for left barbule (eval): [%f, %f]\n", masking.left_range().hmin(), masking.left_range().hmax());
        }

        if (masking.right_range().hmin() > masking.right_range().hmax()){
            printf("Wrong h order for right barbule (eval): [%f, %f]\n", masking.right_range().hmin(), masking.right_range().hmax());
        }
        */

        //Barbule H range can be update directly.
        //m_bsdfs[2]->hmin = masking.left_range().hmin(); m_bsdfs[2]->hmax = masking.left_range().hmax();
        //m_bsdfs[3]->hmin = masking.right_range().hmin(); m_bsdfs[3]->hmax = masking.right_range().hmax();

        if (debugging_bsdf){
            result += left_barbule_color * masking.left_range().rate();
            result += right_barbule_color * masking.right_range().rate();
        }else{
            result += m_bsdfs[2]->eval(bRec_left_barbule, measure) * masking.left_range().rate();
            result += m_bsdfs[3]->eval(bRec_right_barbule, measure) * masking.right_range().rate();
        }

        //result += m_bsdfs[2]->eval(bRec_left_barbule, measure);
        //result += m_bsdfs[3]->eval(bRec_right_barbule, measure);

    }else{

        //To combine the barb and barbules BSDF in a logical manner, we require the the masking weights. 
        //For this reason, we only use the barb term if we do not enable the masking

        result += m_bsdfs[1]->eval(bRec_barb, measure);
    }

    /*
    if (bRec.typeMask == EAll){
        printf("bRec = %s\n", bRec.toString().c_str());
        printf("measure = %d\n", measure);
        printf("result = %s\n", result.toString().c_str());
        printf("barb_angle = %f\n", barb_angle);
    }
    */

    //return result;
    return toSpectrum(result, eval_max);
}

// Evaluate the particular material depending of the geometry texture mask
Spectrum Feather::eval(const BSDFSamplingRecord &bRec, EMeasure measure) const
{
    /*
    printf("eval: its.geoFrame.s %s\n", bRec.its.geoFrame.s.toString().c_str());
    printf("eval: its.geoFrame.t %s\n", bRec.its.geoFrame.t.toString().c_str());
    printf("eval: its.geoFrame.n %s\n", bRec.its.geoFrame.n.toString().c_str());
    printf("eval: its.shFrame.s %s\n", bRec.its.shFrame.s.toString().c_str());
    printf("eval: its.shFrame.t %s\n", bRec.its.shFrame.t.toString().c_str());
    printf("eval: its.shFrame.n %s\n", bRec.its.shFrame.n.toString().c_str());
    */

    /*
    if (bRec.typeMask == EAll){
        printf("bRec = %s\n", bRec.toString().c_str());
        printf("measure = %d\n", measure);
    }
    */

    //printf("typeMask = %d\n", bRec.typeMask);

    //typeMask options:
    /*
    printf("ENull = %d\n", ENull);
    printf("EDiffuseReflection = %d\n", EDiffuseReflection);
    printf("EDiffuseTransmission = %d\n", EDiffuseTransmission);
    printf("EGlossyReflection = %d\n", EGlossyReflection);
    printf("EGlossyTransmission = %d\n", EGlossyTransmission);
    printf("EDeltaReflection = %d\n", EDeltaReflection);
    printf("EDeltaTransmission = %d\n", EDeltaTransmission);
    printf("EDelta1DReflection = %d\n", EDelta1DReflection);
    printf("EDelta1DTransmission = %d\n", EDelta1DTransmission);
    printf("EAnisotropic = %d\n", EAnisotropic);
    printf("ESpatiallyVarying = %d\n", ESpatiallyVarying);
    printf("ENonSymmetric = %d\n", ENonSymmetric);
    printf("EFrontSide = %d\n", EFrontSide);
    printf("EBackSide = %d\n", EBackSide);
    printf("EUsesSampler = %d\n", EUsesSampler);
    */

    //ETypeCombinations
    /*
    printf("EReflection = %d\n", EReflection);
    printf("ETransmission = %d\n", ETransmission);
    printf("EDiffuse = %d\n", EDiffuse);
    printf("EGlossy = %d\n", EGlossy);
    printf("ESmooth = %d\n", ESmooth);
    printf("EDelta = %d\n", EDelta);
    printf("EDelta1D = %d\n", EDelta1D);
    printf("EAll = %d\n", EAll);
    */

    //printf("bRec.eta %f\n", bRec.eta);

    /*
    printf("bRec.mode = %d\n", bRec.mode);
    printf("bRec.component = %d\n", bRec.component);
    printf("bRec.sampledType = %d\n", bRec.sampledType);
    printf("bRec.sampledComponent = %d\n", bRec.sampledComponent);

    printf("ERadiance = %d\n", ERadiance);
    printf("EImportance = %d\n", EImportance);
    printf("ETransportModes = %d\n", ETransportModes);
    */

    //EMeasure
    /*
    printf("EInvalidMeasure = %d\n", EInvalidMeasure);
    printf("ESolidAngle = %d\n", ESolidAngle);
    printf("ELength = %d\n", ELength);
    printf("EArea = %d\n", EArea);
    printf("EDiscrete = %d\n", EDiscrete);
    */

    //Feather texture geometry: semantic segmentation of feather's hierarchical structure
    Float texture_rgb[3] = {0.0f, 0.0f, 0.0f};
    Float color_rgb[3] = {0.0f, 0.0f, 0.0f};

    isFeather(bRec.its, texture_rgb, color_rgb);    
    
    //Spectrum opacity = Spectrum(texture_rgb);
    //Spectrum opacity = getOpacity(texture_rgb);//Old approach
    Spectrum opacity = m_masking_texture->eval(bRec.its,true);// New version Dario

    if (measure == ESolidAngle){
        //return m_bsdfs[1]->eval(bRec, ESolidAngle) * opacity;
        //return featherEval(bRec, ESolidAngle, texture_rgb[2]) * opacity;

        if (  opacity.getLuminance() < 0.05){
            return Spectrum(0.0f);

        }else{
            return featherEval(bRec, ESolidAngle, texture_rgb[2]) * opacity;

        }
    }
    else if (measure == EDiscrete && std::abs(1-dot(bRec.wi, -bRec.wo)) < DeltaEpsilon)
        return Spectrum(1.0f) - opacity;
    else
        return Spectrum(0.0f);
}

// Dario Lanza
Float Feather::pdfTransmittanceEvent(const BSDFSamplingRecord &bRec, EMeasure measure) const {
    bool sampleTransmission = bRec.typeMask & ENull && (bRec.component == -1 );

    //bool sampleNested = bRec.component == -1 || bRec.component < getComponentCount()-1;
    if(sampleTransmission){
        return 1.0f;
    }

    return 0.0f;
    
    /*Float prob = m_opacity->eval(bRec.its).getLuminance();
    if (measure == ESolidAngle) {
        if (!sampleNested)
            return 0.0f;
        Float result = m_nestedBSDF->pdf(bRec, ESolidAngle);
        if (sampleTransmission)
            result *= prob;
        return result;
    } else
    if (measure == EDiscrete && std::abs(1-dot(bRec.wi, -bRec.wo)) < DeltaEpsilon) {
        if (!sampleTransmission)
            return 0.0f;
        if (!sampleNested)
            return 1.0f;
        else
            return 1-prob;
    } else {
        return 0.0f;
    }*/
}

// Compute the pdf of f(wi, wo), the blue channel of the feather texture (texture_rgb) is used to discriminate the corresponding material to be evaluated
Float Feather::featherPdf(const BSDFSamplingRecord &bRec, EMeasure measure, Float texture_mask) const{
    
    // Evaluate the particular material depending of the geometry texture mask
    
    if (texture_mask <= shaft_mask)
    { // Rachis

        Float rachis_pdf = m_bsdfs[0]->pdf(bRec, measure);
        //return rachis_pdf;
        return std::max(Float(0.0f), rachis_pdf);
    }

    //-------------------------------Coordinate system transformations -----------------------------------
    BSDFSamplingRecord bRec_barb(bRec), bRec_left_barbule(bRec), bRec_right_barbule(bRec);
        
    Float barb_angle = getBarbAngle(bRec.its);
    Frame barb_f,left_barbule_f,right_barbule_f;
    setLocalFrames(barb_f, left_barbule_f,right_barbule_f, barb_angle);
    
    Vector temp_wi = bRec.wi;
    Vector temp_wo = bRec.wo;
    bRec_left_barbule.wi = left_barbule_f.toLocal(temp_wi);
    bRec_left_barbule.wo = left_barbule_f.toLocal(temp_wo);
    
    temp_wi = bRec.wi;
    temp_wo = bRec.wo;
    bRec_right_barbule.wi = right_barbule_f.toLocal(temp_wi);
    bRec_right_barbule.wo = right_barbule_f.toLocal(temp_wo);

    temp_wi = bRec.wi;
    temp_wo = bRec.wo;
    bRec_barb.wi = barb_f.toLocal(temp_wi);
    bRec_barb.wo = barb_f.toLocal(temp_wo);

    Float result = 0.0f;

    if (use_masking){

        Feather_masking::Masking masking = getMasking(barb_f, left_barbule_f,right_barbule_f,bRec.its.wi);

        //printf("Barb range size: %d\n", masking.barb_ranges().size());

        for (auto barb_range : masking.barb_ranges()){
            
            //mutex_feather.lock();
            //m_bsdfs[1]->hmin = barb_range.hmin(); m_bsdfs[1]->hmax = barb_range.hmax();

            //printf("Barb H = [%f, %f]\n", m_bsdfs[1]->hmin, m_bsdfs[1]->hmax);

            /*
            if (barb_range.hmin() > barb_range.hmax()){
                printf("Wrong h order for Barb (pdf): [%f, %f]. Barb range size = %d\n", barb_range.hmin(), barb_range.hmax(), masking.barb_ranges().size());
            }
            */
            
            result += m_bsdfs[1]->pdf(bRec_barb, measure) * barb_range.rate();
            //mutex_feather.unlock();
        }

        //Barbule H range can be update directly.
        //m_bsdfs[2]->hmin = masking.left_range().hmin(); m_bsdfs[2]->hmax = masking.left_range().hmax();
        //m_bsdfs[3]->hmin = masking.right_range().hmin(); m_bsdfs[3]->hmax = masking.right_range().hmax();

        if (masking.left_range().hmin() > masking.left_range().hmax()){
            printf("Wrong h order for left barbule (pdf): [%f, %f]\n", masking.left_range().hmin(), masking.left_range().hmax());
        }

         if (masking.right_range().hmin() > masking.right_range().hmax()){
            printf("Wrong h order for right barbule (pdf): [%f, %f]\n", masking.right_range().hmin(), masking.right_range().hmax());
        }

        result += m_bsdfs[2]->pdf(bRec_left_barbule, measure) * masking.left_range().rate();
        result += m_bsdfs[3]->pdf(bRec_right_barbule, measure) * masking.right_range().rate();


        //----------------------Probability density functions validations
        // It happens mostly for the first Rachis BSDF
        
        /*
        int K = 4;
        for (int k = 1; k < K;k++){
            
            Float pdf_test = m_bsdfs[k]->pdf(bRec_left_barbule, measure); 

            if (pdf_test > 1.0){
                printf("The pdf for the BSDF [%d] = %f, is larger than 1\n", k, pdf_test);
            }
        }
        */

    }else{

        //To combine the barb and barbules BSDF in a logical manner, we require the the masking weights. 
        //For this reason, we only use the barb term if we do not enable the masking

        result += m_bsdfs[1]->pdf(bRec_barb, measure);
    }

    //return result;
    return std::max(Float(0.0f), result);
}

// Compute the probability of sampling wo given wi (BSDF f(wi, wo))
// Similar to eval but using the corresponding pdf functions
Float Feather::pdf(const BSDFSamplingRecord &bRec, EMeasure measure) const
{        
    Float texture_rgb[3] = {0.0f, 0.0f, 0.0f};
    Float color_rgb[3] = {0.0f, 0.0f, 0.0f};

    isFeather(bRec.its, texture_rgb, color_rgb);

    bool sampleTransmission = bRec.typeMask & ENull
            && (bRec.component == -1 || bRec.component == getComponentCount()-1);
    bool sampleNested = bRec.component == -1 || bRec.component < getComponentCount()-1;

    //Float prob = texture_rgb[2];

    //Spectrum opacity = getOpacity(texture_rgb);//Old approach
    Spectrum opacity = m_masking_texture->eval(bRec.its,true);//New approach Dario
    Float prob = opacity.getLuminance();
    
    if (measure == ESolidAngle) {
        if (!sampleNested)
            return 0.0f;
        
        //Float result = m_bsdfs[1]->pdf(bRec, ESolidAngle);
        Float result = featherPdf(bRec, ESolidAngle, texture_rgb[2]);
        
        if (sampleTransmission)
            result *= prob;

        return result;

    } else if (measure == EDiscrete && std::abs(1-dot(bRec.wi, -bRec.wo)) < DeltaEpsilon) {
        if (!sampleTransmission)
            return 0.0f;
        if (!sampleNested)
            return 1.0f;
        else
            return 1-prob;
    } else {
        return 0.0f;
    }
}

// Importance sampling weight --> Sample the BSDF and divided by the pdf of the sample
Spectrum Feather::sample(BSDFSamplingRecord &bRec, const Point2 &sample) const
{
    printf("This function is never called\n");

    Float pdf;
    return Feather::sample(bRec, pdf, sample);
}

Spectrum Feather::sampleTransmittanceEvent(BSDFSamplingRecord &bRec, Float &pdf, const Point2 &_sample) const
{
    Point2 sample(_sample);
    Spectrum result(0.0f);
  
    // Spectrum opacity = Spectrum(0.0f);//m_opacity->eval(bRec.its);
    // Float prob = opacity.getLuminance();

    /*
    if (pdf <= 0){
        printf("sampleTransmittanceEvent: Negative o 0 pdf (transmittance). This should be impossible\n");
    }
    */

    bRec.wo = -bRec.wi;
    bRec.eta = 1.0f;
    bRec.sampledComponent = 3;
    
    bRec.sampledType = ENull;

    //printf("sampleTransmittanceEvent pdf = %f\n", pdf);

    result = (Spectrum(1.0f)) / pdf;

    return result;
}

Spectrum Feather::featherSample(BSDFSamplingRecord &bRec, Float &pdf, const Point2 &_sample, Float texture_mask) const{

    if (texture_mask <= shaft_mask)
    { // Rachis
        Point2 m_sample(_sample);

        Spectrum result = m_bsdfs[0]->sample(bRec, pdf, m_sample);

        return toSpectrum(result, sample_max);
    }

    // Coordinate system
    Float barb_angle = getBarbAngle(bRec.its);
    BSDFSamplingRecord bRec_barb(bRec), bRec_left_barbule(bRec), bRec_right_barbule(bRec);    
    Frame barb_f,left_barbule_f,right_barbule_f;

    Float left_barbule_transparency, right_barbule_transparency;
    flatland::Vec2 barb_dir, left_barbule_dir, right_barbule_dir;

   // printf("bRec shFram s %s t %s n %s ",bRec.its.geoFrame.s.toString().c_str(),bRec.its.geoFrame.t.toString().c_str(),bRec.its.geoFrame.n.toString().c_str());
    Vector local_x(1,0,0);
    Vector local_y(0,1,0);
    Vector local_z(0,0,1);
    
    local_x = bRec.its.geoFrame.toWorld(local_x);
    local_y = bRec.its.geoFrame.toWorld(local_y);
    local_z = bRec.its.geoFrame.toWorld(local_z);

    //printf("local_x %s y %s z %s \n",local_x.toString().c_str(),local_y.toString().c_str(),local_z.toString().c_str());
    setLocalFrames(barb_f, left_barbule_f,right_barbule_f, barb_angle);

    //printf("\n");
    bRec_left_barbule.wi = left_barbule_f.toLocal(bRec.wi); 
    bRec_left_barbule.wo = left_barbule_f.toLocal(bRec.wo);
    
    bRec_right_barbule.wi = right_barbule_f.toLocal(bRec.wi);
    bRec_right_barbule.wo = right_barbule_f.toLocal(bRec.wo);
    
    bRec_barb.wi = barb_f.toLocal(bRec.wi); 
    bRec_barb.wo = barb_f.toLocal(bRec.wo);

    if (!use_masking){
        Point2 m_sample(_sample);

        Spectrum result = m_bsdfs[1]->sample(bRec_barb, pdf, m_sample);

        //Dario: Transforming back to world coordinate system (we do the same for barbule left and barbule right)
        bRec.wo = barb_f.toWorld(bRec_barb.wo);
        bRec.wi = barb_f.toWorld(bRec_barb.wi);
        bRec.eta = 1.0f;
        bRec.sampledType = EGlossyReflection;
        bRec.sampledComponent = 0;

        return toSpectrum(result, sample_max);
    }

    //---------------Masking using weights and H range
    
    Feather_masking::Masking masking = getMasking(barb_f, left_barbule_f,right_barbule_f,bRec.its.wi);
    int barb_range_N = masking.barb_ranges().size();//Number of barb interval events

    DiscreteDistribution m_pdf;
    std::vector<Float> m_weights; 

    m_pdf = DiscreteDistribution(barb_range_N+3);//Barb range events + Left barbule event + Right Barbule event + Transmission event

    for (auto barb_range : masking.barb_ranges()){
        m_pdf.append(barb_range.rate());
        m_weights.push_back(barb_range.rate());
    }

    m_pdf.append(masking.left_range().rate());
    m_weights.push_back(masking.left_range().rate());

    m_pdf.append(masking.right_range().rate());
    m_weights.push_back(masking.right_range().rate());

    m_pdf.append(masking.transmittance_rate());
    m_weights.push_back(masking.transmittance_rate());

    m_pdf.normalize();
    Spectrum result(0.0f);

    //Choose a component based on the normalized weights
    Point2 m_sample(_sample);

    size_t m_idx = m_pdf.sampleReuse(m_sample.x); // entry == m_idx
    bool transmittance_event = true;

    if (m_idx < barb_range_N){

        //printf(" Barb ");

        //This is uggly I know, but for some strange reason we are using a list for the Barb Range
        auto barb_range_it = masking.barb_ranges().begin();
        std::advance(barb_range_it, m_idx);

        /*
        if ((*barb_range_it).hmin() > (*barb_range_it).hmax()){
            printf("Wrong h order for Barb (sample): [%f, %f]. Barb range size = %d\n", (*barb_range_it).hmin(), (*barb_range_it).hmax(), masking.barb_ranges().size());
        
            Spectrum test;
            test[2] = 1.0;

            return test;
        }else{
            return Spectrum();
        }
        */

        //mutex_feather.lock();
        if ((*barb_range_it).hmin() > (*barb_range_it).hmax()){
            printf("Wrong h order for Barb (sample): [%f, %f]. Barb range size = %d\n", (*barb_range_it).hmin(), (*barb_range_it).hmax(), masking.barb_ranges().size());
        }

        //m_bsdfs[1]->hmin = (*barb_range_it).hmin(); m_bsdfs[1]->hmax = (*barb_range_it).hmax();
        
        if (pdf < 0){
            printf("feather BSDF (sample, barb): pdf = %f\n", pdf);
        }

        if(debugging_bsdf){
            result = Spectrum(barb_color);
        }else{
            result = m_bsdfs[1]->sample(bRec_barb, pdf, m_sample);
        }

        //mutex_feather.unlock();

        bRec.eta = 1.0f;
        bRec.wo = barb_f.toWorld(bRec_barb.wo);
        bRec.sampledType = EGlossyReflection;
        //bRec.sampledType = EDiffuseTransmission;
        bRec.sampledComponent = 0;
        transmittance_event = false;

        if (pdf < 0){
            printf("feather BSDF (sample): barb pdf is negative\n");
        }
    }
            
    if (m_idx == barb_range_N){
        //printf(" Barbules LEFT ");
        //printf(" %f %f %f ",bRec_left_barbule.wi[0],bRec_left_barbule[1],bRec_left_barbule[2]);
        
        if (masking.left_range().hmin() > masking.left_range().hmax()){
            printf("Wrong h order for left barbule (sample): [%f, %f]\n", masking.left_range().hmin(), masking.left_range().hmax());
        }

        m_bsdfs[2]->hmin = masking.left_range().hmin(); m_bsdfs[2]->hmax = masking.left_range().hmax();
        
        if(debugging_bsdf){
            result = Spectrum(left_barbule_color);
        }else{
            result = m_bsdfs[2]->sample(bRec_left_barbule, pdf, m_sample);
        }

        bRec.eta = 1.0f;
        bRec.wo = left_barbule_f.toWorld(bRec_left_barbule.wo);
        // bRec.wo = bRec_left_barbule.wo;

        bRec.sampledComponent = 1; // --> change this for the eval/sampling
        bRec.sampledType = EGlossyReflection;
        // bRec.sampledType = EDiffuseTransmission;
        // bRec.sampledType = EBarbuleLeft;

        transmittance_event = false;

        if (pdf < 0){
            printf("feather BSDF (sample): barbule left pdf is negative\n");
        }
    }

    if (m_idx == (barb_range_N + 1)){
        // printf(" Barbules RIGHT ");
        
        if (masking.right_range().hmin() > masking.right_range().hmax()){
            printf("Wrong h order for right barbule (sample): [%f, %f]\n", masking.right_range().hmin(), masking.right_range().hmax());
        }

        m_bsdfs[3]->hmin = masking.right_range().hmin(); m_bsdfs[3]->hmax = masking.right_range().hmax();
        
        if(debugging_bsdf){
            result = Spectrum(right_barbule_color);
        }else{
            result = m_bsdfs[3]->sample(bRec_right_barbule, pdf, m_sample);
        }

        bRec.eta = 1.0f;
        bRec.wo = right_barbule_f.toWorld(bRec_right_barbule.wo);
        bRec.sampledComponent = 2;
        bRec.sampledType = EGlossyReflection;
        //bRec.sampledType = EDiffuseTransmission;

        transmittance_event = false;

        if (pdf < 0){
            printf("feather BSDF (pdf): barbule right pdf is negative\n");
        }
    }
    
    if(transmittance_event){
        // printf(" Transmittance ");
        result = sampleTransmittanceEvent(bRec, pdf, m_sample);
    }
    
    if (result.isZero()) 
        return result;
       
    //printf("Feather BSDF (sample): pdf=%f\n", pdf);

    result *= m_weights[m_idx] * pdf;
    pdf *= m_pdf[m_idx];

    if (pdf <= 0){
        printf("Feather BSDF (sample): Negative o 0 pdf, pdf = %f, m_pdf = %f. This should be impossible\n", pdf, m_pdf[m_idx]);
        pdf = 0.0f;
        return Spectrum(0.0f);
    }

    //printf("Sampling clamping maximum %f\n", sample_max);

    Spectrum weight = result / pdf;

    return toSpectrum(weight, sample_max);//Clamp weights [0, sample_max]. Linqi uses sample_max = 2
}

// Importance sampling weight --> Sample the BSDF and divided by the pdf of the sample
// NB: Similar to the other method but the pdf is provided
Spectrum Feather::sample(BSDFSamplingRecord &bRec, Float &pdf, const Point2 &_sample) const
{

    //printf("Feather (sample)\n");

    Float texture_rgb[3] = {0.0f, 0.0f, 0.0f};
    Float color_rgb[3] = {0.0f, 0.0f, 0.0f};

    isFeather(bRec.its, texture_rgb, color_rgb);    

    Point2 sample(_sample);
    Spectrum result(0.0f);
    
    //Spectrum opacity = Spectrum(texture_rgb);
    //Float prob = texture_rgb[2];

    //Spectrum opacity = getOpacity(texture_rgb);//Old approach
    Spectrum opacity = m_masking_texture->eval(bRec.its,true);//New approach Dario
    Float prob = opacity.getLuminance();

    //printf("Feather BSDF (sample): initial pdf = %f, prob = %f\n", pdf, prob);

    bool sampleTransmission = bRec.typeMask & ENull
        && (bRec.component == -1 || bRec.component == getComponentCount()-1);
    bool sampleNested = bRec.component == -1 || bRec.component < getComponentCount()-1;

    if (sampleTransmission && sampleNested) {
        if (sample.x < prob) {
            Float invProb = 1.0f / prob;
            sample.x *= invProb;
            //return m_bsdfs[1]->sample(bRec, sample) * opacity * invProb;
            result = featherSample(bRec, pdf, sample, texture_rgb[2]) * opacity * invProb;
            pdf *= prob;

        } else {
            bRec.wo = -bRec.wi;
            bRec.eta = 1.0f;
            bRec.sampledComponent = getComponentCount()-1;
            bRec.sampledType = ENull;
            pdf = 1.0f - prob;
            result = (Spectrum(1.0f) - opacity) / (1-prob);
        }
    } else if (sampleTransmission) {
        bRec.wo = -bRec.wi;
        bRec.eta = 1.0f;
        bRec.sampledComponent = getComponentCount()-1;
        bRec.sampledType = ENull;
        pdf = 1.0f;
        result = Spectrum(1.0f) - opacity;
    } else if (sampleNested) {
        //return m_bsdfs[1]->sample(bRec, sample) * opacity;

       result = featherSample(bRec, pdf, sample, texture_rgb[2]) * opacity;
    }

    return result;
}

//Adolfo transformation suggestions
void Feather::getVectorMeasurements(Float phi_i, Float theta_i, Float phi_o, Float theta_o, Vector &wi, Vector &wo){

    if(phi_i == 0.0f && phi_o == 0.0f){
        phi_i = degToRad(phi_i);theta_i = degToRad(theta_i); phi_o = degToRad(phi_o);theta_o = degToRad(theta_o); 

        wi = Vector3f(-sin(theta_i),0.0f, cos(theta_i));
        wo = Vector3f(-sin(theta_o),0.0f, cos(theta_o));
    }

    if(phi_i == 0.0f && phi_o == 180.0f){
        phi_i = degToRad(phi_i);theta_i = degToRad(theta_i); phi_o = degToRad(phi_o);theta_o = degToRad(theta_o); 

        wi = Vector3f(-sin(theta_i), 0.0f,cos(theta_i));
        wo = Vector3f(sin(theta_o), 0.0f, cos(theta_o));
    }

    if(phi_i == 90.0f && phi_o == 0.0f){
        phi_i = degToRad(phi_i);theta_i = degToRad(theta_i); phi_o = degToRad(phi_o);theta_o = degToRad(theta_o); 

        wi = Vector3f(0.0f,-sin(theta_i), cos(theta_i));
        wo = Vector3f(0.0f,-sin(theta_o), cos(theta_o));
    }

    if(phi_i == 90.0f && phi_o == 180.0f){
        phi_i = degToRad(phi_i);theta_i = degToRad(theta_i); phi_o = degToRad(phi_o);theta_o = degToRad(theta_o); 

        wi = Vector3f(0.0f, -sin(theta_i),cos(theta_i));
        wo = Vector3f(0.0f, sin(theta_o),cos(theta_o));
    }
}

// Given the incident and outgoing angles, this function record the reflectance of our feather BSDF.
// To Do: Fiber BCSDF needs to be wavelength dependent. Current implementation is RGB-based.
Spectrum Feather::measureBRDF(Float phi_i, Float theta_i, Float phi_o, Float theta_o){

    //printf("Measuring BRDF: phi_i = %f, theta_i = %f, phi_o = %f, theta_o = %f, hit_rachis = %d, barb_angle = %f\n", phi_i, theta_i, phi_o, 
    //      theta_o, hit_rachis, virtual_barb_angle);

    //return Spectrum(0.0f);

    // In a general case, there are three possible scenearios: hit background, hit rachis and hit barb. For measurements, we probably only care
    // about the barb case

    // Create Mitsuba BSDF query: BSDFSamplingRecord and measure 
    Vector wi, wo;
    mitsuba::Intersection its;

    /*
    printf("measureBRDF: its.geoFrame.s %s\n", its.geoFrame.s.toString().c_str());
    printf("measureBRDF: its.geoFrame.t %s\n", its.geoFrame.t.toString().c_str());
    printf("measureBRDF: its.geoFrame.n %s\n", its.geoFrame.n.toString().c_str());
    printf("measureBRDF: its.shFrame.s %s\n", its.shFrame.s.toString().c_str());
    printf("measureBRDF: its.shFrame.t %s\n", its.shFrame.t.toString().c_str());
    printf("measureBRDF: its.shFrame.n %s\n", its.shFrame.n.toString().c_str());
    */

    its.geoFrame.s = Vector(1.0, 0.0, 0.0);
    its.geoFrame.t = Vector(0.0, 1.0, 0.0);
    its.geoFrame.n = Vector(0.0, 0.0, 1.0);
    its.shFrame.s = Vector(1.0, 0.0, 0.0);
    its.shFrame.t = Vector(0.0, 1.0, 0.0);
    its.shFrame.n = Vector(0.0, 0.0, 1.0);

    // Common spherical coordinates to cartessian coordinates transformation
    
    //printf("ferrero setup: %d\n", ferrero_setup);

    if (ferrero_setup){
        phi_o = phi_i + phi_o;// relative to phi_i! 
    }
    
    phi_i = degToRad(phi_i);theta_i = degToRad(theta_i); phi_o = degToRad(phi_o);theta_o = degToRad(theta_o); 

    // Mitsuba spherical coordinate system
    //wi = Vector3f(cos(phi_i) * sin(theta_i), sin(phi_i) * sin(theta_i), cos(theta_i));
    //wo = Vector3f(cos(phi_o) * sin(theta_o), sin(phi_o) * sin(theta_o), cos(theta_o));

    // Latlong plots
    wo = Vector3f(cos(phi_i) * sin(theta_i), sin(phi_i) * sin(theta_i), cos(theta_i));
    wi = Vector3f(cos(phi_o) * sin(theta_o), sin(phi_o) * sin(theta_o), cos(theta_o));

    //getVectorMeasurements(phi_i, theta_i, phi_o, theta_o, wi, wo);//Adolfo's transformations

    BSDFSamplingRecord bRec = BSDFSamplingRecord(its, wi, wo, ERadiance);
    EMeasure measure = ESolidAngle;

    //printf("bRec.its.wi = %s\nwi = %s\n", bRec.its.wi.toString().c_str(), wi.toString().c_str());

    //printf("phi_i = %f, phi_o = %f, theta_i = %f, theta_o = %f, wi = %s, wo = %s\n", phi_i, phi_o, theta_i, theta_o, wi.toString().c_str(), wo.toString().c_str());

    //printf("bRec = %s\n", bRec.toString().c_str());

    // Rachis
    if (hit_rachis)
    {
        //printf("Measuring material directly\n");

        Spectrum rachis_spectrum = m_bsdfs[0]->eval(bRec, measure);
        return rachis_spectrum;
    }

    //-------------------------------Coordinate system transformations -----------------------------------
    //printf("Setting the coordinate system: barb_angle = %f\n", virtual_barb_angle);

    BSDFSamplingRecord bRec_barb(bRec), bRec_left_barbule(bRec), bRec_right_barbule(bRec);
        
    //Float barb_angle = virtual_barb_angle;
    Float barb_angle = virtual_barb_angle - 90;
    Frame barb_f,left_barbule_f,right_barbule_f;
    setLocalFrames(barb_f, left_barbule_f,right_barbule_f, barb_angle);
    
    Vector temp_wi = bRec.wi;
    Vector temp_wo = bRec.wo;
    bRec_left_barbule.wi = left_barbule_f.toLocal(temp_wi);
    bRec_left_barbule.wo = left_barbule_f.toLocal(temp_wo);
    
    temp_wi = bRec.wi;
    temp_wo = bRec.wo;
    bRec_right_barbule.wi = right_barbule_f.toLocal(temp_wi);
    bRec_right_barbule.wo = right_barbule_f.toLocal(temp_wo);

    temp_wi = bRec.wi;
    temp_wo = bRec.wo;
    bRec_barb.wi = barb_f.toLocal(temp_wi);
    bRec_barb.wo = barb_f.toLocal(temp_wo);

    //printf("Barb: phi_i = %f, phi_o = %f, theta_i = %f, theta_o = %f, wi = %s, wo = %s\n", phi_i, phi_o, theta_i, theta_o, bRec_barb.wi.toString().c_str(), bRec_barb.wo.toString().c_str());
    //printf("Left barbule: phi_i = %f, phi_o = %f, theta_i = %f, theta_o = %f, wi = %s, wo = %s\n", phi_i, phi_o, theta_i, theta_o, bRec_left_barbule.wi.toString().c_str(), bRec_left_barbule.wo.toString().c_str());
    //printf("Right barbule: phi_i = %f, phi_o = %f, theta_i = %f, theta_o = %f, wi = %s, wo = %s\n", phi_i, phi_o, theta_i, theta_o, bRec_right_barbule.wi.toString().c_str(), bRec_right_barbule.wo.toString().c_str());

    Spectrum result(0.0f);

    // Debugging

    //result += m_bsdfs[1]->eval(bRec, measure);

    //------------------------------ Fiber BCSDF evaluations using the masking term ---------------------
    
    if (use_masking){

        //printf("Masking computation\n");

        //What is the difference between wi and bRec.its.wi? 
        //Feather_masking::Masking masking = getMasking(barb_f, left_barbule_f, right_barbule_f, bRec.its.wi);
        Feather_masking::Masking masking = getMasking(barb_f, left_barbule_f, right_barbule_f, wi);//Masking
        //Feather_masking::Masking masking = getMasking(barb_f, left_barbule_f, right_barbule_f, wo);//Shadowing

        //printf("Barb range size: %d\n", masking.barb_ranges().size());

        int count = 0;

        for (auto barb_range : masking.barb_ranges()){
            m_bsdfs[1]->hmin = barb_range.hmin(); m_bsdfs[1]->hmax = barb_range.hmax();

            //printf("Barb H = [%f, %f]\n", m_bsdfs[1]->hmin, m_bsdfs[1]->hmax);

            result += m_bsdfs[1]->eval(bRec_barb, measure) * barb_range.rate();

            count++;
        }

        //Barbule H range can be update directly.
        m_bsdfs[2]->hmin = masking.left_range().hmin(); m_bsdfs[2]->hmax = masking.left_range().hmax();
        m_bsdfs[3]->hmin = masking.right_range().hmin(); m_bsdfs[3]->hmax = masking.right_range().hmax();

        result += m_bsdfs[2]->eval(bRec_left_barbule, measure) * masking.left_range().rate();
        result += m_bsdfs[3]->eval(bRec_right_barbule, measure) * masking.right_range().rate();

        // Only a specific barbule test
        //result = m_bsdfs[2]->eval(bRec_left_barbule, measure);
        //result += m_bsdfs[3]->eval(bRec_right_barbule, measure);

    }else{
        result += m_bsdfs[1]->eval(bRec_barb, measure);
        //result += m_bsdfs[2]->eval(bRec_left_barbule, measure);
        //result += m_bsdfs[3]->eval(bRec_right_barbule, measure);
    }

    //printf("BSDF radiance = %s\n", result.toString().c_str());

    /*
    if (bRec.typeMask == EAll){
        printf("bRec = %s\n", bRec.toString().c_str());
        printf("measure = %d\n", measure);
        printf("result = %s\n", result.toString().c_str());
        printf("barb_angle = %f\n", barb_angle);
    }
    */
    
    return result; 
}

// TO DO: update this function
void Feather::addChild(const std::string &name, ConfigurableObject *child)
{

    // Here you can access the textures defined inside the xml configuration file
    if (child->getClass()->derivesFrom(MTS_CLASS(Texture)))
    {

        if (name == "feather_texture")
            m_texture = static_cast<Texture *>(child);

        if (name == "color_texture")
            m_color = static_cast<Texture *>(child);

        if (name == "feather_masking_texture")
            m_masking_texture = static_cast<Texture *>(child);
    }
    else
    {

        if (child->getClass()->derivesFrom(MTS_CLASS(BSDF)))
        {
            BSDF *bsdf = static_cast<BSDF *>(child);
            m_bsdfs.push_back(bsdf);
            bsdf->incRef();
        }
        else
        {
            BSDF::addChild(name, child);
        }
    }
}

// TO DO: update this function
std::string Feather::toString() const
{
    std::ostringstream oss;
    oss << "Feather[" << endl
        << "  id = \"" << getID() << "\"," << endl
        << "  featherTexture = " << indent(m_texture->toString()) << "," << endl
        << "  colorTexture = " << indent(m_color->toString()) << "," << endl
        << "  barbs_rgb = " << barbs_rgb.toString() << "," << endl
        << "  rachis_rgb = " << m_color_rachis.toString() << endl
        << "]";
    return oss.str();
}

MTS_NAMESPACE_END