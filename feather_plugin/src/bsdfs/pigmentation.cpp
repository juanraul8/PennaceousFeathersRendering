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

    Original code: https://sites.cs.ucsb.edu/~lingqi/project_page/fur2/index.html
*/

#include "pigmentation.h"

#include <mitsuba/render/scene.h>

MTS_NAMESPACE_BEGIN

// --------------------[Yan et al. 2017] coordinate system transformation

//Local coordinates to longitudinal-azimuthal parameterization for fur fibers.
// See Figure 2: https://sites.cs.ucsb.edu/~lingqi/publications/paper_fur2.pdf

/*
inline Vector4 getFiberAngles(const BSDFSamplingRecord &bRec){

    Vector3 wi = bRec.its.toWorld(bRec.wi);
    Vector3 wo = bRec.its.toWorld(bRec.wo);

    //Vector3 wi = bRec.wi;
    //Vector3 wo = bRec.wo;
   
    //For scenes using the hair intersection shape, it only works the geoFrame. Otherwise, both geoFrame and shFrame area fine
    const Vector3 &u = bRec.its.geoFrame.s;
    const Vector3 &v = bRec.its.geoFrame.t;
    const Vector3 &w = bRec.its.geoFrame.n;


    //printf ("u (%f, %f, %f)\n", u.x, u.y, u.z);
    //printf ("v (%f, %f, %f)\n", v.x, v.y, v.z);
    //printf ("w (%f, %f, %f)\n", w.x, w.y, w.z);

    Float theta_i = M_PI / 2.0f - acosSafe(dot(u, wi));
    Float theta_r = M_PI / 2.0f - acosSafe(dot(u, wo));

    //Float theta_i = asinSafe(dot(u, wi));
    //Float theta_r = asinSafe(dot(u, wo));

    //NB: Clamp to avoid numerical issues
    theta_i = clamp(theta_i, degToRad(-89.5), degToRad(89.5));
    theta_r = clamp(theta_r, degToRad(-89.5), degToRad(89.5));

    Float phi_i = acosSafe(dot(normalizeSafe(wi - dot(u, wi) * u), v)) * (dot(w, wi) > 0 ? 1 : -1);
    Float phi_r = acosSafe(dot(normalizeSafe(wo - dot(u, wo) * u), v)) * (dot(w, wo) > 0 ? 1 : -1);

    //phi_i = phi_i + M_PI;
    //phi_r = phi_r + M_PI;

    return Vector4(theta_i, theta_r, phi_i, phi_r);
}
*/

// TODO (Adrian Jarabo): Get it out of Pigmentation; into an aux functions library.
/*
void Pigmentation::toLocalFiberCoordinates(const BSDFSamplingRecord &bRec, Float &phi_i, Float &phi_o, Float &theta_i, Float &theta_o) const
{
    // [Yan et al. 2017] transformation
    if (debug_hair)
    {
        Vector4 angles = getFiberAngles(bRec);
        theta_i = angles[0];
        theta_o = angles[1];
        phi_i = angles[2];
        phi_o = angles[3];
        // phi = phi_o - phi_i; // NB: order does not seem to matter
    }
    else
    {
        // Our transformation
        Vector2 param_angle = getAngleParameterization(getCorrectFrame(bRec, bRec.wi, use_feather_coords));
        theta_i = param_angle[0];
        phi_i = param_angle[1];

        param_angle = getAngleParameterization(getCorrectFrame(bRec, bRec.wo, use_feather_coords));
        theta_o = param_angle[0];
        phi_o = param_angle[1];
        // phi = phi_o - phi_i;

        Float m_eps = 1.e-2;
        if (!use_feather_coords){
            //Validate again Yan et al. 2017 transformations. (This makes sense only for the furball scene!)
            Vector4 angles = getFiberAngles(bRec);
            Float theta_i_r = angles[0];
            Float theta_o_r = angles[1];
            Float phi_i_r = angles[2];
            Float phi_o_r = angles[3];
            if ( theta_i_r - theta_i > m_eps || phi_i_r - phi_i > m_eps){
                printf(">> Input (Yan et al.) [%f %f] vs Input (Ours) [%f %f] %s \n",theta_i_r,phi_i_r,theta_i,phi_i,bRec.wi.toString().c_str());
            }
            if ( (theta_o_r - theta_o > m_eps || phi_o_r - phi_o > m_eps) && (bRec.wo[0] > Epsilon || bRec.wo[1] > Epsilon || bRec.wo[2] > Epsilon) ){
                printf(">> Output (Yan et al.) [%f %f] vs Output (Ours) [%f %f] %s \n",theta_o_r,phi_o_r,theta_o,phi_o,bRec.wo.toString().c_str());
            }
        }
    }

    theta_i = clamp(theta_i, degToRad(-89.5), degToRad(89.5));
    theta_o = clamp(theta_o, degToRad(-89.5), degToRad(89.5));

    // regularizePhi(phi);//NB: Original code, but it does not seem to affect the results
}

// The main motivation of creating our own transformation is avoiding the toWorld() transformation from [Yan et al. 2017]

inline Vector getCorrectFrame(const BSDFSamplingRecord &bRec, Vector input_vector, bool use_feather_coords){
    
    if (use_feather_coords){
        return input_vector;//Feather BSDF do not require this transformation
    }
    
    const Frame geoFrame = bRec.its.geoFrame;
    const Frame shFrame = bRec.its.shFrame;
    Vector tmpvector = shFrame.toWorld(input_vector);
    tmpvector = geoFrame.toLocal(tmpvector);//Hair shape works only with the geometric frame
    return tmpvector;
    
}

inline Vector2 getAngleParameterization(Vector m_vector){
    //Float theta_i = M_PI / 2.0f - acosSafe(m_vector.z); 
    Float theta_i = M_PI / 2.0f - acosSafe(m_vector.x); 
    
    
    //Vector yz_plane(0.0f,m_vector.y , m_vector.z);
    //yz_plane = normalize(yz_plane);
    //Float phi_i = acosSafe(yz_plane.z ) * (yz_plane.y > 0 ? 1 : -1); 

    //Float phi_i = atan2(m_vector.y, m_vector.z);//Adrian Jarabo
    Float phi_i = atan2(m_vector.z, m_vector.y);//Dario Lanza
    
    return Vector2(theta_i, phi_i);
}

inline Vector getVectorFromParameterizedAngle(BSDFSamplingRecord bRec, Float theta_o, Float phi_o){
    
    Float ul = sin(theta_o);    
    Float vl = cos(theta_o) * sin(phi_o);
    Float wl = cos(theta_o) * cos(phi_o);

    //return Vector(ul, vl, wl);//Adrian Jarabo (it works well for furball)
    return Vector(ul,vl,wl);//Dario Lanza
}
*/

Pigmentation::Pigmentation(const Properties &props) : BSDF(props) {
    
    // If color texture is given, read and convert it into absorption on the fly
    colorTexture = new ConstantSpectrumTexture(props.getSpectrum("color", Spectrum(1.0f)));

    // Color test
    //Float RGB_sigma[3] = {0.19, 0.24, 0.25};
    Float RGB_sigma[3] = {0.19, 0.24, 0.9};
    fur_color = Spectrum(RGB_sigma);

    // Concretation test
    //Float eu = 8.0;
    //Float pheo = 0.0;
    //fur_color = SigmaAFromConcentration(eu, pheo);

    //-----------------Extracting paramters from the XML scene file
    // TO DO: Transform roughness parameter to radians. In this way, the xml input can be provided in degrees (more intuitive)

    beta_m = props.getFloat("beta_m", 0.055);
    beta_n = props.getFloat("beta_n", 0.055);
    absorptionFactor = props.getFloat("absorptionFactor", 4.0);

    //Masking parameters (constructor)
    hmin = -1.0;hmax = 1.0;

    //Color parameters
    // Color scale: scale by colorScale (if > 0)
    // Color gamma: gamma by colorScale (if > 0)
    colorScale = props.getFloat("colorScale", 2.0);
    colorGamma = props.getFloat("colorGamma", 1.0);
    useColorTexture = props.getBoolean("useColorTexture", false);

    //Cuticle parameters
    alpha = props.getFloat("alpha", 0.0);
    layers = props.getFloat("layers", 0.5);

    //alpha = 0.0f;//Debugging for feathers!
    alpha = -degToRad(alpha);// Tilt should always be negative for physical correctness.

    //Cortex parameters
    use_absorption_cortex = props.getBoolean("use_absorption_cortex", false);
    a_outer = props.getFloat("a_outer", 1.0);
    b_outer = props.getFloat("b_outer", 1.0);
    eta_outer = props.getFloat("eta_outer", 1.55);
    absorption_outer = props.getVector("absorption_outer", Vector3(0.0f, 0.0f, 0.0f));

    //printf("Outer Ellipse a = %f, b = %f\n", a_outer, b_outer);

    //Medulla parameters
    use_diffuse_medulla = props.getBoolean("use_diffuse_medulla", false);
    a_inner = props.getFloat("a_inner", 1.0);
    b_inner = props.getFloat("b_inner", 1.0);
    
    printf("a_inner = %f, b_inner = %f\n", a_inner, b_inner);
    printf("use_diffuse_medulla = %d\n", use_diffuse_medulla);

    Float cx = props.getFloat("center_inner_x", 0.0);
    Float cy = props.getFloat("center_inner_y", 0.0);

    medulla_center = Point2d(cx, cy);
    
    eta_inner = props.getFloat("eta_inner", 1.55);
    //absorption_inner = props.getFloat("absorption_inner", 0.0f);
    use_absorption_medulla = props.getBoolean("use_absorption_medulla", false);
    absorption_inner = props.getVector("absorption_inner", Vector3(0.0f, 0.0f, 0.0f));
    scattering_inner = props.getFloat("scattering_inner", 0.5f);
    medulla_ratio = props.getFloat("medulla_ratio", 0.0);
    medulla_reflectance = props.getSpectrum("medulla_reflectance", Spectrum(1.0f));//Default value is a white medulla

    printf("Inner Ellipse a = %f, b = %f\n", a_inner, b_inner);

    // Thin film parameters
    use_thin_film = props.getBoolean("use_thin_film", false);
    thin_film_IOR = props.getFloat("thin_film_IOR", 1.33);
    thickness = props.getFloat("thickness", thickness);//nanometres

    // Juan Raul Padron Griffe
    if( use_thin_film ){
        SLog(EWarn, "pigmentation: Thin Film Interference is still experimental; recommended to use 'use_thin_film=false' instead.\n");
    }

    //General parameters

    fiber_max_path = props.getInteger("fiber_max_path", 3);
    use_feather_coords = props.getBoolean("use_feather_coords", false);
    use_default_sampling = props.getBoolean("use_default_sampling", false);
    debug_hair = props.getBoolean("debug_hair", false);
    use_h_range = props.getBoolean("use_h_range", false);
    eval_max = props.getFloat("eval_max", 1000.0f);
    sample_max = props.getFloat("sample_max", 2.0f);
    fixed_path = props.getString("fixed_path", "");//Options: "R", "TT", "TRT", ... "" Russian Roulette exploration

    // Juan Raul: Use debug_hair = true for render with hair scenes like straight hair or curly hair
    if( debug_hair ){
        SLog(EWarn, "pigmentation: Using Yan et al.'s original coordinate system; recommended to use 'debug_hair=false' for feather scenes instead.\n");
    }

    // Adrian Jarabo
    if( !use_default_sampling ){
        SLog(EWarn, "pigmentation: Importance sampling is still experimental; recommended to use 'use_default_sampling=true' instead.\n");
    }

    // Juan Raul: If this number is 10 or more, the our current solution crashes with a segmentation fault (core dumped) error
    if( fiber_max_path >= 10){
        SLog(EWarn, "pigmentation: fiber_max_path larger than 9 lobes might cause memory issues; recommended to use fiber_max_path=5 .\n");
    }

    //Lobe selection (mostly debugging)
    
    /*
    std::string active = props.getString("active_lobes", "111");

    printf("Active lobes (binary) = %s\n", active.c_str());

    for( int i = 0; i < active.length(); i++)
    {
        active_lobes[i] = active[i] == '1';
    }

    //printf("Active lobes [%d, %d, %d]\n", active_lobes[0], active_lobes[1], active_lobes[2], active_lobes[3]);
    */

    Spectrum sigma_a;
    color_input_type = props.getInteger("color_input_type", 0);

    // Absorption from the concentration of the pigments
    // Juan Raul [TO DO]: Extend these methods for Absorption
    if (color_input_type == 1){

        printf("Absorption coefficient from eumelanin and pheomelanin concentrations\n");

        /*Float eu = 0.3;
        Float pheo = 0.0;
        Spectrum fur_color = SigmaAFromConcentration(eu, pheo);*/

        eu = props.getFloat("eumelanin", -1.0); //Eumulanin concentration
        pheo = props.getFloat("pheomelanin", -1.0); //Eumulanin concentration
        
        printf("eumelanin = %f, pheomelanin = %f\n", eu, pheo);

        sigma_a = SigmaAFromConcentration(eu, pheo);
        sigma_a.toLinearRGB(absorption_outer.x, absorption_outer.y, absorption_outer.z);
    }

    // Absorption from direct color provided by the user
    if (color_input_type == 2){
        
        printf("Absorption from direct reflectance\n");
        fur_color = props.getSpectrum("fur_color", Spectrum(0.0));

        fur_color = Spectrum(0.0);
        fur_color[0] = 0.18; fur_color[1] = 0.26; fur_color[2] = 0.34;
        beta_n = 0.3f;

        fur_color.toSRGB(fur_color2.x, fur_color2.y, fur_color2.z);
        fur_color = Spectrum(0.0f);
        fur_color[0] = fur_color2.x;fur_color[1] = fur_color2.y; fur_color[2] = fur_color2.z;

        sigma_a = SigmaAFromReflectance(fur_color, beta_n);
        sigma_a.toLinearRGB(absorption_outer.x, absorption_outer.y, absorption_outer.z);
        //fur_color.toLinearRGB(fur_color2.x, fur_color2.y, fur_color2.z);
    }

    printf("absorption_outer = (%f, %f, %f)\n", absorption_outer.x, absorption_outer.y, absorption_outer.z);

    //std::srand(83);//Helpful to debug variable initializations and render consistency

    //solver->setParameters(eta_outer, absorption_outer, alpha, beta_m, beta_n, medulla_ratio, absorption_inner, scattering_inner, layers);

    //fiber_path_tracer = new FiberTracer();

    //--------------------------Unit tests and validations---------------------------------------
    //beta_n = 0.08726646;// Azimuthal roughness of the cuticle (0.138, 0.120, 0.02) 
    //PigmentationSolver fur_solver(eta_outer, Vector3(0.0), alpha, beta_m, beta_n, 0.0, 0.0, 0.0, 0.5);
    //PigmentationSolver fur_solver(eta_outer, Vector3(0.0), alpha, beta_m, beta_n, 0.0, Vector3(0.0), 0.0, 0.5);

    //bool valid;

    //Prediction angles validation
    //valid = fur_solver.ComputeAzimuthalAngles(true, 1.55);//Circular cross section
    //valid = fur_solver.ComputeAzimuthalAngles(1.55);//Elliptical cross section

    //Fresnel values validation
    //valid = fur_solver.ComputeRadiance(true, 1.55, 1.0, 1.0);//Circular cross section
    //valid = fur_solver.ComputeRadiance(false, 1.55, 1.0, 1.0);//Elliptical cross section

    //Azimuthal scattering prediction validation (it might take like a minute for the azimuthal cross section of the BCSDF and a 1 degree resolution for all angles)
    //a = 1.0; b = 1.0;
    //valid = fur_solver.ComputeAzimuthalScattering(a, b);//Elliptical cross section

    //Full scattering prediction validation (it might take around two minutes for the full BCSDF evaluation for a 1 degree resolution for all angles)
    //a = 1.5; b = 1.0;
    //valid = fur_solver.ComputeFullScattering(a, b);//Elliptical cross section

    printf("Initialization Done\n");
}

//Destructor is important to release the memory allocated by the pointer variables
Pigmentation::~Pigmentation(){
    
    // protected destructors!
    //delete this->colorTexture;
    
    //delete this->fiber_path_tracer;
    //fiber_path_tracer = NULL;
}

Pigmentation::Pigmentation(Stream *stream, InstanceManager *manager) : BSDF(stream, manager) {
    colorTexture = static_cast<Texture *>(manager->getInstance(stream));

    useColorTexture = stream->readBool();
    medulla_ratio = stream->readFloat();
    eta_outer = stream->readFloat();
    a_outer = stream->readFloat();
    b_outer = stream->readFloat();
    eta_inner = stream->readFloat();
    a_inner = stream->readFloat();
    b_inner = stream->readFloat();
    eta_outer = stream->readFloat();
    beta_m = stream->readFloat();
    beta_n = stream->readFloat();
    alpha = stream->readFloat();
    layers = stream->readFloat();
    absorptionFactor = stream->readFloat();
    absorption_outer.x = stream->readFloat();
    absorption_outer.y = stream->readFloat();
    absorption_outer.z = stream->readFloat();
    
    absorption_inner.x = stream->readFloat();
    absorption_inner.y = stream->readFloat();
    absorption_inner.z = stream->readFloat();

    scattering_inner = stream->readFloat();

    colorScale = stream->readFloat();
    colorGamma = stream->readFloat();
    fiber_max_path = stream->readInt();

    configure();
}

void Pigmentation::configure() {
    
    /* Verify the input parameter and fix them if necessary */
    colorTexture = ensureEnergyConservation(colorTexture, "color", 1.0f);

    m_components.clear();
    
    if (colorTexture->getMaximum().max() > 0)
        m_components.push_back(EDiffuseReflection | EFrontSide
            | (colorTexture->isConstant() ? 0 : ESpatiallyVarying));
        m_usesRayDifferentials = colorTexture->usesRayDifferentials();

    BSDF::configure();
}

Spectrum Pigmentation::getDiffuseReflectance(const Intersection &its) const {
    return colorTexture->eval(its);
}

inline Vector3 Pigmentation::absorptionFromColorTexture() const{

    printf ("Absorption from color texture\n");

    return Vector3(0.0);

    //Spectrum color = colorTexture->eval(bRec.its, false);//Original Lin Qi code
    
    /*
    Spectrum color = fur_color;

    Float color_r, color_g, color_b;
    color.toLinearRGB(color_r, color_g, color_b);

    color_r = clamp(std::pow(color_r, Float(1.0 / colorGamma)) * colorScale, Float(1.0 / 256.0), Float(0.99));
    color_g = clamp(std::pow(color_g, Float(1.0 / colorGamma)) * colorScale, Float(1.0 / 256.0), Float(0.99));
    color_b = clamp(std::pow(color_b, Float(1.0 / colorGamma)) * colorScale, Float(1.0 / 256.0), Float(0.99));

    Float medulla_occlusion = std::pow(medulla_ratio, 1.5);
    Float absorption_r = (-log(color_r) - absorptionFactor * medulla_occlusion * absorption_inner) / (absorptionFactor - absorptionFactor * medulla_occlusion);
    Float absorption_g = (-log(color_g) - absorptionFactor * medulla_occlusion * absorption_inner) / (absorptionFactor - absorptionFactor * medulla_occlusion);
    Float absorption_b = (-log(color_b) - absorptionFactor * medulla_occlusion * absorption_inner) / (absorptionFactor - absorptionFactor * medulla_occlusion);

    absorption_r = clamp(absorption_r, 0.01, MAX_ABSORPTION_OUTER);
    absorption_g = clamp(absorption_g, 0.01, MAX_ABSORPTION_OUTER);
    absorption_b = clamp(absorption_b, 0.01, MAX_ABSORPTION_OUTER);

    Vector3 absorption = Vector3(absorption_r, absorption_g, absorption_b);

    return absorption; 
    */
}

inline std::vector<Ellipse> Pigmentation::getFiberScene(Vector3 absorption) const {

    //printf("Creating the fiber scene (implicit cylinders)\n");

    /*
    //Some scenarios could be:
    Ellipse inner_ellipse = Ellipse(a_inner, b_inner, eta_outer, eta_outer);//Setting a unit circle with the same IOR as the outer circle
    
    //Non-concentric case
    //inner_ellipse.center.x = 0.1; inner_ellipse.center.y = 0.2;
    */

    std::vector<Ellipse> ellipses;

    // Setting the ouside medium
    Float eta_outside = 1.0;//Air (this could be a parameter too)

    // -------------------Setting the cortex material

    //Float eta_outer = 1.0; //Outer ellipse (testing case, the result should be a black image)
    //a = 1.0, b = 1.0;//Setting a unit circle
    Ellipse outer_ellipse = Ellipse(a_outer, b_outer, eta_outside, eta_outer);//Setting a unit circle
    
    if (use_absorption_cortex){

        outer_ellipse.absorption_mat = true;
        //outer_ellipse.absorption = Vector3(0.0);//No absorption scenario
        //outer_ellipse.absorption = absorption;
        outer_ellipse.absorption = absorption_outer;
    }

    //Diffuse test
    //outer_ellipse.diffuse_mat = true;

    //Thin film interference test
    if (use_thin_film){

        //printf("Using thin film interference material\n");

        outer_ellipse.thin_film_mat = true;
        outer_ellipse.thickness = this->thickness;
        outer_ellipse.thin_film_IOR = this->thin_film_IOR;
    }

    //outer_ellipse.thin_film_mat = true;

    //Float e = outer_ellipse.computeEccentricity();

    ellipses.push_back(outer_ellipse);

    // -------------------Setting the diffuse medulla material
    // a = 0 or b = 0 are degenerate cases
    if (use_diffuse_medulla && a_inner > 0.0f && b_inner > 0.0f){
        
        //printf("Using diffuse medulla for the scene\n");

        Ellipse inner_ellipse = Ellipse(a_inner, b_inner, eta_outer, eta_inner);

        inner_ellipse.center = this->medulla_center;
        inner_ellipse.diffuse_mat = true;

        //medulla_reflectance.toLinearRGB(inner_ellipse.diffuse_reflectance.x, inner_ellipse.diffuse_reflectance.y, inner_ellipse.diffuse_reflectance.z);
        medulla_reflectance.toSRGB(inner_ellipse.diffuse_reflectance.x, inner_ellipse.diffuse_reflectance.y, inner_ellipse.diffuse_reflectance.z);

        //printf("Inner Ellipse2 = (%f, %f, %f)\n", inner_ellipse.diffuse_reflectance.x, inner_ellipse.diffuse_reflectance.y, inner_ellipse.diffuse_reflectance.z);
        //e = inner_ellipse.computeEccentricity();

        ellipses.push_back(inner_ellipse);
    }

    // -------------------Setting the absorbing medulla material
    // a = 0 or b = 0 are degenerate cases
    if (use_absorption_medulla && a_inner > 0.0f && b_inner > 0.0f){
        
        //printf("Using absorbing medulla for the scene\n");

        Ellipse inner_ellipse = Ellipse(a_inner, b_inner, eta_outer, eta_inner);

        inner_ellipse.absorption_mat = true;
        //inner_ellipse.absorption = Vector3(0.0);//No absorption scenario
        inner_ellipse.absorption = absorption_inner;

        ellipses.push_back(inner_ellipse);
    }

    return ellipses;
    
}

inline Float Pigmentation::sampleH(Float h_eps, Float &pdf_h, Sampler *sampler) const{

    Float h = 0.0;
    Float rand_number = 0.0f;

    if (sampler == NULL){
        rand_number = randFloat();
    }else{

        //printf("FiberTracer (rayTrace) using sampler\n");

        rand_number = sampler->next1D();
    }

    //NB: we use eps to avoid the numerical singularities. of our fiber tracer when h = -1 or h = 1. Traditionally, h range is [-1, 1]
    if (use_h_range){
        
        //printf("H range form masking\n");

        h = (rand_number * (hmax-hmin) + hmin); //Far field approximation relying on Monte Carlo Integration h = [hmin, hmax] from Masking 
        pdf_h = 1.0f / (hmax - hmin);
    }else{
        //printf("Default H range\n");
        
        h = (2.0f - 2.0f * h_eps) * rand_number - (1.0f - h_eps); //Far field approximation relying on Monte Carlo Integration h = [-1 + eps, 1 +eps]
        pdf_h = 1.0f / (2.0f - 2.0f * h_eps);
    }

    /*
    if (pdf_h <= 0){
        printf("Negative o 0 pdf: pdf_h = %f, h = [%f. %f]. This should be impossible\n", pdf_h, hmin, hmax);
    }
    */

    //h = 0.0f;
    //pdf_h = 1.0f;

    /*
    if (hmin > -1.0 || hmax < 1.0){
        printf("H = [%f, %f], h = %f\n", hmin, hmax, h);  
    }
    */
    
    //printf("H = [%f, %f], h = %f\n", hmin, hmax, h);

    return h;
}

// Spectrum to RGB safe transformation
// NB: Taking of numerical issues: (1) NAN --> Black (2) Clipping to [0.0, max_value]. For eval is usually 1000, while for sample_max is usually 2.0
inline Spectrum toSpectrum(const Vector3f &value, Float max_value)
{
    Vector3f new_value;

    if (std::isnan(value[0]) || std::isnan(value[1]) || std::isnan(value[2])){

        printf("Pigmentation BCSDF (NAN cases): %f\n", max_value);

        return Spectrum(0.0f);
    }

    if (std::isinf(value[0]) || std::isinf(value[1]) || std::isinf(value[2])){
        
        printf("Pigmentation BCSDF (infinity cases): %f\n", max_value);

        return Spectrum(0.0f);
    }
    
    //Clamp weights [0, max_value]
    for (int i = 0; i < 3; i++){
        new_value[i] = std::min(std::max(Float(0.0f), value[i]), max_value);
    }

    Spectrum retVal(0.0f);
    retVal.fromLinearRGB(new_value[0], new_value[1], new_value[2]);

    return retVal;
}
    
Spectrum Pigmentation::eval(const BSDFSamplingRecord &bRec, EMeasure measure) const{

    //printf("absorption_outer = (%f, %f, %f)\n", absorption_outer.x, absorption_outer.y, absorption_outer.z);

    /*
    if (fiber_path_tracer == NULL){
        return Spectrum(0.0f); 
    }
    */

    // Debugging sampling
    //return Spectrum(0.0);

    //Debugging new H range values

    if (hmin < -1.0 || hmax > 1.0){
        
        printf("Pigmentation (eval): Invalid H range = [%f, %f]\n", hmin, hmax); 

        Spectrum test = Spectrum(0.0f);
        test[2] = 1.0;

        return test; 
    }

    //return Spectrum(0.0f);
    
    /*
    if (hmin > hmax){
        Spectrum test;
        test[0] = 1.0f; test[1] = 0.0f;test[2] = 0.0f;

        return test;
    }else{
        return Spectrum(0.0f);
    }
    */

    /*
    printf("Pigmentation: its.geoFrame.s %s\n", bRec.its.geoFrame.s.toString().c_str());
    printf("Pigmentation: its.geoFrame.t %s\n", bRec.its.geoFrame.t.toString().c_str());
    printf("Pigmentation: its.geoFrame.n %s\n", bRec.its.geoFrame.n.toString().c_str());
    printf("Pigmentation: its.shFrame.s %s\n", bRec.its.shFrame.s.toString().c_str());
    printf("Pigmentation: its.shFrame.t %s\n", bRec.its.shFrame.t.toString().c_str());
    printf("Pigmentation: its.shFrame.n %s\n", bRec.its.shFrame.n.toString().c_str());
    */

    //printf("Calling BSDF eval\n");

    Vector3 absorption = absorption_outer;
    
    if (useColorTexture && medulla_ratio < 1) {

        printf("Using color texture\n");

        absorption = absorptionFromColorTexture();
    }

    //---------------Sampling angles----------------------------------------------
    //Local coordinates to longitudinal-azimuthal parameterization for fur fibers.
    Float theta_i, theta_r, phi_i, phi_r;

    toLocalFiberCoordinates(bRec, phi_i, phi_r, theta_i, theta_r, debug_hair, use_feather_coords);

    //---------------Sampling h-------------------------
    Float pdf_h = 1.0;

    Float h = sampleH(h_eps, pdf_h, bRec.sampler);

    //printf("h = %f\n", h);

    //--------------------Elliptical Fiber BSDF (3D Approach)--------------------------------------
    // Fiber tracing of a  3D cylinder as a numerical implicit solution

    //NB: medulla_ratio corresponds to kappa and it is important for absorption even for the unscattered lobes

    //theta_i = degToRad(0);//No inclination scenario
    //absorption = Vector3(0.0, 0.0, 0.0); // No absorption sceneario

    Query query = Query(phi_i, phi_r, theta_i, theta_r, h, fiber_max_path, beta_n, beta_m);
    query.ellipses = this->getFiberScene(absorption);
    query.sampler = bRec.sampler;

    Vector3 bsdfCosine = solver.computeBCSDF(query);

    //printf("PDF h: %f\n", pdf_h);

    // bsdfCosine /= pdf_h;//Juan Raul: Why should we divide by a pdf inside a eval? 
    
    //printf("Final line of BSDF eval\n");

    return toSpectrum(bsdfCosine, eval_max);
}


Float Pigmentation::pdf(const BSDFSamplingRecord &bRec, EMeasure measure) const {
    
    // Debugging sampling
    // return 0.0f;

    /*
    if (fiber_path_tracer == NULL){
        return 0.0f; 
    }
    */

    // Diffuse reflection pdf
    /*
    if (!(bRec.typeMask & EDiffuseReflection)){

        printf("Diffuse reflection sample? \n");

        return 0.0;
    }
    */

    // Adrian Jarabo
    // Literally, the worst possible sampling. But at least is unbiased!
    if (use_default_sampling)
    {
        return warp::squareToUniformSpherePdf();
    }

    Vector3 absorption = absorption_outer;
    
    if (useColorTexture && medulla_ratio < 1) {

        absorption = absorptionFromColorTexture();
    }

    //---------------Sampling angles----------------------------------------------
    //Local coordinates to longitudinal-azimuthal parameterization for fur fibers.
    Float theta_i, theta_r, phi_i, phi_r;

    toLocalFiberCoordinates(bRec, phi_i, phi_r, theta_i, theta_r, debug_hair, use_feather_coords);

    //---------------Sampling h-------------------------    
    Float pdf_h = 1.0f;;
    Float h = sampleH(h_eps, pdf_h, bRec.sampler);

    //printf("h = %f\n", h);

    //--------------------Elliptical Fiber BSDF (3D Approach)--------------------------------------

    // Fiber tracing of a  3D cylinder as a numerical implicit solution
    Query query = Query(phi_i, phi_r, theta_i, theta_r, h, fiber_max_path, beta_n, beta_m);
    query.ellipses = this->getFiberScene(absorption);
    query.sampler = bRec.sampler;

    Float pdf = solver.pdfBCSDF(query);

    // pdf approach similar to Fur BCSDF
    return std::max(Float(0.0f), pdf);
    //return std::max(Float(0.0f), pdf * pdf_h);
}

Spectrum Pigmentation::sample(BSDFSamplingRecord &bRec, Float &pdf, const Point2 &sample) const {

    /*
    if (hmin > hmax){
        Spectrum test;
        test[0] = 1.0f;test[1] = 0.0f;test[2] = 0.0f;

        return test;
    }else{
        return Spectrum(0.0f);
    }
    */

    /*
    if (fiber_path_tracer == NULL){
        pdf = 0.0f;
        return Spectrum(0.0f); 
    }
    */

    // Diffuse reflection sample
    /*
    if (!(bRec.typeMask & EDiffuseReflection))
    {
        printf("sample: Diffuse reflection sample? \n");

        pdf = 0.0f;
        return Spectrum(0.0f);
    }
    */

    // Adrian Jarabo
    // Literally, the worst possible sampling. But at least is unbiased!
    if (use_default_sampling)
    {
        //printf("Default sampling\n");
        Point2 sample2 = Point2(randFloat(), randFloat());//Why do we use sample2 instead of sample? To avoid correlations? 
        
        bRec.wo = warp::squareToUniformSphere(sample);
        //bRec.wo = warp::squareToUniformSphere(sample2);
        bRec.eta = 1.0f;
        bRec.sampledComponent = 0;
        bRec.sampledType = EDiffuseReflection;
        pdf = warp::squareToUniformSpherePdf();
        
        Spectrum w = this->eval(bRec, ESolidAngle) / pdf;

        return w;
    }

    Vector3 absorption = absorption_outer;
    
    if (useColorTexture && medulla_ratio < 1) {

        printf("sample: Diffuse reflection sample? \n");
        absorption = absorptionFromColorTexture();
    }

    bRec.eta = 1.0f;//Missing this line produce grid artifacts for multiple scattering (more than 10 bounces) and importance sampling strategy
    bRec.sampledComponent = 0;//Yan et al. uses 0 for R, TT and TRT (cortex) and 1 for TTs and TRTs (medulla)
    bRec.sampledType = EDiffuseReflection;

    //---------------Sampling angles----------------------------------------------
    //Local coordinates to longitudinal-azimuthal parameterization for fur fibers.
    Float theta_i, theta_r, phi_i, phi_r;

    toLocalFiberCoordinates(bRec, phi_i, phi_r, theta_i, theta_r, debug_hair, use_feather_coords);

    /*
    printf("\n>>> %s - %s\n", bRec.wi.toString().c_str(), bRec.wo.toString().c_str());
    Vector t = getVectorFromParameterizedAngle(bRec, theta_i, phi_i);
    printf(">> [%f, %f] %s vs %s\n", phi_i, theta_i, bRec.wi.toString().c_str(), t.toString().c_str());
    */

    //---------------Sampling h-------------------------
    Float pdf_h = 1.0f;
    Float h = sampleH(h_eps, pdf_h, bRec.sampler);

    //printf("h = %f\n", h);

    // -------------General sample function -------------

    // Trying to isolate the bug. If it works with the following lines, then we know the issue is the sample function sampleBCSDF 
    bool use_general_sample = false;

    if (use_general_sample){
        Query query_aux = Query(phi_i, phi_r, theta_i, theta_r, h, fiber_max_path, beta_n, beta_m);
        query_aux.ellipses = this->getFiberScene(absorption);
        query_aux.sampler = bRec.sampler;

        pdf = solver.pdfBCSDF(query_aux);

        if (pdf <= 0.0){
            printf("Feather (sample), general sample function: Negative pdf = %f\n", pdf);
        }

        Spectrum weight_aux = eval(bRec, ESolidAngle) / pdf;

        Vector3 weight_res;

        weight_res[0] = weight_aux[0];weight_res[1] = weight_aux[1];weight_res[2] = weight_aux[2];

        return toSpectrum(weight_res, sample_max);
    }

    //--------------------Elliptical Fiber BSDF (3D Approach)--------------------------------------
    // Fiber tracing of a  3D cylinder as a numerical implicit solution

    // Sample an outgoing direction and compute the corresponding weights.

    Float phi_o, theta_o;
    Float pdfS;
    Vector3 weightS;

    Query query = Query(phi_i, phi_r, theta_i, theta_r, h, fiber_max_path, beta_n, beta_m);
    query.ellipses = this->getFiberScene(absorption);
    query.sampler = bRec.sampler;

    int status = solver.sampleBCSDF(query, phi_o, theta_o, pdfS, weightS);

    if (pdfS <= 0){
        printf("pigmentation BCSDF(sample): Negative or 0 pdf S = %f\n", pdfS);

        pdf = 0.0f;
        return Spectrum(0.0f);
    }

    // No valid sampling. 
    // This case should be impossible in our case
    if (status == -1){
        printf("No valid sampling\n");

        pdf = 0.0f;
        return Spectrum(0.0f);
    }

    regularizePhi(phi_o);
    regularizeTheta(theta_o);

    theta_o = clamp(theta_o, degToRad(-89.5), degToRad(89.5));

    // Convert angles to outgoing direction.
    // This should be an inline function to be more readable
    
    //[Yan et al. 2017] transformation
    if (debug_hair){
        const Vector3 &u = bRec.its.geoFrame.s;
        const Vector3 &v = bRec.its.geoFrame.t;
        const Vector3 &w = bRec.its.geoFrame.n;

        Float ul = sin(theta_o);
        Float vl = cos(theta_o) * cos(phi_o);
        Float wl = cos(theta_o) * sin(phi_o);
        bRec.wo = bRec.its.toLocal(u * ul + v * vl + w * wl);
    }else{


        Vector t = getVectorFromParameterizedAngle(bRec,theta_o,phi_o);
        //printf("Base woT %f %f %f \n",t[0],t[1],t[2]);
        
        // Adrian Jarabo: why do we need this? 
        if (!use_feather_coords){
        
            // printf("Base woT %f %f %f \n",t[0],t[1],t[2]);
            t = bRec.its.geoFrame.toWorld(t);
            //printf("GeoFrame woT %f %f %f \n",t[0],t[1],t[2]);
        
            t = bRec.its.shFrame.toLocal(t);
            //printf("Local woT %f %f %f \n",t[0],t[1],t[2]);
            const Vector3 &u = bRec.its.geoFrame.s;
            const Vector3 &v = bRec.its.geoFrame.t;
            const Vector3 &w = bRec.its.geoFrame.n;

            Float ul = sin(theta_o);
            Float vl = cos(theta_o) * cos(phi_o);
            Float wl = cos(theta_o) * sin(phi_o);
            Vector vec_ref = bRec.its.toLocal(u * ul + v * vl + w * wl);
            
            if((t- vec_ref).length() > 1.e-3 )
            {
                printf(">> %s vs %s\n",t.toString().c_str(), vec_ref.toString().c_str());
            } 
        }

        bRec.wo = t;
    }

    if (pdf_h <= 0){
        printf("Pigmentation (sample). Negative o 0 pdf H: %f. This should be impossible\n", pdf_h);
    }


    // Adrian Jarabo: Coordinate system validation
    if(debug_hair)
    {
        
        Vector4 angles = getFiberAngles(bRec);
        theta_i = angles[0];
        theta_r = angles[1];
        phi_i = angles[2];
        phi_r = angles[3];

        const Vector3 &u = bRec.its.geoFrame.s;
        const Vector3 &v = bRec.its.geoFrame.t;
        const Vector3 &w = bRec.its.geoFrame.n;

        Float ul = sin(theta_r);
        Float vl = cos(theta_r) * cos(phi_r);
        Float wl = cos(theta_r) * sin(phi_r);
        Vector wr = bRec.its.toLocal(u * ul + v * vl + w * wl);

        if((bRec.wo-wr).length() > 1.e-3 )
        {
            printf(">> [%f, %f] vs [%f, %f] %s vs %s\n", phi_o, theta_o, phi_r, theta_r, 
                bRec.wo.toString().c_str(), wr.toString().c_str());
        }

    }

    // Calculate the final sampling weight and pdf.
    Vector3 weight = weightS;
    // Vector3 weight = weightS / pdf_h;

    // Compute the pdf
    pdf = pdfS;
    //pdf = pdfS * pdf_h;
    
    /*
    printf("Pigmentation::sample value: %s\n", weightS.toString().c_str());
    printf("Sampled Direction: %f, %f -> %s   ", phi_o, theta_o, bRec.wo.toString().c_str());
    printf("PDF: %f\n\n", pdf);
    */
    
    return toSpectrum(weight, sample_max);
}

Spectrum Pigmentation::sample(BSDFSamplingRecord &bRec, const Point2 &sample) const {
    
    printf("This sample function is never called\n");

    //Diffuse reflection sample
    if (!(bRec.typeMask & EDiffuseReflection))
    {
        return Spectrum(0.0f);
    }

    Float pdf;
    return Pigmentation::sample(bRec, pdf, sample);
    //return this->sample(bRec, pdf, sample) / pdf;
}

void Pigmentation::addChild(const std::string &name, ConfigurableObject *child) {
    if (child->getClass()->derivesFrom(MTS_CLASS(Texture)) && name == "color") {
        colorTexture = static_cast<Texture *>(child);
    } else {
        BSDF::addChild(name, child);
    }
}

void Pigmentation::serialize(Stream *stream, InstanceManager *manager) const {
    BSDF::serialize(stream, manager);

    manager->serialize(stream, colorTexture.get());

    stream->writeBool(useColorTexture);
    stream->writeFloat(medulla_ratio);
    stream->writeFloat(eta_outer);
    stream->writeFloat(a_outer);
    stream->writeFloat(b_outer);
    stream->writeFloat(eta_inner);
    stream->writeFloat(a_inner);
    stream->writeFloat(b_inner);
    stream->writeFloat(beta_m);
    stream->writeFloat(beta_n);
    stream->writeFloat(alpha);
    stream->writeFloat(layers);
    stream->writeFloat(absorptionFactor);
    stream->writeFloat(absorption_outer.x);
    stream->writeFloat(absorption_outer.y);
    stream->writeFloat(absorption_outer.z);
    
    //stream->writeFloat(absorption_inner);  
    stream->writeFloat(scattering_inner);

    stream->writeFloat(absorption_inner.x);
    stream->writeFloat(absorption_inner.y);
    stream->writeFloat(absorption_inner.z);

    stream->writeFloat(colorScale);
    stream->writeFloat(colorGamma);
}

Float Pigmentation::getRoughness(const Intersection &its, int component) const {
    return std::numeric_limits<Float>::infinity();
}

std::string Pigmentation::toString() const {
    std::ostringstream oss;
    oss << "SmoothDiffuse[" << endl
        << "  id = \"" << getID() << "\"," << endl
        << "  reflectance = " << indent(colorTexture->toString()) << endl
        << "]";
    return oss.str();
}

// Color from eumelanin and pheomelanin concentrations --> Page 34
Spectrum Pigmentation::SigmaAFromConcentration(Float ce, Float cp) {
    
    Float sigma_a[3];
    Float eumelaninSigmaA[3] = {0.419f, 0.697f, 1.37f};
    Float pheomelaninSigmaA[3] = {0.187f, 0.4f, 1.05f};
    
    for (int i = 0; i < 3; ++i)
        sigma_a[i] = (ce * eumelaninSigmaA[i] + cp * pheomelaninSigmaA[i]);
    
    Spectrum result(0.0f);

    return Spectrum(sigma_a);
}

// Color from RGB color --> Page 35
Spectrum Pigmentation::SigmaAFromReflectance(const Spectrum &c, Float beta_n) {
    
    Spectrum sigma_a;
    
    for (int i = 0; i < 3; ++i){
        sigma_a[i] = Sqr(std::log(c[i]) /
                    (5.969f - 0.215f * beta_n + 2.532f * Sqr(beta_n) -
                    10.73f * std::pow(beta_n, 3) + 5.574f * std::pow(beta_n, 4) +
                    0.245f * std::pow(beta_n, 5)));
    }
    
    return sigma_a;
}

MTS_NAMESPACE_END