#pragma once
#include <mitsuba/core/warp.h>
#include <mitsuba/core/plugin.h>
#include <typeinfo>
#include <queue>

#include "fiber_math_utils.h"

MTS_NAMESPACE_BEGIN       

struct Query {

    std::vector<Ellipse> ellipses;
    Sampler *sampler = NULL;
    Float beta_m = 0.055f;
    Float beta_n = 0.055f;
    Float phi_i = 0.0f;
    Float phi_o = 0.0f;
    Float theta_i = 0.0f;
    Float theta_o = 0.0f;
    Float h = 0.0f;
    int max_path = 3;

    // Default constructor. It is really important to initialize each variable to avoid grid artifacts!
    Query(Float phi_i, Float phi_o, Float theta_i, Float theta_o, Float h, int max_path, Float beta_n, Float beta_m){
        this->phi_i = phi_i;
        this->phi_o = phi_o;
        this->theta_i = theta_i;
        this->theta_o = theta_o;
        this->h = h;
        this->max_path = max_path;
        this->beta_n = beta_n;
        this->beta_m = beta_m;
    }
}; 

struct Lobe {

    Ray ray = Ray();//Final Ray collected during the path tracing
    Vector3f radiance = Vector3(1.0f, 1.0f, 1.0f);//The radiance carried along the path (we should call this variable throughput)
    int path_length = 0;//The lenght of the travelled path
    string path = "";//Path name represented as string (e.g. "TT", "TRT")
    int N_M = 0;//The number of interaction with the medulla
    Float pdf = 1.0f;// The pdf of the final lobe
    Float weight = 1.0f;//Weight of the pdf

    // Default constructor. It is really important to initialize each variable to avoid grid artifacts!
    Lobe(){
        
        ray = Ray();

        path_length = 0;
        path = "";
        radiance = Vector3f(1.0f);
        N_M = 0;
        pdf = 1.0f;
        weight = 1.0f;
    }
};  

// ---------------------------------Handy math functions

inline void mult3(Vector3f &v0, const Vector3f &v1)
{
    v0[0] *= v1[0];
    v0[1] *= v1[1];
    v0[2] *= v1[2];
}   

// The main goal of this class is to obtain the neccesary variables to compute the gaussian/logistic distribution for the
// fur fiber lobes. For instance, the classical lobes for the a single dielectric cylinder are R, TT, TRT, ..
class FiberTracer {

    public:

    //Float RGB_wavelengths[3] = {650.0, 510.0, 475.0};//Red, Green and Blue channels (Diego Bielsa's Implementation)
    Float RGB_wavelengths[3] = {580.0, 550.0, 450.0};//Red, Green and Blue channels (Belcour's Implementation)

    FiberTracer() {
        //std::srand(83);//Helpful to debug variable initializations and render consistency
    }

    // Compute the reflected vector around the normal.
    // In priciple similar to Mitsuba reflection code: https://github.com/mitsuba-renderer/mitsuba/blob/master/src/libcore/util.cpp 
    Vector3f static our_reflect(Vector3f wi, Vector3f n) {

        Vector3f R = 2.0f * dot(wi, n) * n - wi;
        R = normalizeSafe(R);

        return R;
    }


    // Check if the incident angle and reflected angle are the same with respect to the normal
    bool static validateReflectionLaw(Vector3f n, Vector3f I_dir, Vector3f R_dir){

        Float eps = 1e-3;//Tolerance error

        Float cosThetaI = dot(-I_dir, n);//Reverse vector direction to measure the angle properly
        Float cosThetaR = dot(R_dir, n);

        //printf("Assert Reflection Law %f = %f\n", cosThetaI, cosThetaR);

        Float error = abs(cosThetaI - cosThetaR); 

        if (error > eps){
            printf("Assert Reflection Law %f = %f\n", cosThetaI, cosThetaR);
        }

        return error <= eps;
    }

    //Compute the refraction from the incident vector given the ratio of the index of refraction at the interface.
    Vector3f static our_refract(Vector3f wi, Vector3f n, Float eta_ratio) {
        //Mitsuba refraction code: https://github.com/mitsuba-renderer/mitsuba/blob/master/src/libcore/util.cpp
        //Code: similar to Python reference (feather-rendering/scripts/debugging/raytracer2D_nested.py)

        // No refraction, continue the path (remember this function is called with -wi)
        if (eta_ratio == 1.0){
            return -wi;
        }

        Float cos_theta_i = dot(normalizeSafe(wi), normalizeSafe(n));//Incident angle
        //Float cos_theta_i = dot(wi, n);//Incident angle
        eta_ratio = 1.0/eta_ratio;//Remember eta_ratio is eta2/eta1

        Float sign = 1.0;

        //printf("cosThetaI = %f\n", cos_theta_i);

        //Avoid numerical issues related to H = -1
        Float eps = 1e-9;
        
        if (abs(cos_theta_i) < eps){

            cos_theta_i = 0.0;
            //printf("Numerical related issue!");
        }

        // Inside the medium
        if (cos_theta_i < 0.0){
            sign = -1.0;

            //eta_ratio = 1.0/eta_ratio;//This is actually computed outside of this function
            //printf("Inside the medium\n");
        }

        //Using Snell's law, calculate the squared sine of the angle between the normal and the transmitted ray
        Float c = eta_ratio * eta_ratio * (1.0 - cos_theta_i*cos_theta_i);

        //Total internal reflection case
        if (c > 1.0){

            //printf("Total internal reflection case\n");
            return Vector3f(0.0);
        }

        c = std::sqrt(1.0 - c);

        Vector3f R = (eta_ratio * cos_theta_i - sign * c) * n - eta_ratio * wi;
        
        R = normalizeSafe(R);

        return R;
    }

    // Check if the computed T ray direction respects the Snell's Law (https://en.wikipedia.org/wiki/Snell%27s_law)
    bool static validateSnellLaw(Vector3f n, Vector3f I_dir, Vector3f T_dir, double eta_ratio){
        
        double eps = 1e-4;//Tolerance error (1e-6, 1e-3, 0.5) 
        Float sign = -1.0f;

        double I_angle = acos(dot(I_dir, n));
      
         // Inside the medium
        if (I_angle < 0.0){
            sign = 1.0f;
            printf("Snell Law: Inside the medium\n");
        }

        sign = 1.0f;

        //Snell law formulas
        double T_angle = acos(dot(T_dir, sign * n));

        double term1 = sin (I_angle) / sin(T_angle);
        double term2 = eta_ratio;

        double error = abs(term1 - term2); 

        if (error > eps){
            printf("Assert Snell Law %f = %f\n", term1, term2);
        }

        return error <= eps;
    }

    //NB: fresnelDielectric: https://github.com/mitsuba-renderer/mitsuba/blob/cfeb7766e7a1513492451f35dc65b86409655a7b/src/libcore/util.cpp 
    //NB: Fresnel equations: https://en.wikipedia.org/wiki/Fresnel_equations
    double static fresnelDieletric(Vector3f n, Vector3f I_dir, Vector3f T_dir, double eta1, double eta2){

        if (eta1 == eta2){
            
            //printf("Fresnel: Only transmission case\n");
            return 0.0;
        }
        
        double eta_ratio = eta1 / eta2;

        double cos_theta_i = dot(normalizeSafe(I_dir), normalizeSafe(n));//Incident angle
        double cos_theta_t = dot(normalizeSafe(T_dir), normalizeSafe(n));//Transmission angle

        //printf("cosThetaI = %f\n", cos_theta_i);

        //Numerical issues related to H = -1
        double eps = 1e-9;

        if (abs(cos_theta_i ) < eps){
            cos_theta_i = 0.0;
        }

        //Using Snell's law, calculate the squared sine of the angle between the normal and the transmitted ray
        double c = eta_ratio * eta_ratio * (1.0 - cos_theta_i*cos_theta_i);

        //Total internal reflection case
        if (c > 1.0){

            //printf("Total internal reflection case\n");
            return 1.0;
        }

        c = std::sqrt(1.0 - c);

        cos_theta_i = abs(cos_theta_i);
        cos_theta_t = abs(cos_theta_t);

        //NB: This computation is assuming eta1 = 1.0 (air)
        double Rs = (eta1 * cos_theta_i - eta2 * cos_theta_t) / (eta1 * cos_theta_i + eta2 * cos_theta_t);//Reflectance of s-polarized light
        double Rp = (eta1 * cos_theta_t - eta2 * cos_theta_i) / (eta1 * cos_theta_t + eta2 * cos_theta_i);//Reflectance of p-polarized light

        return 0.5 * (Rs * Rs + Rp * Rp);//Reflectance of unpolarized light ("natural light")
    }

    // Random number between [0, 1]
    inline Float static randFloat() {
        return rand() / (Float) RAND_MAX;
    }

    // Sampling a random diffuse direction based on Cosine Weighted Hemisphere Sampling
    // Experimental code: we should probably check this carefully later
    Vector3f static sampleDiffuseDirection(Sampler *sampler){

        //printf("Sampling diffuse direction \n");

        Point2 sample;
        
        //Using the sampler is probably more effective
        if (sampler == NULL){
            //printf("Sampler hasn't been initialized!\n");
            sample = Point2f(randFloat(), randFloat());
        }else{
            sample = Point2f(sampler->next1D(), sampler->next1D());
        }

        // Different sampling strategies
        Vector3f d = warp::squareToCosineHemisphere(sample);
        //Vector3f d = warp::squareToUniformHemisphere(sample);

        d = normalizeSafe(d);

        return d;
    }

    // Reflection coefficient (s-polarized)
    Float static rs(float eta1, float eta2, float cosI, float cosR) {
        return (eta1 * cosI - eta2 * cosR) / (eta1 * cosI + eta2 * cosR);
    }
     
    // Reflection coefficient (p-polarized)
    Float static rp(float eta1, float eta2, float cosI, float cosR) {
        return (eta2 * cosI - eta1 * cosR) / (eta1 * cosR + eta2 * cosI);
    }
     
    // Transmission coefficient (s-polarized)
    Float static ts(float eta1, float eta2, float cosI, float cosR) {
        return 2.0 * eta1 * cosI / (eta1 * cosI + eta2 * cosR);
    }
     
    // Transmission coefficient (p-polarized)
    Float static tp(float eta1, float eta2, float cosI, float cosR) {
        return 2.0 * eta1 * cosI / (eta1 * cosR + eta2 * cosI);
    }

    // Classical Thin Film reflectance for two layers: external medium, thin film layer and internal medium
    // The idea is replacing the original Fresnel equation with the result of this function (Reflection)
    // Reference: https://gamedev.net/tutorials/_/technical/graphics-programming-and-theory/thin-film-interference-for-computer-graphics-r2962/  
    Float static ThinFilmReflectance(Float cos0, Float lambda, Float thickness, Float eta0, Float eta1, Float eta2){

        // Compute cosines
        // sin1 and sin2 are derived from the Snell Law
        // Basic trigonometry identity: cos² + sin² = 1 
        //Float cos0 = 1;

        Float sin1 = (eta0 / eta1) * (eta0 / eta1) * (1.0 - cos0 * cos0);
        
        if(sin1 > 1.0){
            // Total internal reflection
            return 1.0;
            //return Vector2(1.0, 0.0);
        }

        Float cos1 = sqrt(1.0 - sin1);

        Float sin2 = (eta1 / eta2) * (eta1 / eta2) * (1.0 - cos1 * cos1);
    
        if (sin2 > 1.0){
            // Total internal reflection
            return 1.0;
            //return Vector2(1.0, 0.0);
        }

        Float cos2 = sqrt(1.0 - sin2);

        // Compute phi in order to perform the reflection phase change

        Float d10 = (eta1 >= eta0) ? 0.0 : M_PI;
        Float d12 = (eta1 >= eta2) ? 0.0 : M_PI;
        Float delta = d10 + d12;

        Float phi = (2.0 * M_PI / lambda) * (2.0 * eta1 * thickness * cos1) + delta;

        // Once we obtain the terms of the fresnel equation, we compute each corresponding polarization case
        // (S-polarized light and P-polarized light)

        // Fresnel equations

        Float B_sPolariced = ts(eta0, eta1, cos0, cos1) * ts(eta1, eta2, cos1, cos2);
        Float B_pPolariced = tp(eta0, eta1, cos0, cos1) * tp(eta1, eta2, cos1, cos2);
        Float a_sPolariced = rs(eta1, eta0, cos1, cos0) * rs(eta1, eta2, cos1, cos2);
        Float a_pPolariced = rp(eta1, eta0, cos1, cos0) * rp(eta1, eta2, cos1, cos2);

        Float T_s = B_sPolariced * B_sPolariced / (a_sPolariced * a_sPolariced - 2.0 * a_sPolariced * cos(phi) + 1.0);
        Float T_p = B_pPolariced * B_pPolariced / (a_pPolariced * a_pPolariced - 2.0 * a_pPolariced * cos(phi) + 1.0);
        Float T = (T_s + T_p) / 2.0;

        // Energy conservation
        Float conservation_ratio = (eta2 * cos2) / (eta0 * cos0);

        Float I_T = conservation_ratio * T;
        I_T = min(1.0, max(0.0, 0.0 + I_T));

        Float I_R = min(1.0, max(0.0, 1.0 - I_T));

        return I_R;
        //return Vector2(I_R, I_T);//If we want to work with both!
    }

    // Evaluate transmittance term using the Beer-Lambert Law
    Vector3f static evalTransmittance(Float t, const Vector3f &mu_t)
    {

        /*
        Vector3f transmittance;
        for (int i = 0; i < SPECTRUM_SAMPLES; ++i)
            transmittance[i] = mu_t[i] != 0
                                   ? math::fastexp(mu_t[i] * (-t))
                                   : (Float)1.0f;
        */

        return Vector3f(math::fastexp(-mu_t[0] * t), math::fastexp(-mu_t[1] * t), math::fastexp(-mu_t[2] * t));
        //return exp3(Vector3(-t * mu_t));
    }

     // Implicit path tracing of a nested set of ellipses
    std::vector<Lobe> static rayTrace(Query query, Ray init_ray){

        // Key: collect k final lobes and then perform Russian Roulette. Splitting is helpful to reduce variance, while the Russian Roulette is helpful to avoid
        // the bias of the solution
        // As the medulla is diffuse, the final path would always end in T except for the R lobe. The path can contain only two Ts, one for the cortex and one for the medulla
        // Note: Lets compute all the terms and leave Mitsuba sampling to cancel the respective terms, it might create numerical problems but the implementation is clenar


        //Only for debugging purposes
        //Float phi_i = query.phi_i;
        //Float theta_i = query.theta_i;
        //Float h = query.h;

        vector<Ellipse> ellipses = query.ellipses;
        int max_n_lobes = query.max_path;
        Sampler *sampler = query.sampler;


        //Initial parameters
        Float ray_eps = 1e-6;
        int n_lobes = 0;

        Ray ray = Ray(init_ray);
        Lobe lobe; 
        lobe.ray = Ray(init_ray);

        std::vector<Lobe> lobes;
        lobes.reserve(max_n_lobes + 1);

        int max_iterations = 100;
        int iterations = 0;

        while (n_lobes <= max_n_lobes && iterations < max_iterations ) {
            //----------------- Ray intersection

            ray = lobe.ray;
            iterations++;

            /*
            if (iterations > 20){
                printf("FiberTracer (rayTraceNew): potential loop (number of iterations = %d) = %s\n", iterations, lobe.path.c_str());
            }
            */

            // Current ray:
            ray.o = ray.o + ray_eps * ray.d;//Avoid self intersection
            //ray.d = normalizeSafe(ray.d);//Not neccesary

            Float tmin = std::numeric_limits<Float>::infinity();//default case: no intersection
            int active_idx = -1;//Index of the intersected geometry

            //printf("Ellipses size: %d\n", ellipses.size());

            for (int i = 0; i < ellipses.size();i++){

                Float t = ellipses[i].rayIntersect(ray);

                // Some validations of the ray intersection (rarely, but it happens)
                if (t < 0){
                    printf("FiberTracer (rayTraceNew): t should never be negative = %f\n", t);
                }

                //assert(!(t <= 0));

                //printf("t = %f\n", t);

                if (t < tmin){
                    tmin = t;
                    active_idx = i;
                }
            }

            //printf("Intersection ellipse %d from %d for lobe %s\n", active_idx, ellipses.size(), lobe.path.c_str());
            //printf("Ellipse 1: diffuse = %d vs Ellipse 2 diffuse =%d\n", ellipses[0].diffuse_mat, ellipses[1].diffuse_mat);

            //No intersection with any geometry, just leave the fiber!
            if(active_idx == -1){

                //printf("Fiber Tracer (Degenerate case): %s\n", lobe.path.c_str());

                break;
            }

            // Check that the initial ray is not inside
            /*
            if (lobe.path_length == 0 && cos_theta_i <= 0.0f){
                //printf("Degenerate case: theta_i = %f, h = %f\n", theta_i, h);

                //break;

                printf("Distance travelled: %f - %f\n", tmin, cos_theta_i);

                lobe.radiance = Vector3f(0.0f);
                lobe.pdf = 0.0f;
            }
            */

            //------------------Update the intersection point and query the normal of that point
            Ellipse active_ellipse = ellipses[active_idx];

            Point3f o = ray.o + tmin * ray.d;//Intersection point
            //Point3f o = lobe.ray.o + tmin * lobe.ray.d;//Intersection point
            Vector3f normal = active_ellipse.getNormal(o);//Normal at the intersection point

            //printf("ray = %s, tmin = %f, ray.o (%f, %f, %f)\n", lobe.path.c_str(), tmin, ray.o.x, ray.o.y, ray.o.z);

            //printf("ray = %s, origin = (%f, %f, %f), normal = (%f, %f, %f)\n", lobe.path.c_str(), o.x, o.y, o.z, 
            //       normal.x, normal.y, normal.z);

            //Material properties at the interaction
            double eta1 = active_ellipse.eta_exterior;
            double eta2 = active_ellipse.eta_interior;

            Float eta_ratio = eta2 / eta1; // Outside medium: eta2 > eta1 --> eta_ratio > 1
            //Float eta_ratio = eta1 / eta2;

            Float cos_theta_i = dot(normalizeSafe(-ray.d), normalizeSafe(normal));

            // Degenerate case:  the initial ray is already inside
            // This situation can arise due to numerical singularities. For instance, H = -1 or H = 1.
            // assert(!(lobe.path_length == 0 && cos_theta_i <= 0.0f));
            if (lobe.path_length == 0 && cos_theta_i <= 0.0f){
                //printf("Degenerate case: theta_i = %f, h = %f\n", theta_i, h);
                //printf("Distance travelled: %f - %f\n", tmin, cos_theta_i);

                lobe.radiance = Vector3f(0.0f);
                lobe.pdf = 0.0f;

                lobes[0] = lobe;

                break;
            }

            //printf("path (%s) length = %d, inside = %d\n", lobe.path.c_str(), lobe.path_length, cos_theta_i < 0.0);

            //---------------Absorbing media (Beer-Lambert Law, Attenuation) --> absorbing cortex
            // If distance travelled is 0.0 or the absorption coefficient is (0.0, 0.0, 0.0), we do not have attenuation (good for debugging purposes)
            // Dario Lanza: Absorption coeffient should be computed always except for the initial ray. 
            
            Vector3 absorption = active_ellipse.absorption;//Absorption coefficient

            if (tmin < 0.0){
                printf("rayTraceNew: tmin should not be negative %f\n", tmin);
            }

            // Initial ray should not be consider!
            if (lobe.path_length > 0){

                mult3(lobe.radiance, evalTransmittance(tmin, absorption));
            }
            
            //---------------Diffuse material (Diffuse Structural Color Approximation) --> diffuse medulla
            // The main idea is sampling a ray direction based on Cosine-Weighted Hemisphere Sampling
            // Remember INV_PI is the normalization factor for diffuse materials together with the cosine term 
            if (active_ellipse.diffuse_mat)
            {
                //printf("Diffuse material = (%f, %f, %f)\n", active_ellipse.diffuse_reflectance.x, active_ellipse.diffuse_reflectance.y, active_ellipse.diffuse_reflectance.z);

                //active_ellipse.diffuse_reflectance = Vector3f(1.0, 1.0, 1.0);//default value is white (debugging)
                active_ellipse.diffuse_reflectance = normalize(active_ellipse.diffuse_reflectance);

                //Transform the sampled direction to the local intersection frame (tangent space)
                Vector t(0, 0, 1);// Fiber is aligned along the tangent
                Vector s = cross(normal, t);

                Vector wo = sampleDiffuseDirection(sampler);

                lobe.ray = Ray(o, wo[0] * t + wo[1] * s + wo[2] * normal, Float(0.0));               
                lobe.path += "d";//Inner geometry diffuse reflection (common scenario)
                lobe.path_length++;
                
                Float cos_theta_d = dot(normalizeSafe(lobe.ray.d), normalizeSafe(normal));

                assert(fabs(warp::squareToCosineHemispherePdf(wo) - dot(lobe.ray.d, normal)/M_PI) < 1.e-3);

                /*
                if (cos_theta_d > 0.0){
                    printf("rayTraceNew: diffuse ray should never enter the medulla: %s\n", lobe.path.c_str());
                }
                */

                // Note: the code seems to enter this conditional many times, but the lobe paths does not seem to imply the ray entered the medulla.
                // A ray inside the medulla should be something like TddT 

                if (cos_theta_d < -1.e-3){
                    printf("rayTraceNew: diffuse ray should never enter the medulla (cos_theta_d = %f): %s\n", cos_theta_d, lobe.path.c_str());
                }

                // We always do the explicit computation and let the pdf canceled the correspoding terms
                //Vector3 diffuse_radiance = INV_PI * active_ellipse.diffuse_reflectance * max(cos_theta_d, 0.0f);
                //diffuse_radiance /= warp::squareToCosineHemispherePdf(wo);//Cancels out INV_PI * max (cos_theta_d, 0.0f);

                Vector3 diffuse_radiance = active_ellipse.diffuse_reflectance;

                mult3(lobe.radiance, diffuse_radiance);
                lobe.pdf *= warp::squareToCosineHemispherePdf(wo);//Update the pdf

                lobe.N_M++;

            }else{
                //---------------Dielectric material (aka normal hair BCSDF)

                //Reflection (specular material)
                Ray R_ray;

                R_ray.o = o;
                R_ray.d = our_reflect(-ray.d, normal);
                //R_ray.d = mitsuba::reflect(-lobe.ray.d, normal);
            
                //bool valid = validateReflectionLaw(normal, lobe.ray.d, R_ray.d);
                assert(validateReflectionLaw(normal, ray.d, R_ray.d));

                //Transmission (specular material)
                Ray T_ray;

                T_ray.o = o;
                T_ray.d = our_refract(-ray.d, normal, eta_ratio);
                //T_ray.d = our_refract(-lobe.ray.d, normal, eta_ratio);
                //T_ray.d = mitsuba::refract(ray.d, normal, eta_ratio);

                //bool valid = validateSnellLaw(normal, lobe.ray.d, T_ray.d, eta_ratio);

                //assert(validateSnellLaw(normal, -lobe.ray.d, T_ray.d, eta_ratio));

                //if (!valid){
                //    printf("Snell Law is not valid for path_length = %d, h = %f, theta_i = %f, eta_ratio = %f\n", lobe.path_length, h, theta_i, eta_ratio);
                //}

                //Fresnel dielectric
                Float F_dielectric = fresnelDieletric(normal, ray.d, T_ray.d, eta1, eta2);//3D initial ray

                //Special case (no refraction, the ray continues the straight path)
                if (eta1 == eta2){
                    F_dielectric = 0.0f;
                }

                // Reflection path
                Lobe R_lobe = lobe;

                R_lobe.ray = Ray(R_ray);
                R_lobe.pdf *= F_dielectric; 
                R_lobe.radiance = Vector3(R_lobe.radiance.x * F_dielectric, R_lobe.radiance.y * F_dielectric, R_lobe.radiance.z * F_dielectric);
                R_lobe.path += "R";
                R_lobe.path_length++;

                // Transmission path
                Lobe T_lobe = lobe;

                T_lobe.ray = Ray(T_ray);
                T_lobe.pdf *= (1.0f - F_dielectric); 
                T_lobe.radiance = Vector3(T_lobe.radiance.x * (1.0f - F_dielectric), T_lobe.radiance.y * (1.0f - F_dielectric), T_lobe.radiance.z * (1.0f - F_dielectric));
                T_lobe.path += "T";
                T_lobe.path_length++;

                // First lobe is always the specular reflected lobe
                if (n_lobes == 0){
                    
                    //(Dario Lanza) just swap the eta ratio once. (we assume that there is only 1 interface)
                    // Inside the cortex, we should change the IORs 
                    eta_ratio = 1.0/eta_ratio;

                    lobe = T_lobe;
                    lobes[n_lobes] = R_lobe;
                    n_lobes++;

                    continue;
                }

                // Splitting: collecting all the lobes a long the way of the single path
                if (n_lobes < max_n_lobes){

                    //Dario Lanza, if Total Internal Reflection then we skip an iteration
                    if(T_lobe.ray.d == Vector3f(0,0,0)){
                        lobe = R_lobe;
                        continue;
                    }

                    lobe = R_lobe;
                    lobes[n_lobes] = T_lobe;
                    n_lobes++;

                }else{ // Russian roulette as termination criteria

                    //printf("Running russian_roulette for lobe %s\n", lobe.path.c_str());

                    Float rand_number = 0.0f;

                    if (sampler == NULL){
                        rand_number = randFloat();
                    }else{

                        //printf("FiberTracer (rayTrace) using sampler\n");

                        rand_number = sampler->next1D();
                    }

                    if (rand_number < F_dielectric){

                        lobe = R_lobe;

                    }else{

                        lobes[n_lobes] = T_lobe;
                        break;
                    }

                }
            }
        }

        return lobes;
    }


};  // class FiberTracer

MTS_NAMESPACE_END