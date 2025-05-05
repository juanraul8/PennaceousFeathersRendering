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
    NB: Only near field implementation
*/

#pragma once

#if !defined(__PIGMENTATION_H)
#define __PIGMENTATION_H

#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/texture.h>
#include <mitsuba/hw/basicshader.h>
#include <mitsuba/core/warp.h>
#include "mitsuba/render/scene.h"
#include "../integrators/fiber_path/fiber_path.cpp"
#include <vector>
#include <map>
using namespace std;

//Elliptical BSDF test
#include "fiber_tracer.cpp"

#include <random>

MTS_NAMESPACE_BEGIN

// ================ Pigmentation util functions ================

inline Float Sqr(Float v) { 
    return v * v; 
}

Vector3 normalizeSafe(Vector3 v){
    if (v.length() == 0.0)
        return v;
    else
        return normalize(v);
}

const Float EPS = 1e-6f;

inline Float clamp(Float x, Float a, Float b) {
    if (x < a) return a;
    if (x > b) return b;
    return x;
}

inline int clamp(int x, int a, int b) {
    if (x < a) x = a;
    if (x > b) x = b;
    return x;
}

inline Vector3 clamp(Vector3 x, Float a, Float b) {
    for (int i = 0; i < 3; i++)
        x[i] = clamp(x[i], a, b);
    return x;
}

inline Vector3 max(Vector3 x, Float a) {
    for (int i = 0; i < 3; i++)
        x[i] = std::max(x[i], a);
    return x;
}

inline Vector3 max(Float a, Vector3 x) {
    return max(x, a);
}

inline Vector3 min(Vector3 x, Float a) {
    for (int i = 0; i < 3; i++)
        x[i] = std::min(x[i], a);
    return x;
}

inline Vector3 min(Float a, Vector3 x) {
    return min(x, a);
}

inline Vector3 discardLarger(Vector3 x, Float a) {
    for (int i = 0; i < 3; i++)
        if (x[i] > a){
            x[i] = 0.0f;
        }
    return x;
}

inline Vector3 discardLarger(Float a, Vector3 x) {
    return discardLarger(x, a);
}

//Gaussian function: mean 0, nonnegative std
inline Float G(Float x, Float sigma) {
    if (sigma <= 0.0) {
        //printf("In function G: (sigma = %f) <= 0!\n", sigma);
        return 0.0f;
    }
    return 1.0f / (sqrt(2.0f * M_PI) * sigma) * exp(-0.5f * (x / sigma) * (x / sigma));
}

inline Float GSafe(Float x, Float sigma) {
    if (sigma <= 0.0) {
        return 0.0f;
    }

    Float spreadAngle = sigma / M_PI * 180.0;

    Float gl = 1.0f / (sqrt(2.0f * M_PI) * sigma) * exp(-0.5f * (x / sigma) * (x / sigma));
    Float hf = 1.0f / M_PI;

    if (spreadAngle < 30.0) {
        return gl;
    } else if (spreadAngle > 52.0) {
        return hf;
    } else {
        Float a = (52.0 - spreadAngle) / (52.0 - 30.0);
        return gl * a + hf * (1.0 - a);
    }
}

inline Float laplacian(Float x, Float b) {
    return 1.0 / (2.0 * b) * exp(-fabs(x) / b);
}

inline Float laplacian(Float x, Float mu, Float b) {
    return laplacian(x - mu, b);
}

//Float randFloat() {
//    thread_local std::mt19937 generator(std::random_device{}());
//    std::uniform_real_distribution<Float> distribution(Float(0.0f), Float(1.0f));
//    return distribution(generator);
//}

inline Float randFloat() {
    return rand() / (Float) RAND_MAX;
}

inline Float sampleGaussian(Float mu, Float sigma, Float u1, Float u2) {
    u1 = clamp(u1, EPS, 1.0 - EPS);
    u2 = clamp(u2, EPS, 1.0 - EPS);
    const double two_pi = 2.0 * 3.14159265358979323846;
    double z0, z1;
    z0 = sqrt(-2.0 * log(u1)) * cos(two_pi * u2);
    z1 = sqrt(-2.0 * log(u1)) * sin(two_pi * u2);
    return z0 * sigma + mu;
}

// Float degToRad(Float x) {
//     return x / 180.0 * M_PI;
// }

// Float radToDeg(Float x) {
//     return x / M_PI * 180.0;
// }

inline Float fresnel(Float theta, Float n1, Float n2) {
    if (n1 == n2) return 0.0;

    theta = clamp(theta, -M_PI / 2.0 + EPS, M_PI / 2.0 - EPS);

    Float cosThetaI = cos(theta);
    Float sinThetaI2 = 1.0 - cosThetaI * cosThetaI;
    Float sinThetaT2 = n1 * n1 * sinThetaI2 / (n2 * n2);
    if (sinThetaT2 > 1.0)
        return 1.0;
     
    Float rSFirstTerm = n1 * cosThetaI;
    Float rSSecondTerm = n2 * sqrt(1.0 - sinThetaT2);
    Float rS = (rSFirstTerm - rSSecondTerm) / (rSFirstTerm + rSSecondTerm);
    rS *= rS;
 
    Float rPFirstTerm = n1 * sqrt(1.0 - sinThetaT2);
    Float rPSecondTerm = n2 * cosThetaI;
    Float rP = (rPFirstTerm - rPSecondTerm) / (rPFirstTerm + rPSecondTerm);
    rP *= rP;
 
    return (rS + rP) / 2.0;
}

inline Float fresnelCuticle(Float theta, Float n1, Float n2, Float layers) {
    if (n1 == n2) return 0.0;

    //printf("layers = %f\n", layers);

    theta = clamp(theta, -M_PI / 2.0 + EPS, M_PI / 2.0 - EPS);

    Float cosThetaI = cos(theta);
    Float sinThetaI2 = 1.0 - cosThetaI * cosThetaI;
    Float sinThetaT2 = n1 * n1 * sinThetaI2 / (n2 * n2);
    if (sinThetaT2 > 1.0) //Total internal reflection case
        return 1.0;

    Float rSFirstTerm = n1 * cosThetaI;
    Float rSSecondTerm = n2 * sqrt(1.0 - sinThetaT2);
    Float rS = (rSFirstTerm - rSSecondTerm) / (rSFirstTerm + rSSecondTerm);
    rS *= rS;

    Float rPFirstTerm = n1 * sqrt(1.0 - sinThetaT2);
    Float rPSecondTerm = n2 * cosThetaI;
    Float rP = (rPFirstTerm - rPSecondTerm) / (rPFirstTerm + rPSecondTerm);
    rP *= rP;

    //Different from fresnel function!
    Float phi_1 = 0.5 * (rS + (1.0 - rS) * (1.0 - rS) * rS / (1.0 - rS * rS)) +
                  0.5 * (rP + (1.0 - rP) * (1.0 - rP) * rP / (1.0 - rP * rP));

    Float phi_m = layers * phi_1 / (1.0 + (layers - 1.0) * phi_1);

    return phi_m;
}

// For non-separable R lobe. All in radians.
// A Fiber Scattering Model With Non-Separable Lobes: http://www.eugenedeon.com/wp-content/uploads/2015/03/nonsephair2014slides.pdf (Slide 35)
// If alpha is 0, then the longitudinal outgoing angle is -theta_i (the expected output for a smooth dieletric cylinder)
inline Float theta_cone_R(Float theta_i, Float alpha, Float phi) {
    Float u = cos(phi / 2.0) * cos(alpha) * cos(theta_i) + sin(alpha) * sin(theta_i);
    Float t = sin(theta_i) - 2.0 * sin(alpha) * u;
    return -asinSafe(t);
}

// For non-separable R lobe. All in radians.
// A Fiber Scattering Model With Non-Separable Lobes: http://www.eugenedeon.com/wp-content/uploads/2015/03/nonsephair2014slides.pdf (Slide 36)
inline Float beta_cone_R(Float beta, Float phi) {
    return beta * sqrt(2.0) * cos(phi / 2.0);
}

inline void regularizeTheta(Float &theta) {
    if (theta > M_PI / 2.0f)
        theta = M_PI - theta;
    else if (theta < -M_PI / 2.0f)
        theta = -M_PI - theta;

    theta = clamp(theta, degToRad(-89.5), degToRad(89.5));
}

inline Float congjugateTheta(Float theta) {
    if (theta > 0.0f)
        theta = M_PI - theta;
    else
        theta = -M_PI - theta;
    return theta;
}

inline void regularizePhi(Float &phi) {
    while (phi > M_PI) phi -= 2.0f * M_PI;
    while (phi < -M_PI) phi += 2.0f * M_PI;
}

inline Vector3 pow3(const Vector3 &v, Float p) {
    return Vector3(pow(v[0], p), pow(v[1], p), pow(v[2], p));
}

inline Vector3 mul3(const Vector3 &v1, const Vector3 &v2) {
    return Vector3(v1[0] * v2[0], v1[1] * v2[1], v1[2] * v2[2]);
}

bool cmdOptionExists(char** begin, char** end, const std::string& option) {
    return std::find(begin, end, option) != end;
}

char* getCmdOption(char **begin, char **end, const std::string & option) {
    char **itr = std::find(begin, end, option);
    if (itr != end && ++itr != end) {
        return *itr;
    }
    return 0;
}

template<typename T>
T getCmdValue(int argc, char *argv[], const std::string & option, T defaultValue) {
    char *s = getCmdOption(argv, argv + argc, option);
    return s? T(atof(s)) : defaultValue;
}

// ================ New functions ===================

inline Float Logistic(Float x, Float s) {
    x = std::abs(x);
    return std::exp(-x / s) / (s * Sqr(1 + std::exp(-x / s)));
}

inline Float LogisticCDF(Float x, Float s) {
    return 1 / (1 + std::exp(-x / s));
}

inline Float TrimmedLogistic(Float x, Float s, Float a, Float b) {
    return Logistic(x, s) / (LogisticCDF(b, s) - LogisticCDF(a, s));
}

// ================ Pigmentation Solver class ================

// For dual scattering
const Float MAX_ABSORPTION_OUTER = 2.4f;
//const int NUM_SIGMA_SAMPLES = 128;

Vector3 operator*(const Vector3& lhs, const Vector3& rhs) {
    return Vector3(lhs[0] * rhs[0], lhs[1] * rhs[1], lhs[2] * rhs[2]);
}

void operator*=(Vector3 &lhs, const Vector3& rhs) {
    lhs[0] *= rhs[0];
    lhs[1] *= rhs[1];
    lhs[2] *= rhs[2];
}

Vector3 operator/(const Vector3& lhs, const Vector3& rhs) {
    return Vector3(lhs[0] / rhs[0], lhs[1] / rhs[1], lhs[2] / rhs[2]);
}

void operator/=(Vector3 &lhs, const Vector3& rhs) {
    lhs[0] /= rhs[0];
    lhs[1] /= rhs[1];
    lhs[2] /= rhs[2];
}

Float maxComponent(const Vector3& v) {
    return std::max(std::max(v[0], v[1]), v[2]);
}

class PigmentationSolver {

public:

    PigmentationSolver() {

        //std::srand(83);//Helpful to debug variable initializations and render consistency
    }


    // Paper section 4.2 (Unscattered Lobes): Distribution term (Dp) of the Azimuthal scattering function (Np)
    // The term is approximated as a Gaussian (Eqn. 7).
    // NB: Why is it computed in -1, 0 and 1?
    Float Du(Float phi, Float sigma) const {
        Float retVal = 0.0;
        
        for (int i = -1; i <= 1; i++)
            retVal += G(phi + i * 2.0 * M_PI, sigma);

        //retVal += G(phi, sigma);

        /*int count = 0;
        for (Float t = -1.0; t <= 1.0; t+=0.01){
            retVal += G(phi + t * 2.0 * M_PI, sigma);

            count++;
        }

        //printf("count = %d\n", count);
        retVal /= count;*/


        return retVal;
    }

    // ------------------------------------Unit Tests and Validations---------------------------------------------------------------------

    inline static Float computeRoughness(Float beta_n, int p){

        if (p <= 0){

            printf("Wrong p to compute the update of the roughness\n");

        }

        // Yan et al. 2017 roughness (specular lobes looks softer)
        // Warning: produce grid artifacts with multiple cores, especially if you include a medulla
        //Float beta_n_p_sqr = beta_n * beta_n * p;//Eqn 16
        //Float beta_n_p = std::sqrt(beta_n_p_sqr);

        // PBRT roughness (specular lobes look too strong )
        //Float beta_n_p = 0.626657069f * (0.265f * beta_n + 1.194f * beta_n * beta_n + 5.372f * std::pow(beta_n, 22.0));

        //Khungurn et al. 2017 Fig 5 uses same beta_n for all lobes.
        Float beta_n_p = beta_n;

        return beta_n_p;
    }

    inline static Ray getInitialRay(const std::vector<Ellipse> &ellipses, Float phi_i, Float theta_i, Float h)
    {
        Ray init_ray;

        // The outer ellipse is used to define the fiber azimuthal radius parametrization (h)
        Float outer_a = ellipses[0].a;
        Float outer_b = ellipses[0].b;

        // Adrian Jarabo: Just assume that the ray sampling is centered on the outermost ellipse; that is the obvious
        // choice for an outer ellipse bounding all others.
        Point2d center = ellipses[0].center;

        Float offset = (outer_a + outer_b) / cos(theta_i); // Offset to make sure the ray is in fact outside of the elliptical cylinder

        // Spherical coordinates (Linqi/Marschner scenario, elevation interpretation)
        // I think the negative d is wrong (see the python path tracing viewer)
        init_ray.d = -Vector3f(cos(phi_i) * cos(theta_i), sin(phi_i) * cos(theta_i), sin(theta_i));
        //init_ray.d = Vector3f(cos(phi_i) * cos(theta_i), sin(phi_i) * cos(theta_i), sin(theta_i));

        // Setting the origin with respect to the h (fiber radius)
        init_ray.o = Point3f(-h * outer_a * sin(phi_i) * cos(theta_i) + center.x, h * outer_b * cos(phi_i) * cos(theta_i) + center.y, sin(theta_i));

        // Setting the ray outside of the cross section.
        init_ray.o = init_ray.o - offset * init_ray.d;

        // init_ray.o = Point3f(-h * outer_a * sin(phi_i) + center.x, h * outer_b * cos(phi_i), 0.);
        // init_ray.o = init_ray.o - (outer_a + outer_b) * init_ray.d;

        return init_ray;
    }

    //General BCSDF computation, where the cylinders can be dieletric, absorbing and diffuse.
    static Vector3 computeBCSDF(Query query){

        //printf("Computing new eval: Fiber Tracer max path depth %d\n", max_path);

        vector<Ellipse> ellipses = query.ellipses;
        Float phi_i = query.phi_i;
        Float phi_o = query.phi_o;
        Float theta_i = query.theta_i;
        Float theta_o = query.theta_o;
        Float h = query.h;
        int max_path = query.max_path;
        Float beta_m = query.beta_m;
        Float beta_n = query.beta_n;

        //------------------------ Sanity check and initialization ------------------------

        if (ellipses.empty()){
             printf("computeBCSDF: Ellipses has not been initialized!\n");
             return Vector3(0.0, 0.0, 0.0);
        }

        // Avoid singularity across grazing angles
        if (fabs(theta_i) == M_PI / 2)
        {
            printf("computeBCSDF: singular case, theta_i grazing angles!\n");
            return Vector3(0.0f);
        }

        // printf("Theta_i: %f\n", theta_i);

        //------------------------ Path Tracing (Splitting with Russian roulette termination)
        // Note: an static array (classical or std::array) defined with a macro might be a more effient solution. It would require to compile every time the
        // code for difrent number of lobes

        Ray init_ray = getInitialRay(ellipses, phi_i, theta_i, h);
          
        std::vector<Lobe> lobes;
        int max_n_lobes = max_path;

        FiberTracer fiber_tracer;

        lobes = fiber_tracer.rayTrace(query, init_ray);

        //printf("computeBCSDFNew: Fiber tracer executed correctly\n");

        //------------------------ Evaluate the BCSDF (Sum of the lobes)

        // TO DO
        // Add code to handle the specific fixed path. Helpful for debugging specific lobes like R, TT and TRT
        Vector3 bsdf = Vector3(0.0f, 0.0f, 0.0f);

        for (int k = 0; k < max_n_lobes + 1;k++){

            Lobe lobe = lobes[k];

            //printf("computeBCSDFNew: Evaluating lobe %s\n", lobe.path.c_str());

            // No intersection scenario. Maybe it happens due to numerical issues
            if (lobe.path_length == 0){
                
                //printf("computeBCSDF (DEGENERATE CASE): no intersection scenario for fiber tracer : %s (k = %d)\n", lobe.path.c_str(), k);
                //return Vector3(0.0, 1.0, 1.0);
                continue;
            }

            // Dario Lanza: Angles from outgoing ray of the lobe
            Vector2 out_ray_angles = FiberTracer2BCSDF(lobe.ray.d);
            Float Theta = out_ray_angles[0];
            Float Phi = out_ray_angles[1];

            // Sample lobes: Gaussian sampling
            Float beta_n_p = computeRoughness(beta_n, lobe.path_length); //Eqn 16 (azimuthal roughness per lobe)

            Float query_phi = Phi - phi_o;//Marschner/LinQi
            regularizePhi(query_phi);// Remap azimuthal query to [-pi,pi]

            Float N_D = 0.0f;

            //I guess we query also query_phi-2pi and query_phi+2pi to consider the case where the query phi is outside of the valid range [-pi, pi]
            for (int i = -1; i <= 1; i++){
                N_D += G(query_phi + i * 2.0 * M_PI, beta_n_p);
            }

            Float query_theta = Theta - theta_o;
            regularizeTheta(query_theta);//Not really neccesary

            Float M_D = G(query_theta, beta_m);//Our Longitudinal Gaussian

            bsdf += M_D * N_D * lobe.radiance;//Full computation
            //bsdf += lobe.radiance;//Ignoring gaussian detectors
        }

        // Note: Division by cosine as the eval of a BCSDF f() / cosine^2. Since Mitsuba assumes BSDF return f()*cosine(), the division cancels one of the 
        // cosine

        bsdf /= cos(theta_o);

        return bsdf;
    }

    // --------------------------------------------Sampling functions--------------------------------------------------------------------------

    static Float pdfBCSDF(Query query){
        
        //printf("Computing new pdf\n");

        vector<Ellipse> ellipses = query.ellipses;
        Float phi_i = query.phi_i;
        Float phi_o = query.phi_o;
        Float theta_i = query.theta_i;
        Float theta_o = query.theta_o;
        Float h = query.h;
        int max_path = query.max_path;
        Float beta_m = query.beta_m;
        Float beta_n = query.beta_n;

        //------------------------ Sanity check and initialization ------------------------

        if (ellipses.empty()){
             printf("pdfBCSDF: Ellipses has not been initialized!\n");
             return -1;
        }

        // Avoid singularity across grazing angles
        if (fabs(theta_i) == M_PI / 2)
        {
            printf("pdfBCSDF: singular case, theta_i grazing angles!\n");
            return -1;
        }

        // printf("Theta_i: %f\n", theta_i);

        //------------------------ Path Tracing (Splitting with Russian roulette termination)

        Ray init_ray = getInitialRay(ellipses, phi_i, theta_i, h);
          
        std::vector<Lobe> lobes;
        int max_n_lobes = max_path;
        FiberTracer fiber_tracer;

        lobes = fiber_tracer.rayTrace(query, init_ray);

        //printf("pdfBCSDFNew: Fiber tracer executed correctly\n");

        //------------------------ 
        // We coumpute the CDF of each lobe
        std::vector<Float> cdf;
        cdf.reserve(max_n_lobes + 1);

        //CDF computation
        Float sum_pdf = 0.0f;
        for (int i = 0; i < max_n_lobes + 1; i++){
            sum_pdf += lobes[i].radiance.length();
        }

        // No lobe is selected
        if (sum_pdf == 0.0) {
            
            printf("pdfBCSDF CDF is 0 case: h = %f, theta_i = %f, phi_i = %f\n", h, theta_i, phi_i);
            return -1;
        }
            
        //Normalize CDF
        for (int i = 0; i < max_n_lobes+1; i++){

            cdf[i] =  lobes[i].radiance.length() / sum_pdf;
        }

        //------------------------ Compute the total probability of the lobes (sum)
        Float pdf = 0.0f;

        for (int k = 0; k < max_n_lobes + 1;k++){

            Lobe lobe = lobes[k];

            // No intersection scenario. Maybe it happens due to numerical issues
            if (lobe.path_length == 0){
                
                printf("pdfBCSDF: DEGENERATE CASE: no intersection scenario for lobe %s\n", lobe.path.c_str());
                return -1;
            }

            // Dario Lanza: Angles from outgoing ray of the lobe
            Vector2 out_ray_angles = FiberTracer2BCSDF(lobe.ray.d);
            Float Theta = out_ray_angles[0];
            Float Phi = out_ray_angles[1];

            // Sample lobes: Gaussian sampling
            Float beta_n_p = computeRoughness(beta_n, lobe.path_length); //Eqn 16 (azimuthal roughness per lobe)

            // Implementation similar to Fur BCSDF
            
            //Longitudinal profile
            //Float query_theta = theta_o - (-theta_i);
            Float query_theta = Theta - theta_o;
            regularizeTheta(query_theta);

            Float pdfM_D = G(query_theta, beta_m);
            //pdfM = 1.0;//Debugging
            
            //Float query_phi = Phi - phi;
            Float query_phi = Phi - phi_o;//Elliptical cross section should be like this
            regularizePhi(query_phi);

            Float pdfN_D = G(query_phi, beta_n_p);
            //pdfN = 1.0;

            //-------------------------Uniform sampling
            // NB: if we do not know exactly how to compute it (similar to scattered lobes from Fur BCSDF)
            bool use_uniform_sampling = false;

            if (use_uniform_sampling){
                pdfM_D = 1.0 / M_PI;
                pdfN_D = 1.0 / (2.0 * M_PI);
            }

            pdf += pdfM_D * pdfN_D * cdf [k];
            //pdf += cdf [k];//Ignoring gaussian detectors

        }

        if (pdf <= 0.0){
            //printf("pdfBCSDF: Negative pdf. pdf = %f, cos(theta_o)=%f\n", pdf, cos(theta_o));
            return -1;
        }

        //Solid angle instead of polar coordinates would be without dividing by cos(theta_o)
        pdf /= cos(theta_o);//cos(theta_o) might be negative is the angle is not between [-pi/2, pi/2] 

       return pdf;
    }

   static int sampleBCSDF(Query query, Float &phi_o, Float &theta_o, Float &pdfS, Vector3 &weightS){

        //printf("Computing new sample\n");

        vector<Ellipse> ellipses = query.ellipses;
        Float phi_i = query.phi_i;
        Float theta_i = query.theta_i;
        Float h = query.h;
        int max_path = query.max_path;
        Sampler *sampler = query.sampler;
        Float beta_m = query.beta_m;
        Float beta_n = query.beta_n;

        //------------------------ Sanity check and initialization ------------------------

        if (ellipses.empty()){
             printf("sampleBCSDF: Ellipses has not been initialized!\n");
             return -1;
        }

        // Avoid singularity across grazing angles
        if (fabs(theta_i) == M_PI / 2)
        {
            printf("sampleBCSDF: singular case, theta_i grazing angles!\n");
            return -1;
        }

        // printf("Theta_i: %f\n", theta_i);

        //------------------------ Path Tracing (Splitting with Russian roulette termination)

        Ray init_ray = getInitialRay(ellipses, phi_i, theta_i, h);
          
        std::vector<Lobe> lobes;
        int max_n_lobes = max_path;
        FiberTracer fiber_tracer;


        lobes = fiber_tracer.rayTrace(query, init_ray);

        //printf("sampleBCSDFNew: Fiber tracer executed correctly\n");

        //------------------------ 
        // We pick one lobe from the collected lobes using a discrete variable where the probability is proportional to their throughput
        // The energy that lobe p carries depends on fresnel interactions, attenuation by the absorbing cortex and the diffuse reflectance of the medulla

        std::vector<Float> cdf;
        cdf.reserve(max_n_lobes + 1);

        //CDF computation
        Float sum_pdf = 0;
        for (int i = 0; i < max_n_lobes + 1; i++){
            sum_pdf += lobes[i].radiance.length();
        }

        // No lobe is selected
        if (sum_pdf == 0.0) {
            
            printf ("sampleBCSDF: Integral CDF is 0 --> No lobe selected!\n");
            //printf("CDF is 0 case: h = %f, theta_i = %f, phi_i = %f\n", h, theta_i, phi_i);

            return -1;
        }
            
        //Normalize CDF
        for (int i = 0; i < max_n_lobes+1; i++){

            cdf[i] =  lobes[i].radiance.length() / sum_pdf;
        }

        int lobe_idx = max_n_lobes;
        Float rand_num = 0.0f;

        if( sampler == NULL )
        {
            //printf("sampleBCSDFNew: using randFloat()\n");
            rand_num = randFloat();
        }
        else
        {
            //printf("sampleBCSDFNew: using sampler\n");
            rand_num = sampler->next1D();
        }
            
        for (int i = 0; i < max_n_lobes+1; i++) {
            rand_num -= cdf[i];
            
            if (rand_num <= 0.0) {
                lobe_idx = i;
                break;
            }
        }

        Lobe lobe = lobes[lobe_idx];
        Float pdf_lobe = cdf[lobe_idx];

        if (pdf_lobe > 1.0){
            printf("sampleBCSDF: pdf lobe > 1 = %f, %s\n", pdf_lobe, lobe.path.c_str());
        }

        if (pdf_lobe <= 0){
            printf("sampleBCSD: Negative or 0 pdf lobe = %f\n", pdf_lobe);

            return -1;
        }

        //------------------------ Evaluate the corresponding lobe

        // No intersection scenario. Maybe it happens due to numerical issues
        if (lobe.path_length == 0){
            
            printf("sampleBCSDF (DEGENERATE CASE): no intersection scenario for sampled lobe: %s\n", lobe.path.c_str());
            return -1;
        }

        // Dario Lanza: Angles from outgoing ray of the lobe
        Vector2 out_ray_angles = FiberTracer2BCSDF(lobe.ray.d);
        Float Theta = out_ray_angles[0];
        Float Phi = out_ray_angles[1]; 

        // Sample lobes: Gaussian sampling
        Float beta_n_p = computeRoughness(beta_n, lobe.path_length); //Eqn 16 (azimuthal roughness per lobe)

        /*
        printf("Lobe radiance = (%f, %f, %f), Lobe radiance (length)= %f vs lobe pdf = %f\n", lobe.radiance.x, lobe.radiance.y, lobe.radiance.z,
            lobe.radiance.length(), lobe.pdf);
        */

        // Implementation similar to Fur BCSDF
        // Gaussian sampling for the query outgoing ray produces energy conservation issues once we introduce a medulla inside the cortex
        // we can apply the uniform sampling strategy only if the lobe interact with the medulla

        /*
        if (lobe.N_M == 0){

            //Gaussian sampling around the predicted outgoing ray
            
            //printf("sampleBCSDF: Using Gaussian sampling for queried outgoing angles\n");
            
            if( sampler == NULL )
            {
                //printf("sampleBCSDFNew: using randFloat()\n");

                theta_o = sampleGaussian(Theta, beta_m, randFloat(), randFloat());// Longitudinal profile
                phi_o = sampleGaussian(Phi, beta_n_p, randFloat(), randFloat());// Azimuthal profile
            }
            else
            {
                //printf("sampleBCSDFNew: using sampler\n");

                // Gaussian sampling
                theta_o = sampleGaussian(Theta, beta_m,  sampler->next1D(), sampler->next1D());// Longitudinal profile   
                phi_o = sampleGaussian(Phi, beta_n_p, sampler->next1D(), sampler->next1D());// Azimuthal profile
            }   
        }else{

            //Uniform sampling around the whole sphere

            //printf("sampleBCSDF: Using uniform sampling for queried outgoing angles\n");

            if( sampler == NULL )
            {
                //printf("sampleBCSDFNew: using randFloat()\n");

                theta_o = asinSafe(randFloat() * 2.0 - 1.0);
                Float phi_o = (randFloat() * 2.0 - 1.0) * M_PI;//[-PI, PI]
            }
            else
            {
                //printf("sampleBCSDFNew: using sampler\n");

                // Uniform sampling
                theta_o = asinSafe(sampler->next1D() * 2.0 - 1.0);
                Float phi_o = (sampler->next1D() * 2.0 - 1.0) * M_PI;//[-PI, PI]
            }

        }

        */

        //Direct sampling

        //theta_o = Theta;
        //phi_o = Phi;

        // Gaussian sampling around predicted outgoing ray
        if( sampler == NULL )
        {
            //printf("sampleBCSDFNew: using randFloat()\n");

            theta_o = sampleGaussian(Theta, beta_m, randFloat(), randFloat());// Longitudinal profile
            phi_o = sampleGaussian(Phi, beta_n_p, randFloat(), randFloat());// Azimuthal profile
        }
        else
        {
            //printf("sampleBCSDFNew: using sampler\n");

            // Gaussian sampling
            theta_o = sampleGaussian(Theta, beta_m,  sampler->next1D(), sampler->next1D());// Longitudinal profile   
            phi_o = sampleGaussian(Phi, beta_n_p, sampler->next1D(), sampler->next1D());// Azimuthal profile
        } 

        //Uniform sampling around the whole sphere

        //printf("sampleBCSDF: Using uniform sampling for queried outgoing angles\n");

        /*
        if( sampler == NULL )
        {
            //printf("sampleBCSDFNew: using randFloat()\n");

            theta_o = asinSafe(randFloat() * 2.0 - 1.0);
            Float phi_o = (randFloat() * 2.0 - 1.0) * M_PI;//[-PI, PI]
        }
        else
        {
            //printf("sampleBCSDFNew: using sampler\n");

            // Uniform sampling
            theta_o = asinSafe(sampler->next1D() * 2.0 - 1.0);
            Float phi_o = (sampler->next1D() * 2.0 - 1.0) * M_PI;//[-PI, PI]
        }
        */

        // It is important to remap theta_o to [-pi/2, pi/2] to avoid the cos(o) < 0 which would implies pdf_s to be negative too.
        regularizePhi(phi_o);
        regularizeTheta(theta_o);
        

        //Longitudinal profile
        //Float query_theta = theta_o - (-theta_i);
        Float query_theta = Theta - theta_o;
        regularizeTheta(query_theta);

        Float pdfM_D = G(query_theta, beta_m);
        //pdfM_D = 1.0;//Debugging
        Float weightM = 1.0;
        
        //Float query_phi = Phi - phi;
        Float query_phi = Phi - phi_o;//Elliptical cross section should be like this
        regularizePhi(query_phi);

        Float pdfN_D = G(query_phi, beta_n_p);
        //pdfN_D = 1.0;
        
        //Vector3 weightN = Vector3(1.0);
        Vector3 weightN = lobe.radiance;
        //Vector3 weightN = max(discardLarger(lobe.radiance, 5.0), 0.0);//Not really helpful

        //-------------------------Uniform sampling
        // NB: if we do not know exactly how to compute it (similar to scattered lobes from Fur BCSDF)
        
        bool use_uniform_sampling = false;

        if (use_uniform_sampling){
            
            //Longitudinal profile
            theta_o = asinSafe(sampler->next1D() * 2.0 - 1.0);
            theta_o = clamp(theta_o, degToRad(-89.5), degToRad(89.5));
            
            //pdfM_D = 0.5 * cos(theta_o);
            pdfM_D = 1.0 / M_PI;

            //weightM /= pdfM_D;

            //Azimuthal profile
            Float phi_o = (sampler->next1D() * 2.0 - 1.0) * M_PI;//[-PI, PI]

            pdfN_D = 1.0 / (2.0 * M_PI);
            //weightN /= pdfN_D;
        }

        //We should probably checking the reason of this division more carefully 
        
        if (lobe.N_M > 1){
            weightM /= pdfM_D;
            weightN /= pdfN_D;
        }

        //Total
        weightS = weightM * weightN / pdf_lobe;
        //weightS = lobe.radiance;

        pdfS = pdfM_D * pdfN_D * pdf_lobe / cos(theta_o);
        //pdfS = pdfM_D * pdfN_D * pdf_lobe * 0.5f/ cos(theta_o);//The probability of the uniform sampling is 1/(b-a), in our case 1/2 (not helping)

        //pdfS = pdf_lobe / cos(theta_o);//Ignoring gaussian detectors
        //pdfS = pdfM_D * pdfN_D * pdf_lobe;//Solid angle instead of polar coordinates (trying stuffs)

        if (pdfS <= 0){
            printf("sampleBCSDF: Negative pdfS. pdfM_D = %f, pdfN_D = %f, pdf_lobe = %f, cos(theta_o)=%f\n", pdfM_D, pdfN_D, pdf_lobe, cos(theta_o));
        }

        //printf("sampleBCSDFNew: Sampling executed correctly\n");
        
        return 0;
    }
};

// ================ Pigmentation BSDF class ================
class Pigmentation : public BSDF {
    
    public:

    // Fur Public Methods
    Pigmentation(const Properties &props);
    Pigmentation(Stream *stream, InstanceManager *manager);
    ~Pigmentation();
    
    void configure();
    void serialize(Stream *stream, InstanceManager *manager) const;

    Float sampleH(Float h_eps, Float &pdf_h, mitsuba::Sampler* sampler) const;

    Spectrum eval(const BSDFSamplingRecord &bRec, EMeasure measure) const;
    Float pdf(const BSDFSamplingRecord &bRec, EMeasure measure) const;
    Spectrum sample(BSDFSamplingRecord &bRec, const Point2 &sample) const;
    Spectrum sample(BSDFSamplingRecord &bRec, Float &pdf, const Point2 &sample) const;
    
    //Additional functions
    static Spectrum SigmaAFromConcentration(Float ce, Float cp);
    static Spectrum SigmaAFromReflectance(const Spectrum &c, Float beta_n);
    Vector3 evalFactorizedLobes(Vector3 absorption, Float theta_i, Float theta_o, Float phi_i, Float phi_o, Float h) const;
    Vector3 evalAzimuthalLobes(Vector3 absorption, Float theta_i, Float theta_o, Float phi_i, Float phi_o, Float h) const; 
    Vector3 absorptionFromColorTexture() const;
    std::vector<Ellipse> getFiberScene(Vector3 absorption) const;
    
    Spectrum getDiffuseReflectance(const Intersection &its) const;
    Float getRoughness(const Intersection &its, int component) const;

    void addChild(const std::string &name, ConfigurableObject *child);
    std::string toString() const;

    MTS_DECLARE_CLASS()
    private:

    // Auxiliary functions
    //void toLocalFiberCoordinates(const BSDFSamplingRecord &bRec, Float &phi_i, Float &phi_o, Float &theta_i, Float &theta_o, bool debug_hair, bool use_feather_coords) const; 
    //void toLocalFiberCoordinates(const BSDFSamplingRecord &bRec, Float &phi_i, Float &theta_i, bool debug_hair, bool use_feather_coords) const;

    // General parameters
    bool active_lobes[5] = {true, true, true, false, false};//Configuration variable to the active lobes
    Float beta_m = 0.3f, beta_n = 0.3f;//Longitudinal roughness and azimuthal roughness of the cuticle respectively (stdev)
    Float absorptionFactor = 1.0f;
    int fiber_max_path = 3;//Maximum number of final lobes

    // Color parameters
    Float colorScale, colorGamma;
    ref<Texture> colorTexture;
    bool useColorTexture;

    //Pigmentation concentration
    Spectrum fur_color;//Cortex
    Vector3 fur_color2;
    int color_input_type = 0;//0: Direct absorption coeff, 1: Eu and Pheo concentrations, 2: RGB color.
    Float eu = 0.0f;// eumelanin concentration
    Float pheo = 0.0f;// pheomelanin concentration

    // Cuticle 
    Float alpha = 0.0;//Scale tilt for cuticle
    Float layers = 0.5;//Layers of the cuticle

    // Outer cylinder (cortex)
    //Float a_outer, b_outer;//Elliptical cross section parameters (this variable is actually defined in bsdf.h)
    Vector3 absorption_outer = Vector3(0.0, 0.0, 0.0);//Absorption coefficient in the cortex
    Float eta_outer = 1.55f;//Index of refraction of the cortex
    bool use_absorption_cortex = false;//Activate or not the arborption of the cortex cylinder

    // Inner cylinder (medulla)
    bool use_diffuse_medulla = false;//Activate or not the medulla cylinder
    bool use_absorption_medulla = false;//Activate or not the arborption of the cortex cylinder
    //Float a_inner, b_inner;//Medulla size parameters (defined in bsdf.h)
    Point2d medulla_center = Point2d(0.0, 0.0);//Center of the medulla
    Float medulla_ratio = 0.0f;//Medullary index (relative radius length)
    Float eta_inner = 1.55f;//Index of refraction of the medulla
    Float scattering_inner = 0.5f;//Scattering coefficient inside the medulla
    Vector3 absorption_inner = Vector3(0.0, 0.0, 0.0);//Absorption coefficient inside the medulla. [Yan et al. 2017] uses a float (grey scale value)
    Spectrum medulla_reflectance;//Ideally this value should be obtained by simulations/measured data

    // Thin Film parameters 
    // For now only outer cylinder, but it could be extended to the inner cylinder too
    bool use_thin_film = false;//Enable or disable the thin film interference
    Float thin_film_IOR;//IOR of the thin film layer (medium layer)
    Float thickness;//Layer of the thin film layer

    //Implicit path tracing variables
    std::vector<Ellipse> ellipses; 
    PigmentationSolver solver;//Auxiliary class for Fiber Tracer solver

    //Debugging variables
    mutable Float phi_min = 1000, phi_max = -1000;
    mutable Float theta_min = 1000, theta_max = -1000;
    bool verbose = false;
    bool use_default_sampling = false;//Enable uniform sphere sampling
    bool debug_hair = false;//true: Use [Yan et al. 2017] transformations
    bool use_feather_coords = false;//For feather, we do not need to transform the input coordinates as we do for hair shape scenes.
    bool use_h_range = false;//Enable custom H range obtained from the masking class 
    Float eval_max = 1000.0f;//Maximum value to clamp the BCSDF evaluation
    Float sample_max = 2.0f;//Maximum value to clamp the weights of the BCSDF sampling
    Float h_eps = 0.05f;//Epsilon of the H range [-1.0 + h_eps, 1.0 - h_eps]
    std::string fixed_path = "";//Explore a specific path for the fiber tracing (e.g. "R", "TT", "TRT"). fixed = "" explore arbitrary paths using Russian Roulette
};

MTS_IMPLEMENT_CLASS_S(Pigmentation, false, BSDF)
MTS_EXPORT_PLUGIN(Pigmentation, "Pigmentation BSDF")

MTS_NAMESPACE_END

#endif /* __FUR_H */