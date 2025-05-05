#pragma once
#include <mitsuba/core/warp.h>
#include <mitsuba/core/vector.h>
#include <mitsuba/core/plugin.h>
#include <typeinfo>
#include <queue>

#define INV_FOURPI 0.07957747154594766788f

// This is uggly I know, but it works for now
// Adrian Jarabo: Just typing "using namespace mitsuba" should do the exact same
typedef mitsuba::Spectrum Spectrum;

typedef mitsuba::Vector2d Vector2d;
typedef mitsuba::Point2d Point2d;

typedef mitsuba::Vector3 Vector3;
typedef mitsuba::Vector3f Vector3f;
typedef mitsuba::Vector3d Vector3d;

typedef mitsuba::Point3d Point3d;
typedef mitsuba::Point3f Point3f;

typedef mitsuba::Float Float;

//---------------------------------------Math util functions!
Float eta_min2 = 100, eta_max2 = -1; 

inline Vector2d normalizeSafe(const Vector2d &v){
    if (v.length() == 0.0)
        return v;
    else
        return mitsuba::normalize(v);
}

Vector3d normalizeSafe(const Vector3d &v){
    if (v.length() == 0.0)
        return v;
    else
        return mitsuba::normalize(v);
}

Vector3f normalizeSafe(const Vector3f &v){
    if (v.length() == 0.0)
        return v;
    else
        return mitsuba::normalize(v);
}

double ourDot(Vector2d a, Vector2d b){

    return a.x*b.x + a.y * b.y; 
}

double clamp(double x, double a, double b) {
    if (x < a) return a;
    if (x > b) return b;

    return x;
}

inline double acosSafe(double cos_theta) {
    double eps = 1e-12;

    return acos(clamp(cos_theta, -1.0 + eps, 1.0 - eps));
}

inline Vector3 exp3(const Vector3 &v) {
    return Vector3(exp(v[0]), exp(v[1]), exp(v[2]));
}

inline Float acosSafe(Float cosTheta) {
    Float eps = 1e-12;

    return acos((Float) clamp(cosTheta, Float(-1.0 + eps), Float(1.0 - eps)));
}

inline Float asinSafe(Float sinTheta) {
    Float eps = 1e-12;

    return asin((Float) clamp(sinTheta, Float(-1.0 + eps), Float(1.0 - eps)));
}

struct OurRay{
    Point2d o;//Origin of the ray        
    Vector2d d;//Direction of the ray
    std::string path = "";//name of the path
    double fresnel = 1.0;//Accumulate fresnel value (radiance)
    double dist_cortex = 0.0;//Distance travelled inside the cortex
    double dist_medulla = 0.0;//Distance travelled inside the medulla
    int N_R = 0;//The number of internal reflections

    //Auxiliary variables
    int active_idx = -1;
    int prev_idx = -1;
};

// Lobes results

struct Lobes {

    OurRay R_ray;
    OurRay TT_ray;
    OurRay TRT_ray;

    std::vector<OurRay> rays;//Rays collected during the ray tracing
    std::vector<Vector3> N;//Azimuthal Scattering for unscattered lobes
    std::vector<double> D;//Azimuthal Distribution term for unscattered lobes
    std::vector<double> M;//Longitudinal Distribution term for unscattered lobes
    std::vector<Vector3> A;//Attenuation term for unscattered lobes
    std::vector<double> Phi;//Exit angle by Path Tracing
    std::vector<double> fresnel;//Final fresnel term for unscattered lobes
};  


//Ellipse struct
struct Ellipse {             
    
    //Geometry
    double a = 1.0;//Mayor axis        
    double b = 1.0;//Minor axis
    Point2d center = Point2d(0.0, 0.0);//Center of the ellipse

    //Materials
    double eta_exterior = 1.0;//Exterior index of refraction
    double eta_interior = 1.0;//Interior index of refraction
    bool specular_mat = false;//Is it the ellipse a dieletric material?

    Vector3f absorption = Vector3f(0.0f);//Absorption coefficient
    bool absorption_mat = false;//Has the ellipse an absorbing media?

    Vector3f diffuse_reflectance = Vector3f(1.0f);//Diffuse reflectance
    bool diffuse_mat;//Is it the ellipse a lambertian material?
    
    bool thin_film_mat = false;//There is thin film interference for this material?
    Float thickness = 300.0;//Thickness of the thin film layer
    Float thin_film_IOR = 1.33;//IOR of the medium layer (thin film layer)

    Ellipse() {
    
        this->a = 1.0;
        this->b = 1.0;
        this->center = Point2d(0.0, 0.0);
        
        this->eta_exterior = 1.0;
        this->eta_interior = 1.0;
        this->absorption = Vector3f(0.0f);
        this->diffuse_reflectance = Vector3f(1.0f);

        absorption_mat = false;
        diffuse_mat = false;
        specular_mat = false;
        thin_film_mat = false;

        this->thickness = 300.0;
        this->thin_film_IOR = 1.33;
    }   

    Ellipse(double a, double b, double eta_exterior, double eta_interior) {
    
        this->a = a;
        this->b = b;
        this->center = Point2d(0.0, 0.0);
        this->eta_exterior = eta_exterior;
        this->eta_interior = eta_interior;

        this->absorption = Vector3f(0.0f);
        this->diffuse_reflectance = Vector3f(1.0f);

        absorption_mat = false;
        diffuse_mat = false;
        specular_mat = false;
        thin_film_mat = false;
    }

    Float computeEccentricity(){    
        Float e = 0.0f;

        if (a > b){
            e = sqrt(1.0 - (b*b)/(a*a));
        }else{
            e = sqrt(1.0 - (a*a)/(b*b));
        }

        printf("Ellipse (eccentricity) a = %f, b = %f, e = %f\n", a, b, e);

        return e;
    }

    void print(){
        printf("Properties: absorption (%d), diffuse (%d), specular (%d), thin film interference (%d)\n", absorption_mat, diffuse_mat, specular_mat, thin_film_mat);
        printf("a = %f, b = %f, eta exterior = %f, eta interior = %f, center (%f, %f)\n Absorption coeff (%f, %f, %f) \n", a, b, eta_exterior, eta_interior, 
              center.x, center.y, absorption.x, absorption.y, absorption.z);
    }   

    // Compute ellipse and ray intersection
    double rayIntersect(OurRay &ray){
        // NB: it can be similar to sphere.cpp ray intersection but with the ellipses equation (mitsuba/src/shapes/sphere.cpp)
        // Numerical issues can arise depending on how -mint- and -maxt- are selected
        // if t is inf --> no intersection!

        //Numerical setup
        double inf = std::numeric_limits<double>::infinity();
        double ray_eps = 1e-3;//1e-6
        double inner_eps = 1e-12;//1e-12

        Point2d o = Point2d(ray.o.x, ray.o.y);
        Vector2d d = Vector2d(ray.d.x, ray.d.y);

        d = normalizeSafe(d);

        //Quadtratic equation for Ellipse-Ray Intersection
        double A = (d.x*d.x) / (a * a) + (d.y*d.y) / (b * b);
        double B = (2.0 * d.x * (o.x - center.x)) / (a * a) + (2.0 * d.y * (o.y - center.y)) / (b * b);
        double C = (o.x - center.x) * (o.x - center.x) / (a * a) + (o.y - center.y) * (o.y - center.y) / (b * b) - 1.0;

        double inner = B * B - 4.0 * A * C;

        // Avoiding numerical issues

        //printf("Ray Intersection (inner) = %.18f\n", inner);
        //inner = 0.0;

        if (abs(inner) < inner_eps){
            inner = 0.0;
        }

        //No solutions
        if (inner < 0.0){

            //printf("Ray Intersection: No solution\n");
            return inf;
        }

        double t0 = (-B - std::sqrt(inner)) / (2.0 * A);
        double t1 = (-B + std::sqrt(inner)) / (2.0 * A);

        if (t0 <= 0) t0 = inf;
        if (t1 <= 0) t1 = inf;

        // Pick nearest intersection
        double t = t0 < t1 ? t0 : t1;

        return t;
    }

    // Compute ellipse and ray intersection
    double rayIntersect(mitsuba::Ray &ray){
        //NB: it can be similar to shere.cpp ray intersection but with the ellipses equation (mitsuba/src/shapes/sphere.cpp)
        //Numerical issues can arise depending on how -mint- and -maxt- are selected
        // if t is inf --> no intersection!

        //Numerical setup
        double inf = std::numeric_limits<double>::infinity();
        double ray_eps = 1e-3;//1e-6
        double inner_eps = 1e-12;//1e-12

        // Cross section information for the ray-ellipse intersection
        Point2d o = Point2d(ray.o.x, ray.o.y);
        Vector2d d = Vector2d(ray.d.x, ray.d.y);

        d = normalizeSafe(d);

        //Quadtratic equation for Ellipse-Ray Intersection
        double A = (d.x*d.x) / (a * a) + (d.y*d.y) / (b * b);
        double B = (2.0 * d.x * (o.x - center.x)) / (a * a) + (2.0 * d.y * (o.y - center.y)) / (b * b);
        double C = (o.x - center.x) * (o.x - center.x) / (a * a) + (o.y - center.y) * (o.y - center.y) / (b * b) - 1.0;

        double inner = B * B - 4.0 * A * C;

        // Avoiding numerical issues

        //printf("Ray Intersection (inner) = %.18f\n", inner);
        //inner = 0.0;

        if (abs(inner) < inner_eps){
            inner = 0.0;
        }

        //No solutions
        if (inner < 0.0){

            //printf("Ray Intersection: No solution\n");
            return inf;
        }

        double t0 = (-B - std::sqrt(inner)) / (2.0 * A);
        double t1 = (-B + std::sqrt(inner)) / (2.0 * A);

        if (t0 < 0) t0 = inf;
        if (t1 < 0) t1 = inf;

        // Pick nearest intersection
        double t = t0 < t1 ? t0 : t1;

        // return t;//2D case
        
        // 3D case for an infinite cylinder needs to consider the inclination of the ray!
        
        Float rho = sqrt(ray.d.x * ray.d.x + ray.d.y * ray.d.y);
        Float theta = atan2(ray.d.z, rho); //elevation interpretation

        // return t/cos(theta);//3D case

        //Adrian Jarabo change
    
        if(abs(ray.d.z) > (1.0 - 1.e-8))
            return inf;

        return t / sqrt(1.0 - ray.d.z*ray.d.z);
    }

    //Tangent vector at the contact point p based on implicit differentiation
    //More details about implicit differentiation: https://www.cuemath.com/calculus/implicit-differentiation/
    Vector2d getTangent(Point2d p){

        Vector2d tangent;

        //Special cases (avoid division by 0)
        //NB: check this carefully later!!!)

        //Vertical tangent
        if (p.y == 0.0){
            
            //printf("Special case: p = (%f, %f)\n", p.x, p.y);

            double sign_x = 1.0;

            if (p.x < 0.0) sign_x = -1.0;

            return Vector2d(0.0, sign_x);
        }

        //Derivate by implicit differentiation
        double c = (b * b) / (a * a);
        double m = -(b * b) / (a * a) * (p.x - center.x) / (p.y - center.y);

        double x = -(a+b);

        if (p.y < 0.0) x = -x;

        double y = m * (x - p.x) + p.y;

        tangent = Vector2d(x - p.x, y - p.y);//Vector between two points: contact point and the derived one
        tangent = normalizeSafe(tangent);

        return tangent;
    }

    Vector2d getNormal(Point2d p){
        //NB: We use the point form from the following reference (https://mathemerize.com/equation-of-normal-to-ellipse-in-all-forms/) 

        //Equation of normal in terms of point of contact --> n = normalizeSafe(p - self.center)

        //Vector2d normal = Vector2d(p.x - center.x, p.y - center.y);
        //normal = normalizeSafe(normal);

        //Implicit Differentiation
        //Another way: https://www.toppr.com/ask/content/concept/normal-to-an-ellipse-207451/ 

        Vector2d tangent = getTangent(p);
        Vector2d normal = Vector2d(tangent.y, -tangent.x);//Rotate tangent vector by -90º (Perpendicular vector)

        //Float length = sqrt(normal.x * normal.x + normal.y * normal.y);
        //printf("Length normal = %f\n", length);

        return normal;
    }

    Vector3f getTangent(Point3f p){

        Vector3f tangent;

        //Special cases (avoid division by 0)
        //NB: check this carefully later!!!)

        //Vertical tangent
        if (p.y == 0.0){
            
            //printf("Special case: p = (%f, %f)\n", p.x, p.y);

            double sign_x = 1.0;

            if (p.x < 0.0) sign_x = -1.0;

            return Vector3f(0.0, sign_x, 0.0);
        }

        //Derivate by implicit differentiation
        double c = (b * b) / (a * a);
        double m = -(b * b) / (a * a) * (p.x - center.x) / (p.y - center.y);

        double x = -(a+b);

        if (p.y < 0.0) x = -x;

        double y = m * (x - p.x) + p.y;

        tangent = Vector3f(x - p.x, y - p.y, 0.0);//Vector between two points: contact point and the derived one
        tangent = normalizeSafe(tangent);

        return tangent;
    }

    Vector3f getNormal(Point3f p){
        //NB: We use the point form from the following reference (https://mathemerize.com/equation-of-normal-to-ellipse-in-all-forms/) 
        //printf("Using 3d normal\n");

        //Implicit Differentiation
        //Another way: https://www.toppr.com/ask/content/concept/normal-to-an-ellipse-207451/ 

        /*
        Point2d p2 = Point2d(p.x, p.y);

        Vector2d tangent = getTangent(p2);
        Vector3f normal = Vector3f(tangent.y, -tangent.x, 0.0f);//Rotate tangent vector by -90º (Perpendicular vector)
        normal = normalizeSafe(normal);
        */

        // Flatland (Adolfo Munoz)
        Vector3f normal = Vector3f(2.0f * (p.x-center.x)/(this->a * this-> a), 2.0f * (p.y-center.y)/(this->b * this -> b), 0.0f);
        normal = normalizeSafe(normal);

        //Float length = sqrt(normal.x * normal.x + normal.y * normal.y);
        //printf("Length normal = %f\n", length);

        return normal;
    }

}; 

MTS_NAMESPACE_BEGIN
// --------------------Our coordinate system transformation

// Previous transformation 
inline Vector2 vector2polar(Vector m_vector){
    float theta = M_PI / 2.0f - acosSafe(m_vector.x);//This should be equivalent to call asin. 
    Float phi = atan2(m_vector.z, m_vector.y);  
    return Vector2(theta, phi);
}

// New coordinate system: important to align the fiber along the proper axis
// In practice this coordinate system transformation is equivalent to flipNormal=true inside the feather underlying geometry (curved plane)
/*
inline Vector2 vector2polar(Vector m_vector){
    Float theta = asinSafe(m_vector.y); 
    Float phi = atan2(m_vector.z, m_vector.x);  
    
    return Vector2(theta, phi);
}
*/

// Angles from outgoing ray of the lobe (old)
inline Vector2 FiberTracer2BCSDF(Vector m_vector){
    Float theta = asinSafe(m_vector.z);//elevation interpretation. AsinSafe avoid numerical issues
    Float phi = atan2(m_vector.y, m_vector.x);  
    
    return Vector2(theta, phi);
}

inline Vector getCorrectFrame(const BSDFSamplingRecord &bRec, Vector input_vector, bool use_feather_coords){
    
    if (use_feather_coords){
        return input_vector;//Feather BSDF do not require this transformation
    }
    
    // The code of Yan et al. transform to World Coords and from that it calculates the local coords sys. 
    // In their code, they use the geoFrame, however by default Mitsuba transform to local coordinate using the shadingFrame,
    // which is not the same to the geoFrame for the hair geometry. Therefore, we need to first go from local to world using the 
    // shading frame, then going from World2Local using the geoFrame. In this way we can compare with their results.  
    const Frame geoFrame = bRec.its.geoFrame;
    const Frame shFrame = bRec.its.shFrame;
    Vector tmpvector = shFrame.toWorld(input_vector);
    return geoFrame.toLocal(tmpvector);//Hair shape works only with the geometric frame
    
}

inline Vector4 getFiberAngles(const BSDFSamplingRecord &bRec){

    Vector3 wi = bRec.its.toWorld(bRec.wi);
    Vector3 wo = bRec.its.toWorld(bRec.wo);

    //Vector3 wi = bRec.wi;
    //Vector3 wo = bRec.wo;
   
    //For scenes using the hair intersection shape, it only works the geoFrame. Otherwise, both geoFrame and shFrame are fine
    const Vector3 &u = bRec.its.geoFrame.s;
    const Vector3 &v = bRec.its.geoFrame.t;
    const Vector3 &w = bRec.its.geoFrame.n;

    /*
    const Vector3 &u = bRec.its.shFrame.s;
    const Vector3 &v = bRec.its.shFrame.t;
    const Vector3 &w = bRec.its.shFrame.n;
    */

    /*
    printf ("u (%f, %f, %f)\n", u.x, u.y, u.z);
    printf ("v (%f, %f, %f)\n", v.x, v.y, v.z);
    printf ("w (%f, %f, %f)\n", w.x, w.y, w.z);
    */

    Float theta_i = M_PI / 2.0f - acosSafe(dot(u, wi));
    Float theta_r = M_PI / 2.0f - acosSafe(dot(u, wo));

    //Float theta_i = asinSafe(dot(u, wi));
    //Float theta_r = asinSafe(dot(u, wo));

    //NB: Clamp to avoid numerical issues
    // theta_i = clamp(theta_i, degToRad(-89.5), degToRad(89.5));
    // theta_r = clamp(theta_r, degToRad(-89.5), degToRad(89.5));

    Float phi_i = acosSafe(dot(normalizeSafe(wi - dot(u, wi) * u), v)) * (dot(w, wi) > 0 ? 1 : -1);
    Float phi_r = acosSafe(dot(normalizeSafe(wo - dot(u, wo) * u), v)) * (dot(w, wo) > 0 ? 1 : -1);

    theta_i = clamp(phi_i, degToRad(-89.5), degToRad(89.5));
    theta_r = clamp(phi_r, degToRad(-89.5), degToRad(89.5));

    //phi_i = phi_i + M_PI;
    //phi_r = phi_r + M_PI;

    return Vector4(theta_i,phi_i, theta_r, phi_r);
}

void inline toLocalFiberCoordinates(const BSDFSamplingRecord &bRec, Float &phi_i, Float &phi_o, Float &theta_i, Float &theta_o, bool debug_hair, bool use_feather_coords)
{
    // [Yan et al. 2017] transformation
    if (debug_hair)
    {
        Vector4 angles = getFiberAngles(bRec);
        theta_i = angles[0];
        phi_i = angles[1];
        theta_o = angles[2];
        phi_o = angles[3];
    }
    else
    {
        // Our transformation
        Vector2 param_angle = vector2polar(getCorrectFrame(bRec, bRec.wi, use_feather_coords));
        theta_i = param_angle[0];
        phi_i = param_angle[1];

        param_angle = vector2polar(getCorrectFrame(bRec, bRec.wo, use_feather_coords));
        theta_o = param_angle[0];
        phi_o = param_angle[1];

        if (!use_feather_coords){
            Float m_eps = 1.e-2;

            //Validate again Yan et al. 2017 transformations. (This makes sense only for the furball scene!)
            Vector4 angles = getFiberAngles(bRec);
            Float theta_i_r = angles[0];
            Float phi_i_r = angles[1];
            Float theta_o_r = angles[2];
            Float phi_o_r = angles[3];
            
            if ( theta_i_r - theta_i > m_eps || phi_i_r - phi_i > m_eps){
                printf(">> Input (Yan et al.) [%f %f] vs Input (Ours) [%f %f] %s \n",theta_i_r,phi_i_r,theta_i,phi_i,bRec.wi.toString().c_str());
            }
            
            if ( (theta_o_r - theta_o > m_eps || phi_o_r - phi_o > m_eps) && (bRec.wo[0] > Epsilon || bRec.wo[1] > Epsilon || bRec.wo[2] > Epsilon) ){
                printf(">> Output (Yan et al.) [%f %f] vs Output (Ours) [%f %f] %s \n",theta_o_r,phi_o_r,theta_o,phi_o,bRec.wo.toString().c_str());
            }
        }
    }

    //Avoiding numerical issues at grazing angles
    theta_i = clamp(theta_i, degToRad(-89.5), degToRad(89.5));
    theta_o = clamp(theta_o, degToRad(-89.5), degToRad(89.5));
}

void inline toLocalFiberCoordinates(const BSDFSamplingRecord &bRec, Float &phi_i, Float &theta_i, bool debug_hair, bool use_feather_coords)
{
    // [Yan et al. 2017] transformation
    if (debug_hair)
    {
        Vector4 angles = getFiberAngles(bRec);
        theta_i = angles[0];
        phi_i = angles[1];
    }
    else
    {
        // Our transformation
        Vector2 param_angle = vector2polar(getCorrectFrame(bRec, bRec.wi, use_feather_coords));
        theta_i = param_angle[0];
        phi_i = param_angle[1];

        Float m_eps = 1.e-2;
        if (!use_feather_coords){
            //Validate again Yan et al. 2017 transformations. (This makes sense only for the furball scene!)
            Vector4 angles = getFiberAngles(bRec);
            Float theta_i_r = angles[0];
            Float phi_i_r = angles[1];
            Float theta_o_r = angles[2];
            Float phi_o_r = angles[3];
            if ( theta_i_r - theta_i > m_eps || phi_i_r - phi_i > m_eps){
                printf(">> Input (Yan et al.) [%f %f] vs Input (Ours) [%f %f] %s \n",theta_i_r,phi_i_r,theta_i,phi_i,bRec.wi.toString().c_str());
            }
        }
    }

    theta_i = clamp(theta_i, degToRad(-89.5), degToRad(89.5));
}

inline Vector getVectorFromParameterizedAngle(BSDFSamplingRecord bRec, Float theta_o, Float phi_o){
    
    Float ul = sin(theta_o);    
    Float vl = cos(theta_o) * sin(phi_o);
    Float wl = cos(theta_o) * cos(phi_o);

    //return Vector(ul, vl, wl);//Adrian Jarabo (it works well for furball)
    return Vector(ul,vl,wl);//Dario Lanza
}

// Dario Lanza: New parametrization is important to locate the fiber along the proper axis
/*
inline Vector getVectorFromParameterizedAngle(BSDFSamplingRecord bRec, Float theta_o, Float phi_o){
    Float ul = cos(theta_o) * cos(phi_o);     
    Float vl = sin(theta_o) * cos(phi_o);
    Float wl = sin(theta_o);
    
    return Vector(ul, vl, wl);
}
*/

//Juan Raul: this function is not using the parameter h. This parameter is key to parametrize the fiber!
inline static Ray getInitialRayTest(const std::vector<Ellipse> &ellipses, Float phi_i, Float theta_i, Float h)
    {
        Ray init_ray;

        // The outer ellipse is used to define the fiber azimuthal radius parametrization (h)
        Float outer_a = ellipses[0].a;
        Float outer_b = ellipses[0].b;

        // Adrian Jarabo: Just assume that the ray sampling is centered on the outermost ellipse; that is the obvious
        // choice for an outer ellipse bounding all others.
        Point2d center = ellipses[0].center;
        
        // Float offset = (outer_a + outer_b) / cos(theta_i); // Offset to make sure the ray is in fact outside of the elliptical cylinder
        Float offset = abs( (outer_a + outer_b) / cos(phi_i)); // Offset to make sure the ray is in fact outside of the elliptical cylinder
        // Float offset = abs( (outer_a + outer_b)); // Offset to make sure the ray is in fact outside of the elliptical cylinder

        // Spherical coordinates (Linqi/Marschner scenario, elevation interpretation)
        // I think the negative d is wrong (see the python path tracing viewer)
        // init_ray.d = -Vector3f(cos(phi_i) * cos(theta_i), sin(phi_i) * cos(theta_i), sin(theta_i));

        // init_ray.d = -Vector3f(cos(phi_i) * cos(theta_i), sin(phi_i) * cos(theta_i), sin(theta_i));
        //D.L coordinate frame sys
        // init_ray.d = -Vector3f(cos(phi_i) * cos(theta_i), sin(phi_i), sin(theta_i)*cos(phi_i));
        init_ray.d = -Vector3f(cos(phi_i)*cos(theta_i), sin(phi_i) * cos(theta_i), sin(theta_i));
        
        // init_ray.d = Vector3f(sin(phi_i) * cos(theta_i), cos(phi_i) * cos(theta_i), sin(theta_i));
        //init_ray.d = Vector3f(cos(phi_i) * cos(theta_i), sin(phi_i) * cos(theta_i), sin(theta_i));

        // Setting the origin with respect to the h (fiber radius)
        // TODO - Fix this
        // init_ray.o = Point3f(-h * outer_a * sin(phi_i) * cos(theta_i) + center.x, h * outer_b * cos(phi_i) * cos(theta_i) + center.y, sin(theta_i));
        // init_ray.o = Point3f(h * outer_a * cos(theta_i) * cos(phi_i) + center.x, h * outer_b * sin(theta_i)*cos(phi_i) + center.y, sin(phi_i));
        
        // use a radius for all 3 coords otherwise we break the polar coord sys
        Float maj_axis = outer_a > outer_b ? outer_a : outer_b;
        Float r =  maj_axis;
        init_ray.o = Point3f( r * cos(phi_i) * cos(theta_i) + center.x, r*sin(phi_i)*cos(theta_i) + center.y, r *sin(theta_i));

        // Setting the ray outside of the cross section.
        init_ray.o = init_ray.o - offset * init_ray.d;

        // init_ray.o = Point3f(-h * outer_a * sin(phi_i) + center.x, h * outer_b * cos(phi_i), 0.);
        // init_ray.o = init_ray.o - (outer_a + outer_b) * init_ray.d;

        return init_ray;
    }


MTS_NAMESPACE_END