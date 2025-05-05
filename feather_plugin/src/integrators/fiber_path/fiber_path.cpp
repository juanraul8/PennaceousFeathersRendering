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

#include <mitsuba/render/scene.h>
#include <mitsuba/render/scenehandler.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/core/statistics.h>

MTS_NAMESPACE_BEGIN

static StatsCounter avgPathLength("Path tracer", "Average path length", EAverage);

/*! \plugin{path}{Path tracer}
 * \order{2}
 * \parameters{
 *     \parameter{maxDepth}{\Integer}{Specifies the longest path depth
 *         in the generated output image (where \code{-1} corresponds to $\infty$).
 *         A value of \code{1} will only render directly visible light sources.
 *         \code{2} will lead to single-bounce (direct-only) illumination,
 *         and so on. \default{\code{-1}}
 *     }
 *     \parameter{rrDepth}{\Integer}{Specifies the minimum path depth, after
 *        which the implementation will start to use the ``russian roulette''
 *        path termination criterion. \default{\code{5}}
 *     }
 *     \parameter{strictNormals}{\Boolean}{Be strict about potential
 *        inconsistencies involving shading normals? See the description below
 *        for details.\default{no, i.e. \code{false}}
 *     }
 *     \parameter{hideEmitters}{\Boolean}{Hide directly visible emitters?
 *        See page~\pageref{sec:hideemitters} for details.
 *        \default{no, i.e. \code{false}}
 *     }
 * }
 *
 * This integrator implements a basic path tracer and is a \emph{good default choice}
 * when there is no strong reason to prefer another method.
 *
 * To use the path tracer appropriately, it is instructive to know roughly how
 * it works: its main operation is to trace many light paths using \emph{random walks}
 * starting from the sensor. A single random walk is shown below, which entails
 * casting a ray associated with a pixel in the output image and searching for
 * the first visible intersection. A new direction is then chosen at the intersection,
 * and the ray-casting step repeats over and over again (until one of several
 * stopping criteria applies).
 * \begin{center}
 * \includegraphics[width=.7\textwidth]{images/integrator_path_figure.pdf}
 * \end{center}
 * At every intersection, the path tracer tries to create a connection to
 * the light source in an attempt to find a \emph{complete} path along which
 * light can flow from the emitter to the sensor. This of course only works
 * when there is no occluding object between the intersection and the emitter.
 *
 * This directly translates into a category of scenes where
 * a path tracer can be expected to produce reasonable results: this is the case
 * when the emitters are easily ``accessible'' by the contents of the scene. For instance,
 * an interior scene that is lit by an area light will be considerably harder
 * to render when this area light is inside a glass enclosure (which
 * effectively counts as an occluder).
 *
 * Like the \pluginref{direct} plugin, the path tracer internally relies on multiple importance
 * sampling to combine BSDF and emitter samples. The main difference in comparison
 * to the former plugin is that it considers light paths of arbitrary length to compute
 * both direct and indirect illumination.
 *
 * For good results, combine the path tracer with one of the
 * low-discrepancy sample generators (i.e. \pluginref{ldsampler},
 * \pluginref{halton}, or \pluginref{sobol}).
 *
 * \paragraph{Strict normals:}\label{sec:strictnormals}
 * Triangle meshes often rely on interpolated shading normals
 * to suppress the inherently faceted appearance of the underlying geometry. These
 * ``fake'' normals are not without problems, however. They can lead to paradoxical
 * situations where a light ray impinges on an object from a direction that is classified as ``outside''
 * according to the shading normal, and ``inside'' according to the true geometric normal.
 *
 * The \code{strictNormals}
 * parameter specifies the intended behavior when such cases arise. The default (\code{false}, i.e. ``carry on'')
 * gives precedence to information given by the shading normal and considers such light paths to be valid.
 * This can theoretically cause light ``leaks'' through boundaries, but it is not much of a problem in practice.
 *
 * When set to \code{true}, the path tracer detects inconsistencies and ignores these paths. When objects
 * are poorly tesselated, this latter option may cause them to lose a significant amount of the incident
 * radiation (or, in other words, they will look dark).
 *
 * The bidirectional integrators in Mitsuba (\pluginref{bdpt}, \pluginref{pssmlt}, \pluginref{mlt} ...)
 * implicitly have \code{strictNormals} set to \code{true}. Hence, another use of this parameter
 * is to match renderings created by these methods.
 *
 * \remarks{
 *    \item This integrator does not handle participating media
 *    \item This integrator has poor convergence properties when rendering
 *    caustics and similar effects. In this case, \pluginref{bdpt} or
 *    one of the photon mappers may be preferable.
 * }
 */
class FiberPathTracer : public MonteCarloIntegrator {
public:

    FiberPathTracer(const Properties &props)
        : MonteCarloIntegrator(props) { 

        printf("Fiber Path Tracing construction\n");

        /*
        std::string fiber_file = "/home/juanraul/Desktop/New/feathers/feather-rendering/scenes/fiber_path_tracing.xml"; 
        std::string dest_file = "test.xml";

        Properties scene_properties;

        FileResolver *resolver = Thread::getThread()->getFileResolver();
        const fs::path scenePath =
            resolver->resolveAbsolute("data/tests/test_bsdf.xml");
        
        ref<Scene> fiber_scene = new Scene(scene_properties);
        fiber_scene->addShape(NULL);*/
        
        //std::string scene_path = "/home/juanraul/Desktop/New/feathers/feather-rendering/scenes/fiber_path_tracing.xml";//It does not work

        FileResolver *resolver = Thread::getThread()->getFileResolver();
        fs::path scene_path = resolver->resolveAbsolute("data/fiber/fiber_path_tracing.xml");

        std::cout << "Processing file: "<< scene_path << std::endl;

        fiber_scene = SceneHandler::loadScene(scene_path);
        
        printf("Test\n");

        fiber_scene->initialize();

        printf("Construction done\n");
    }

    /// Unserialize from a binary data stream
    FiberPathTracer(Stream *stream, InstanceManager *manager)
        : MonteCarloIntegrator(stream, manager) { }

    Spectrum Li(const RayDifferential &r, RadianceQueryRecord &rRec) const {
        
        //Setting the fiber scene
        ref<Scene> scene = fiber_scene;

        //Replace the original scene with the fiber scene
        RadianceQueryRecord new_rRec = RadianceQueryRecord(rRec);
        new_rRec.scene = scene;

        //new_rRec = rRec;//Original scene
        //rRec = RadianceQueryRecord(new_rRec);

        /* Some aliases and local variables */
        //const Scene *scene = rRec.scene;//Original code
        Intersection &its = new_rRec.its;
        RayDifferential ray(r);
        Spectrum Li(0.0f);
        bool scattered = false;

        //printf(scene->m_shapes.size());
        //printf("Number of shapes: %d\n", scene->getShapes().size());

        /* Perform the first ray intersection (or ignore if the
           intersection has already been provided). */
        new_rRec.rayIntersect(ray);
        ray.mint = Epsilon;

        Spectrum throughput(1.0f);
        Float eta = 1.0f;

        //printf("Valid intersection: %d\n", its.isValid());

        while (new_rRec.depth <= m_maxDepth || m_maxDepth < 0) {
            
            //printf("Fiber Path Tracing depth = %d\n", new_rRec.depth);

            if (!its.isValid()) {
                /* No intersection could be found */

                //printf("No intersection!\n");

                break;
            }

            const BSDF *bsdf = its.getBSDF(ray);

            Spectrum test = bsdf->getDiffuseReflectance(its);
            //printf("BSDF diffuse reflectance: (%f, %f, %f)\n", test[0], test[1], test[2]);

            /* Include radiance from a subsurface scattering model if requested */
            if (its.hasSubsurface() && (rRec.type & RadianceQueryRecord::ESubsurfaceRadiance))
                Li += throughput * its.LoSub(scene, new_rRec.sampler, -ray.d, new_rRec.depth);

            if ((new_rRec.depth >= m_maxDepth && m_maxDepth > 0)
                || (m_strictNormals && dot(ray.d, its.geoFrame.n)
                    * Frame::cosTheta(its.wi) >= 0)) {

                /* Only continue if:
                   1. The current path length is below the specifed maximum
                   2. If 'strictNormals'=true, when the geometric and shading
                      normals classify the incident direction to the same side */
                break;
            }

            /* ==================================================================== */
            /*                            BSDF sampling                             */
            /* ==================================================================== */
            //printf("BSDF sampling!\n");

            /* Sample BSDF * cos(theta) */
            Float bsdfPdf;
            BSDFSamplingRecord bRec(its, new_rRec.sampler, ERadiance);
            Spectrum bsdfWeight = bsdf->sample(bRec, bsdfPdf, rRec.nextSample2D());
            if (bsdfWeight.isZero())
                break;

            //printf("Test!\n");
            scattered |= bRec.sampledType != BSDF::ENull;

            /* Prevent light leaks due to the use of shading normals */
            const Vector wo = its.toWorld(bRec.wo);
            Float woDotGeoN = dot(its.geoFrame.n, wo);
            if (m_strictNormals && woDotGeoN * Frame::cosTheta(bRec.wo) <= 0)
                break;

            /* Evaluate BSDF * cos(theta) */
            const Spectrum bsdfVal = bsdf->eval(bRec);
            Li += throughput * bsdfVal;
            //Li += bsdfVal;

            //printf("Trace ray!\n");

            /* Trace a ray in this direction */
            ray = Ray(its.p, wo, ray.time);
            
            //printf("Ray Intersection? %d\n", scene->rayIntersect(ray, its));

            //Here we should save the final rays I guess? 
            if (!scene->rayIntersect(ray, its)) {
                
                //printf("No intersection!\n");
                break;
            }

            /* Keep track of the throughput and relative refractive index along the path */
            throughput *= bsdfWeight;
            eta *= bRec.eta;

            /* ==================================================================== */
            /*                         Indirect illumination                        */
            /* ==================================================================== */
            //printf("Indirect Illumiations!\n");


            /* Set the recursive query type. Stop if no surface was hit by the
               BSDF sample or if indirect illumination was not requested */
            if (!its.isValid() || !(new_rRec.type & RadianceQueryRecord::EIndirectSurfaceRadiance))
                break;
            
            new_rRec.type = RadianceQueryRecord::ERadianceNoEmission;

            if (new_rRec.depth++ >= m_rrDepth) {
                /* Russian roulette: try to keep path weights equal to one,
                   while accounting for the solid angle compression at refractive
                   index boundaries. Stop with at least some probability to avoid
                   getting stuck (e.g. due to total internal reflection) */

                Float q = std::min(throughput.max() * eta * eta, (Float) 0.95f);
                if (new_rRec.nextSample1D() >= q)
                    break;
                throughput /= q;
            }
        }

        /* Store statistics */
        avgPathLength.incrementBase();
        avgPathLength += new_rRec.depth;

        //printf("Li = (%f, %f, %f)\n", Li[0], Li[1], Li[2]);

        return Li;
    }

    inline Float miWeight(Float pdfA, Float pdfB) const {
        pdfA *= pdfA;
        pdfB *= pdfB;
        return pdfA / (pdfA + pdfB);
    }

    void serialize(Stream *stream, InstanceManager *manager) const {
        MonteCarloIntegrator::serialize(stream, manager);
    }

    std::string toString() const {
        std::ostringstream oss;
        oss << "FiberPathTracer[" << endl
            << "  maxDepth = " << m_maxDepth << "," << endl
            << "  rrDepth = " << m_rrDepth << "," << endl
            << "  strictNormals = " << m_strictNormals << endl
            << "]";
        return oss.str();
    }

    ref<Scene> fiber_scene;//Custom fiber scene

    MTS_DECLARE_CLASS()
};

//NB: Avoid MTS_EXPORT_PLUGIN if you use this path tracer with the pigmentation classes
MTS_IMPLEMENT_CLASS_S(FiberPathTracer, false, MonteCarloIntegrator)
//MTS_EXPORT_PLUGIN(FiberPathTracer, "Fiber path tracer");
MTS_NAMESPACE_END
