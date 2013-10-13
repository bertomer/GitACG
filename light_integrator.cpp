/*******************************************************************************
 *  light_integrator.cpp
 *******************************************************************************
 *  Copyright (c) 2013 Alexandre Kaspar <alexandre.kaspar@a3.epfl.ch>
 *  For Advanced Computer Graphics, at the LGG / EPFL
 *
 *        DO NOT REDISTRIBUTE
 ***********************************************/

#include <nori/acg.h>
#include <nori/bsdf.h>
#include <nori/common.h>
#include <nori/integrator.h>
#include <nori/luminaire.h>
#include <nori/sampler.h>
#include <nori/scene.h>
#include <vector>
#include <time.h>

NORI_NAMESPACE_BEGIN

// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// put your group number here!
#define GROUP_NUMBER 13
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

GROUP_NAMESPACE_BEGIN()

/**
 * \brief Simple local illumination integrator
 * using light area sampling
 */
class LightIntegrator : public Integrator {
public:

        LightIntegrator(const PropertyList &propList) {
                Q_UNUSED(propList);
                srand(time(NULL));
        }

        /// Return the mesh corresponding to a given luminaire
        inline const Mesh *getMesh(const Luminaire *lum) const {
                const Mesh *mesh = dynamic_cast<const Mesh *> (lum->getParent());
                if (!mesh) throw NoriException("Unhandled type of luminaire!");
                return mesh;
        }

        /**
         * \brief Directly sample the lights, providing a sample weighted by 1/pdf
         * where pdf is the probability of sampling that given sample
         *
         * \param scene
         * the scene to work with
         *
         * \param lRec
         * the luminaire information storage
         *
         * \param _sample
         * the 2d uniform sample
         *
         * \return the sampled light radiance including its geometric, visibility and pdf weights
         */
        inline Color3f sampleLights(const Scene *scene, LuminaireQueryRecord &lRec, const Point2f &_sample) const {
                Point2f sample(_sample);
                const std::vector<Luminaire *> &luminaires = scene->getLuminaires();

                if (luminaires.size() == 0)
                        throw NoriException("LightIntegrator::sampleLights(): No luminaires were defined!");

                // TODO Implement the following steps
                // and take care of using the good G, V terms to work with the Li method below






                // 3. Compute geometry term G and visibility term on the luminaire's side (no cos(w) of the mesh side)
                // as well as the pdf of that point being found
                // use Mesh::pdf to get the probability of choosing the point in Mesh::samplePosition

                // 4. Return radiance emitted from luminaire multiplied by the appropriate terms G, V ...

                // 1. Choose one luminaire at random
                int randomIdx = rand() % luminaires.size();
                lRec.luminaire = luminaires.at(randomIdx);

                // 2. Sample the position on the luminaire mesh
                // using Mesh::samplePosition(const Point2d &sample, Point3f &p, Normal3f &n)
                Mesh *m = (Mesh*) luminaires.at(randomIdx)->getParent();

                // lRec.p is the sample point on the mesh
                m->samplePosition(sample, lRec.p, lRec.n);
                lRec.pdf = m->pdf();

                // lRec.d = direction vector from ref to p
                lRec.d = lRec.p - lRec.ref;
                lRec.dist = lRec.d.norm();
                lRec.d /= lRec.dist;

                Ray3f shadowRay(lRec.ref, lRec.d);

                float bigV=0;
                Intersection its;
                // are we hiting the an obstacle or the luminaire?
                if(scene->rayIntersect(shadowRay, its)){
                    if(its.mesh == m) {
                        bigV = 1;
                    } else {
                        bigV = 0;
                    }
                }

                // cos theta prime prime :) (luminaire's side)
                float cosThetaI = std::max(0.0f, lRec.n.dot(-lRec.d));

                // return V * cos theta'' * fr(x", x', x) / (pdf(x") * dist²)
                return bigV*cosThetaI*lRec.luminaire->eval(lRec) / ((lRec.dist * lRec.dist) * lRec.pdf);
        }

        /**
         * \brief Simple local illumination integration:
         * We cast a ray from the camera, intersects it with the first element
         * in the scene it goes through and from there we directly sample the
         * light's in the scene to compute direct lighting.
         */
        Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray_) const {
                Ray3f ray(ray_);

                /* Find the surface that is visible in the requested direction */
                Intersection its;
                if (!scene->rayIntersect(ray, its))
                        return Color3f(0.0f);

                const Mesh *mesh = its.mesh;
                const BSDF *bsdf = mesh->getBSDF();

                // TODO implement direct lighting using light sampling using
                //      sampleLights(const Scene *, LuminaireQueryRecord &, const Point2d &)
                // which you also have to implement

                /* If we hit a luminaire, use its related color information */

                if (mesh->isLuminaire()) {
                        const Luminaire *luminaire = its.mesh->getLuminaire();
                        LuminaireQueryRecord lRec(luminaire, ray.o, its.p, its.shFrame.n);
                        return luminaire->eval(lRec);
                }

                LuminaireQueryRecord lRec;
                lRec.ref = its.p;
                Color3f sl = sampleLights(scene, lRec, sampler->next2D());

                BSDFQueryRecord bRec(its.toLocal(-ray.d), its.toLocal(lRec.d), ESolidAngle);
                // TODO: bsdf->eval(bRec)...
                // cos theta prime (mesh's side)
                float cosTheta = std::max(0.0f, its.shFrame.n.dot(lRec.d));

                return sl*bsdf->eval(bRec)*cosTheta;
        }

        QString toString() const {
                return QString("LightIntegrator[]");
        }
};

GROUP_NAMESPACE_END

NORI_REGISTER_GROUP_CLASS(LightIntegrator, "light")
NORI_NAMESPACE_END
