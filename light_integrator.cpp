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

                // 1. Choose one luminaire at random


                // 2. Sample the position on the luminaire mesh
                // using Mesh::samplePosition(const Point2d &sample, Point3f &p, Normal3f &n)

                // 3. Compute geometry term G and visibility term on the luminaire's side (no cos(w) of the mesh side)
                // as well as the pdf of that point being found
                // use Mesh::pdf to get the probability of choosing the point in Mesh::samplePosition

                // 4. Return radiance emitted from luminaire multiplied by the appropriate terms G, V ...

                int randomIdx = rand() % luminaires.size();
                lRec.luminaire = luminaires.at(randomIdx);

                Mesh *m = (Mesh*) luminaires.at(randomIdx)->getParent();

                // lRec.p is the sample point on the mesh
                m->samplePosition(sample, lRec.p, lRec.n);
                lRec.pdf = m->pdf();

                // lRec.d = direction vector from ref to p
                lRec.d = lRec.p - lRec.ref;
                lRec.dist = lRec.d.norm();
                lRec.d /= lRec.dist;

                // This shit doesn't work
                Ray3f ray;
                ray.o = lRec.ref;
                ray.d = lRec.d;

                float bigV;
                Intersection its;
                if(scene->rayIntersect(ray, its)){
                    qDebug()<<"intersect with some shit";
                    bigV = 1;
                }

                // cos theta prime prime :)
                float cosThetaI = std::max(0.0f, lRec.n.dot(-lRec.d));

                return bigV*cosThetaI / (lRec.dist * lRec.dist);
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


                // lRec.ref = its.p ??????
                LuminaireQueryRecord lRec;
                lRec.ref = its.p;
                sampleLights(scene, lRec, sampler->next2D());



                return Color3f(0.0f);
        }

        QString toString() const {
                return QString("LightIntegrator[]");
        }
};

GROUP_NAMESPACE_END

NORI_REGISTER_GROUP_CLASS(LightIntegrator, "light")
NORI_NAMESPACE_END
