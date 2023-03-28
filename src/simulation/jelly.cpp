#include "jelly.h"

#include "Eigen/Dense"
#include <iostream>

#include "../util/helper.h"
namespace simulation {
constexpr float g_cdK = 2500.0f;
constexpr float g_cdD = 50.0f;

Jelly::Jelly()
    : particleNumPerEdge(10),
      jellyLength(2.0),
      initialPosition(Eigen::Vector3f(0.0, 0.0, 0.0)),
      springCoefStruct(g_cdK),
      springCoefShear(g_cdK),
      springCoefBending(g_cdK),
      damperCoefStruct(g_cdD),
      damperCoefShear(g_cdD),
      damperCoefBending(g_cdD) {
    particleNumPerFace = particleNumPerEdge * particleNumPerEdge;
    initializeParticle();
    initializeSpring();
}

Jelly::Jelly(const Eigen::Vector3f &a_kInitPos, const float jellyLength, const int numAtEdge, const float dSpringCoef,
           const float dDamperCoef)
    : particleNumPerEdge(numAtEdge),
      jellyLength(jellyLength),
      initialPosition(a_kInitPos),
      springCoefStruct(dSpringCoef),
      springCoefShear(dSpringCoef),
      springCoefBending(dSpringCoef),
      damperCoefStruct(dDamperCoef),
      damperCoefShear(dDamperCoef),
      damperCoefBending(dDamperCoef) {
    particleNumPerFace = numAtEdge * numAtEdge;
    initializeParticle();
    initializeSpring();
}

int Jelly::getParticleNum() const { return static_cast<int>(particles.size()); }

int Jelly::getSpringNum() const { return static_cast<int>(springs.size()); }

int Jelly::getNumAtEdge() const { return particleNumPerEdge; }

unsigned int Jelly::getPointMap(const int a_ciSide, const int a_ciI, const int a_ciJ) {
    int r = -1;

    switch (a_ciSide) {
        case 1:  // [a_ciI][a_ciJ][0] bottom face
            r = particleNumPerFace * a_ciI + particleNumPerEdge * a_ciJ;
            break;
        case 6:  // [a_ciI][a_ciJ][9] top face
            r = particleNumPerFace * a_ciI + particleNumPerEdge * a_ciJ + particleNumPerEdge - 1;
            break;
        case 2:  // [a_ciI][0][a_ciJ] front face
            r = particleNumPerFace * a_ciI + a_ciJ;
            break;
        case 5:  // [a_ciI][9][a_ciJ] back face
            r = particleNumPerFace * a_ciI + particleNumPerEdge * (particleNumPerEdge - 1) + a_ciJ;
            break;
        case 3:  // [0][a_ciI][a_ciJ] left face
            r = particleNumPerEdge * a_ciI + a_ciJ;
            break;
        case 4:  // [9][a_ciI][a_ciJ] ra_ciIght face
            r = particleNumPerFace * (particleNumPerEdge - 1) + particleNumPerEdge * a_ciI + a_ciJ;
            break;
    }

    return r;
}

Particle &Jelly::getParticle(int particleIdx) { return particles[particleIdx]; }

std::vector<Particle> *Jelly::getParticlePointer() { return &particles; }

Spring &Jelly::getSpring(int springIdx) { return springs[springIdx]; }

void Jelly::setSpringCoef(const float springCoef, const Spring::SpringType springType) {
    if (springType == Spring::SpringType::STRUCT) {
        springCoefStruct = springCoef;
        updateSpringCoef(springCoef, Spring::SpringType::STRUCT);
    } else if (springType == Spring::SpringType::SHEAR) {
        springCoefShear = springCoef;
        updateSpringCoef(springCoef, Spring::SpringType::SHEAR);
    } else if (springType == Spring::SpringType::BENDING) {
        springCoefBending = springCoef;
        updateSpringCoef(springCoef, Spring::SpringType::BENDING);
    }
}

void Jelly::setDamperCoef(const float damperCoef, const Spring::SpringType springType) {
    if (springType == Spring::SpringType::STRUCT) {
        damperCoefStruct = damperCoef;
        updateDamperCoef(damperCoef, Spring::SpringType::STRUCT);
    } else if (springType == Spring::SpringType::SHEAR) {
        damperCoefShear = damperCoef;
        updateDamperCoef(damperCoef, Spring::SpringType::SHEAR);
    } else if (springType == Spring::SpringType::BENDING) {
        damperCoefBending = damperCoef;
        updateDamperCoef(damperCoef, Spring::SpringType::BENDING);
    }
}

void Jelly::resetJelly(const Eigen::Vector3f &offset, const float &rotate) {
    float dTheta = util::radians(rotate);  //  change angle from degree to
                                           //  radian

    for (unsigned int uiI = 0; uiI < particles.size(); uiI++) {
        int i = uiI / particleNumPerFace;
        int j = (uiI / particleNumPerEdge) % particleNumPerEdge;
        int k = uiI % particleNumPerEdge;
        float offset_x = (float)((i - particleNumPerEdge / 2) * jellyLength / (particleNumPerEdge - 1));
        float offset_y = (float)((j - particleNumPerEdge / 2) * jellyLength / (particleNumPerEdge - 1));
        float offset_z = (float)((k - particleNumPerEdge / 2) * jellyLength / (particleNumPerEdge - 1));

        Eigen::Vector3f RotateVec(offset_x, offset_y,
                                  offset_z);  //  vector from center of cube to the particle

        Eigen::AngleAxis<float> rotation(dTheta, Eigen::Vector3f(1.0f, 0.0f, 1.0f).normalized());

        RotateVec = rotation * RotateVec;

        particles[uiI].setPosition(initialPosition + offset + RotateVec);
        particles[uiI].setForce(Eigen::Vector3f::Zero());
        particles[uiI].setVelocity(Eigen::Vector3f::Zero());
    }
}

void Jelly::addForceField(const Eigen::Vector3f &force) {
    for (unsigned int uiI = 0; uiI < particles.size(); uiI++) {
        particles[uiI].setAcceleration(force);
    }
}

void Jelly::computeInternalForce() {
    // TODO#2-3: Compute the internal force (including spring force and damper force) for each spring.
    //   1. Read the start-particle and end-particle index from spring.
    //   2. Use `getPosition()` to get particle i's position and use `getVelocity()` get particle i's velocity.
    //   3. Call `computeSpringForce` and `computeDamperForce` to compute spring force and damper force.
    //   4. Compute net internal force and call `addForce` to apply the force onto particles.
    // Note:
    //   1. Direction of the force.
    for (int i = 0; i < springs.size(); i++) {
        int start_id = springs[i].getSpringStartID();
        int end_id = springs[i].getSpringEndID();
        Eigen::Vector3f start_pose = particles[start_id].getPosition();
        Eigen::Vector3f start_vel = particles[start_id].getVelocity();
        Eigen::Vector3f end_pose = particles[end_id].getPosition();
        Eigen::Vector3f end_vel = particles[end_id].getVelocity();

        Eigen::Vector3f spring_force = 
            computeSpringForce(start_pose, end_pose, springs[i].getSpringCoef(), springs[i].getSpringRestLength());
        Eigen::Vector3f damper_force = 
            computeDamperForce(start_pose, end_pose, start_vel, end_vel, springs[i].getDamperCoef());

        particles[start_id].addForce(spring_force);
        particles[start_id].addForce(damper_force);
        particles[end_id].addForce(-spring_force);
        particles[end_id].addForce(-damper_force);
    }

}

Eigen::Vector3f Jelly::computeSpringForce(const Eigen::Vector3f &positionA, const Eigen::Vector3f &positionB,
                                         const float springCoef, const float restLength) {
    // TODO#2-1: Compute spring force given the two positions of the spring.
    //   1. Review "particles.pptx¡¨ from p.9 - p.13
    //   2. The sample below just set spring force to zero
    Eigen::Vector3f delta_x = positionA - positionB;
    Eigen::Vector3f l = delta_x.normalized();
    Eigen::Vector3f spring_force = -springCoef * (delta_x.norm() - restLength) * l;
    return spring_force;
}

Eigen::Vector3f Jelly::computeDamperForce(const Eigen::Vector3f &positionA, const Eigen::Vector3f &positionB,
                                         const Eigen::Vector3f &velocityA, const Eigen::Vector3f &velocityB,
                                         const float damperCoef) {
    // TODO#2-2: Compute damper force given the two positions and the two velocities of the spring.
    //   1. Review "particles.pptx¡¨ from p.9 - p.13
    //   2. The sample below just set damper force to zero
    Eigen::Vector3f delta_x = positionA - positionB;
    Eigen::Vector3f l = delta_x.normalized();
    float delta_v = (velocityA - velocityB).dot(l);
    Eigen::Vector3f damper_force = -damperCoef * delta_v * l;
    return damper_force;
}

void Jelly::initializeParticle() {
    for (int i = 0; i < particleNumPerEdge; i++) {
        for (int j = 0; j < particleNumPerEdge; j++) {
            for (int k = 0; k < particleNumPerEdge; k++) {
                Particle Particle;
                float offset_x = (float)((i - particleNumPerEdge / 2) * jellyLength / (particleNumPerEdge - 1));
                float offset_y = (float)((j - particleNumPerEdge / 2) * jellyLength / (particleNumPerEdge - 1));
                float offset_z = (float)((k - particleNumPerEdge / 2) * jellyLength / (particleNumPerEdge - 1));
                Particle.setPosition(Eigen::Vector3f(initialPosition(0) + offset_x, initialPosition(1) + offset_y,
                                                     initialPosition(2) + offset_z));
                particles.push_back(Particle);

            }
        }
    }
}

void Jelly::connectParticles(const int iParticleID, const int iNeighborID, const Spring::SpringType springType) {
    Eigen::Vector3f SpringStartPos = particles[iParticleID].getPosition();
    Eigen::Vector3f SpringEndPos = particles[iNeighborID].getPosition();
    Eigen::Vector3f Length = SpringStartPos - SpringEndPos;
    float absLength = sqrt(Length[0] * Length[0] + Length[1] * Length[1] + Length[2] * Length[2]);
    float springCoef = 0;
    float damperCoef = 0;
    if (springType == Spring::SpringType::STRUCT) {
        springCoef = springCoefStruct;
        damperCoef = damperCoefStruct;
    }
    else if (springType == Spring::SpringType::BENDING) {
        springCoef = springCoefBending;
        damperCoef = damperCoefBending;
    }
    else if (springType == Spring::SpringType::SHEAR) {
        springCoef = springCoefShear;
        damperCoef = damperCoefShear;
    }
    springs.push_back(Spring(iParticleID, iNeighborID, absLength, springCoef, damperCoef, springType));
}

void Jelly::initializeSpring() {

    // TODO#1: Connect particles with springs.  
    //   1. Consider the type of springs and compute indices of particles which the spring connect to.
    //   2. Compute rest spring length using particle positions.
    //   2. Iterate the particles. Push spring objects into `springs` vector
    // Note:
    //   1. The particles index can be computed in a similar way as below:
    //   ===============================================
    //   0 1 2 3 ... particlesPerEdge
    //   particlesPerEdge + 1 ....
    //   ... ... particlesPerEdge * particlesPerEdge - 1
    //   ===============================================
    // Here is a simple example which connects the structrual springs along z-axis.

    int n_struct = 0;
    int n_bending = 0;
    int n_shear = 0;

    // struct
    for (int i = 0; i < particleNumPerEdge; i++) {
        for (int j = 0; j < particleNumPerEdge; j++) {
            for (int k = 0; k < particleNumPerEdge; k++) {
                int iParticleID = i * particleNumPerFace + j * particleNumPerEdge + k;
                int iNeighborID;
                if (k < particleNumPerEdge - 1) {
                    iNeighborID = iParticleID + 1;
                    connectParticles(iParticleID, iNeighborID, Spring::SpringType::STRUCT);
                    n_struct++;
                }
                if (j < particleNumPerEdge - 1) {
                    iNeighborID = iParticleID + particleNumPerEdge;
                    connectParticles(iParticleID, iNeighborID, Spring::SpringType::STRUCT);
                    n_struct++;
                }
                if (i < particleNumPerEdge - 1) {
                    iNeighborID = iParticleID + particleNumPerFace;
                    connectParticles(iParticleID, iNeighborID, Spring::SpringType::STRUCT);
                    n_struct++;
                }
            }
        }
    }

    // bending
    for (int i = 0; i < particleNumPerEdge; i++) {
        for (int j = 0; j < particleNumPerEdge; j++) {
            for (int k = 0; k < particleNumPerEdge; k++) {
                int iParticleID = i * particleNumPerFace + j * particleNumPerEdge + k;
                int iNeighborID;
                if (k < particleNumPerEdge - 2) {
                    iNeighborID = iParticleID + 2;
                    connectParticles(iParticleID, iNeighborID, Spring::SpringType::BENDING);
                    n_bending++;
                }
                if (j < particleNumPerEdge - 2) {
                    iNeighborID = iParticleID + 2 * particleNumPerEdge;
                    connectParticles(iParticleID, iNeighborID, Spring::SpringType::BENDING);
                    n_bending++;
                }
                if (i < particleNumPerEdge - 2) {
                    iNeighborID = iParticleID + 2 * particleNumPerFace;
                    connectParticles(iParticleID, iNeighborID, Spring::SpringType::BENDING);
                    n_bending++;
                }
            }
        }
    }

    // shear
    for (int i = 0; i < particleNumPerEdge; i++) {
        for (int j = 0; j < particleNumPerEdge; j++) {
            for (int k = 0; k < particleNumPerEdge; k++) {
                int iParticleID = i * particleNumPerFace + j * particleNumPerEdge + k;
                int iNeighborID;
                if (j < particleNumPerEdge - 1 && k < particleNumPerEdge - 1) {
                    iNeighborID = iParticleID + particleNumPerEdge + 1;
                    connectParticles(iParticleID, iNeighborID, Spring::SpringType::SHEAR);
                    n_shear++;
                }
                if (j < particleNumPerEdge - 1 && k > 0) {
                    iNeighborID = iParticleID + particleNumPerEdge - 1;
                    connectParticles(iParticleID, iNeighborID, Spring::SpringType::SHEAR);
                    n_shear++;
                }
                if (i < particleNumPerEdge - 1 && k < particleNumPerEdge - 1) {
                    iNeighborID = iParticleID + particleNumPerFace + 1;
                    connectParticles(iParticleID, iNeighborID, Spring::SpringType::SHEAR);
                    n_shear++;
                }
                if (i < particleNumPerEdge - 1 && k > 0) {
                    iNeighborID = iParticleID + particleNumPerFace - 1;
                    connectParticles(iParticleID, iNeighborID, Spring::SpringType::SHEAR);
                    n_shear++;
                }
                if (i < particleNumPerEdge - 1 && j < particleNumPerEdge - 1) {
                    iNeighborID = iParticleID + particleNumPerFace + particleNumPerEdge;
                    connectParticles(iParticleID, iNeighborID, Spring::SpringType::SHEAR);
                    n_shear++;
                }
                if (i < particleNumPerEdge - 1 && j > 0) {
                    iNeighborID = iParticleID + particleNumPerFace - particleNumPerEdge;
                    connectParticles(iParticleID, iNeighborID, Spring::SpringType::SHEAR);
                    n_shear++;
                }
                if (i < particleNumPerEdge - 1 && j < particleNumPerEdge - 1 && k < particleNumPerEdge - 1) {
                    iNeighborID = iParticleID + particleNumPerFace + particleNumPerEdge + 1;
                    connectParticles(iParticleID, iNeighborID, Spring::SpringType::SHEAR);
                    n_shear++;
                }
                if (i < particleNumPerEdge - 1 && j < particleNumPerEdge - 1 && k > 0) {
                    iNeighborID = iParticleID + particleNumPerFace + particleNumPerEdge - 1;
                    connectParticles(iParticleID, iNeighborID, Spring::SpringType::SHEAR);
                    n_shear++;
                }
                if (i > 0 && j < particleNumPerEdge - 1 && k < particleNumPerEdge - 1) {
                    iNeighborID = iParticleID - particleNumPerFace + particleNumPerEdge + 1;
                    connectParticles(iParticleID, iNeighborID, Spring::SpringType::SHEAR);
                    n_shear++;
                }
                if (i > 0 && j < particleNumPerEdge - 1 && k > 0) {
                    iNeighborID = iParticleID - particleNumPerFace + particleNumPerEdge - 1;
                    connectParticles(iParticleID, iNeighborID, Spring::SpringType::SHEAR);
                    n_shear++;
                }
            }
        }
    }

    std::cout << "[Number of springs]" << std::endl;
    std::cout << "struct: " << n_struct << std::endl;
    std::cout << "bending: " << n_bending << std::endl;
    std::cout << "shear: " << n_shear << std::endl;
}

void Jelly::updateSpringCoef(const float a_cdSpringCoef, const Spring::SpringType a_cSpringType) {
    for (unsigned int uiI = 0; uiI < springs.size(); uiI++) {
        if (springs[uiI].getType() == a_cSpringType) {
            springs[uiI].setSpringCoef(a_cdSpringCoef);
        }
    }
}

void Jelly::updateDamperCoef(const float a_cdDamperCoef, const Spring::SpringType a_cSpringType) {
    for (unsigned int uiI = 0; uiI < springs.size(); uiI++) {
        if (springs[uiI].getType() == a_cSpringType) {
            springs[uiI].setDamperCoef(a_cdDamperCoef);
        }
    }
}
}  //  namespace simulation
