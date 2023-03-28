#include "integrator.h"
#include <iostream>

#include <vector>
namespace simulation {
// Factory
std::unique_ptr<Integrator> IntegratorFactory::CreateIntegrator(IntegratorType type) {
    switch (type) {
        case simulation::IntegratorType::ExplicitEuler:
            return std::make_unique<ExplicitEulerIntegrator>();
        case simulation::IntegratorType::ImplicitEuler:
            return std::make_unique<ImplicitEulerIntegrator>();
        case simulation::IntegratorType::MidpointEuler:
            return std::make_unique<MidpointEulerIntegrator>();
        case simulation::IntegratorType::RungeKuttaFourth:
            return std::make_unique<RungeKuttaFourthIntegrator>();
        default:
            throw std::invalid_argument("TerrainFactory::CreateTerrain : invalid TerrainType");
            break;
    }
    return nullptr;
}

//
// ExplicitEulerIntegrator
//

IntegratorType ExplicitEulerIntegrator::getType() { return IntegratorType::ExplicitEuler; }

void ExplicitEulerIntegrator::integrate(MassSpringSystem& particleSystem) {

    // TODO#4-1: Integrate position and velocity
    //   1. Integrate position using current velocity.
    //   2. Integrate velocity using current acceleration.
    //   3. Clear force
    // Note:
    //   1. You should do this first. Then you can check whether your collision is correct or not.
    //   2. See functions in `particle.h` if you don't know how to access data.
    //   3. Review ¡§ODE_basics.pptx¡¨ from p.15 - p.16
    for (int i = 0; i < particleSystem.getJellyCount(); i++) {
        Jelly* jelly = particleSystem.getJellyPointer(i);
        for (int j = 0; j < jelly->getParticleNum(); j++) {
            Particle& p = jelly->getParticle(j);
            p.setPosition(p.getPosition() + particleSystem.deltaTime * p.getVelocity());
            p.setVelocity(p.getVelocity() + particleSystem.deltaTime * p.getAcceleration());
            p.setForce(Eigen::Vector3f::Zero());
        }
    }
}

//
// ImplicitEulerIntegrator
//

IntegratorType ImplicitEulerIntegrator::getType() { return IntegratorType::ImplicitEuler; }

void ImplicitEulerIntegrator::integrate(MassSpringSystem& particleSystem) {
    // TODO#4-2: Integrate position and velocity
    //   1. Backup original particles' data.
    //   2. Integrate position and velocity using explicit euler to get Xn+1.
    //   3. Compute refined Xn+1 and Vn+1 using (1.) and (2.).
    // Note:
    //   1. Use `MassSpringSystem::computeJellyForce` with modified position and velocity to get Xn+1
    //   2. Review ¡§ODE_basics.pptx¡¨ from p.18 - p.19

    std::vector<Eigen::Vector3f> pose;
    std::vector<Eigen::Vector3f> vel;
    for (int i = 0; i < particleSystem.getJellyCount(); i++) {
        Jelly* jelly = particleSystem.getJellyPointer(i);
        for (int j = 0; j < jelly->getParticleNum(); j++) {
            Particle& p = jelly->getParticle(j);
            pose.push_back(p.getPosition());
            vel.push_back(p.getVelocity());
            // p.setPosition(p.getPosition() + particleSystem.deltaTime * p.getVelocity());
            // p.setVelocity(p.getVelocity() + particleSystem.deltaTime * p.getAcceleration());
            p.setForce(Eigen::Vector3f::Zero());
        }
    }
    particleSystem.computeAllForce();
    int idx = 0;
    for (int i = 0; i < particleSystem.getJellyCount(); i++) {
        Jelly* jelly = particleSystem.getJellyPointer(i);
        for (int j = 0; j < jelly->getParticleNum(); j++) {
            Particle& p = jelly->getParticle(j);
            p.setPosition(pose[idx] + particleSystem.deltaTime * p.getVelocity());
            p.setVelocity(vel[idx] + particleSystem.deltaTime * p.getAcceleration());
            p.setForce(Eigen::Vector3f::Zero());
            idx++;
        }
    }
}

//
// MidpointEulerIntegrator
//

IntegratorType MidpointEulerIntegrator::getType() { return IntegratorType::MidpointEuler; }

void MidpointEulerIntegrator::integrate(MassSpringSystem& particleSystem) {
    // TODO#4-3: Integrate position and velocity
    //   1. Backup original particles' data.
    //   2. Integrate position and velocity using explicit euler to get Xn+1.
    //   3. Compute refined Xn+1 using (1.) and (2.).
    // Note:
    //   1. Use `MassSpringSystem::computeJellyForce` with modified position and velocity to get Xn+1.
    //   2. Review ¡§ODE_basics.pptx¡¨ from p .18 - p .20
    std::vector<Eigen::Vector3f> pose;
    std::vector<Eigen::Vector3f> vel;
    for (int i = 0; i < particleSystem.getJellyCount(); i++) {
        Jelly* jelly = particleSystem.getJellyPointer(i);
        for (int j = 0; j < jelly->getParticleNum(); j++) {
            Particle& p = jelly->getParticle(j);
            pose.push_back(p.getPosition());
            vel.push_back(p.getVelocity());
            p.setPosition(p.getPosition() + 0.5 * particleSystem.deltaTime * p.getVelocity());
            p.setVelocity(p.getVelocity() + 0.5 * particleSystem.deltaTime * p.getAcceleration());
            p.setForce(Eigen::Vector3f::Zero());
        }
    }
    particleSystem.computeAllForce();
    int idx = 0;
    for (int i = 0; i < particleSystem.getJellyCount(); i++) {
        Jelly* jelly = particleSystem.getJellyPointer(i);
        for (int j = 0; j < jelly->getParticleNum(); j++) {
            Particle& p = jelly->getParticle(j);
            p.setPosition(pose[idx] + particleSystem.deltaTime * p.getVelocity());
            p.setVelocity(vel[idx] + particleSystem.deltaTime * p.getAcceleration());
            p.setForce(Eigen::Vector3f::Zero());
            idx++;
        }
    }
}

//
// RungeKuttaFourthIntegrator
//

IntegratorType RungeKuttaFourthIntegrator::getType() { return IntegratorType::RungeKuttaFourth; }

void RungeKuttaFourthIntegrator::integrate(MassSpringSystem& particleSystem) {
    
    // TODO#4-4: Integrate velocity and acceleration
    //   1. Backup original particles' data.
    //   2. Compute k1, k2, k3, k4
    //   3. Compute refined Xn+1 using (1.) and (2.).
    // Note:
    //   1. Use `MassSpringSystem::computeJellyForce` with modified position and velocity to get Xn+1.
    //   2. StateStep struct is just a hint, you can use whatever you want.
    //   3. Review ¡§ODE_basics.pptx¡¨ from p.21
    
    std::vector<Eigen::Vector3f> pose;
    std::vector<Eigen::Vector3f> vel;

    struct StateStep {
        Eigen::Vector3f deltaVel;
        Eigen::Vector3f deltaPos;
    };
    std::vector<StateStep> k1, k2, k3, k4;

    for (int i = 0; i < particleSystem.getJellyCount(); i++) {
        Jelly* jelly = particleSystem.getJellyPointer(i);
        for (int j = 0; j < jelly->getParticleNum(); j++) {
            Particle& p = jelly->getParticle(j);
            pose.push_back(p.getPosition());
            vel.push_back(p.getVelocity());

            StateStep statestep;
            statestep.deltaPos = particleSystem.deltaTime * p.getVelocity();
            statestep.deltaVel = particleSystem.deltaTime * p.getAcceleration();
            k1.push_back(statestep);

            p.setPosition(p.getPosition() + 0.5 * statestep.deltaPos);
            p.setVelocity(p.getVelocity() + 0.5 * statestep.deltaVel);

            p.setForce(Eigen::Vector3f::Zero());
        }
    }

    particleSystem.computeAllForce();
    int idx = 0;
    for (int i = 0; i < particleSystem.getJellyCount(); i++) {
        Jelly* jelly = particleSystem.getJellyPointer(i);
        for (int j = 0; j < jelly->getParticleNum(); j++) {
            Particle& p = jelly->getParticle(j);

            StateStep statestep;
            statestep.deltaPos = particleSystem.deltaTime * p.getVelocity();
            statestep.deltaVel = particleSystem.deltaTime * p.getAcceleration();
            k2.push_back(statestep);

            p.setPosition(pose[idx] + 0.5 * statestep.deltaPos);
            p.setVelocity(vel[idx] + 0.5 * statestep.deltaVel);

            p.setForce(Eigen::Vector3f::Zero());
            idx++;
        }
    }

    particleSystem.computeAllForce();
    idx = 0;
    for (int i = 0; i < particleSystem.getJellyCount(); i++) {
        Jelly* jelly = particleSystem.getJellyPointer(i);
        for (int j = 0; j < jelly->getParticleNum(); j++) {
            Particle& p = jelly->getParticle(j);

            StateStep statestep;
            statestep.deltaPos = particleSystem.deltaTime * p.getVelocity();
            statestep.deltaVel = particleSystem.deltaTime * p.getAcceleration();
            k3.push_back(statestep);

            p.setPosition(pose[idx] + statestep.deltaPos);
            p.setVelocity(vel[idx] + statestep.deltaVel);

            p.setForce(Eigen::Vector3f::Zero());
            idx++;
        }
    }

    particleSystem.computeAllForce();
    idx = 0;
    for (int i = 0; i < particleSystem.getJellyCount(); i++) {
        Jelly* jelly = particleSystem.getJellyPointer(i);
        for (int j = 0; j < jelly->getParticleNum(); j++) {
            Particle& p = jelly->getParticle(j);

            StateStep statestep;
            statestep.deltaPos = particleSystem.deltaTime * p.getVelocity();
            statestep.deltaVel = particleSystem.deltaTime * p.getAcceleration();
            k4.push_back(statestep);

            p.setPosition(pose[idx] + (k1[idx].deltaPos + 2 * k2[idx].deltaPos + 2 * k3[idx].deltaPos + k4[idx].deltaPos) / 6);
            p.setVelocity(vel[idx] + (k1[idx].deltaVel + 2 * k2[idx].deltaVel + 2 * k3[idx].deltaVel + k4[idx].deltaVel) / 6);

            p.setForce(Eigen::Vector3f::Zero());
            idx++;
        }
    }
}
}  // namespace simulation
