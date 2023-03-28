#include "terrain.h"

#include <stdexcept>

#include "../util/helper.h"
#include <iostream>

namespace simulation {
// Factory
std::unique_ptr<Terrain> TerrainFactory::CreateTerrain(TerrainType type) {
    switch (type) {
        case simulation::TerrainType::Plane:
            return std::make_unique<PlaneTerrain>();
        case simulation::TerrainType::Bowl:
            return std::make_unique<BowlTerrain>();
        default:
            throw std::invalid_argument("TerrainFactory::CreateTerrain : invalid TerrainType");
            break;
    }
    return nullptr;
}
// Terrain

Eigen::Matrix4f Terrain::getModelMatrix() { return modelMatrix; }

// Note:
// You should update each particles' velocity (base on the equation in
// slide) and force (contact force : resist + friction) in handleCollision function

// PlaneTerrain //

PlaneTerrain::PlaneTerrain() { modelMatrix = util::translate(0.0f, position[1], 0.0f) * util::scale(60, 1, 60); }

TerrainType PlaneTerrain::getType() { return TerrainType::Plane; }

void PlaneTerrain::updatePose(const double time) {
    position[1] = hole_position[1] = std::sin(time);
    vel[1] = std::cos(time);
}

Eigen::Vector3f PlaneTerrain::getPose() { return position; }

void PlaneTerrain::handleCollision(const float delta_T, Jelly& jelly) {
    constexpr float eEPSILON = 0.01f;
    constexpr float coefResist = 0.8f;
    constexpr float coefFriction = 0.2f;
    // TODO#3-1: Handle collision when a particle collide with the plane terrain.
    //   If collision happens:
    //      1. Directly update particles' velocity
    //      2. Apply contact force to particles when needed
    // Note:
    //   1. There are `jelly.getParticleNum()` particles.
    //   2. See TODOs in `Jelly::computeInternalForce` and functions in `particle.h` if you don't know how to access
    //   data.
    // Hint:
    //   1. Review "particles.pptx¡¨ from p.14 - p.19
    //   1. Use a.norm() to get length of a.
    //   2. Use a.normalize() to normalize a inplace.
    //          a.normalized() will create a new vector.
    //   3. Use a.dot(b) to get dot product of a and b.
    for (int i = 0; i < jelly.getParticleNum(); i++) {
        Particle &particle = jelly.getParticle(i);
        Eigen::Vector3f x = particle.getPosition();
        if ((x - hole_position).norm() >= hole_radius) {
            Eigen::Vector3f v = particle.getVelocity() - vel;
            Eigen::Vector3f vn = v.dot(normal) * normal;
            Eigen::Vector3f vt = v - vn;
            Eigen::Vector3f f = particle.getForce();

            if (normal.dot(x - position) < eEPSILON) {
                if (normal.dot(v) < -eEPSILON) {
                    // collision
                    particle.setVelocity(-coefResist * vn + vt);
                }

                if (normal.dot(f) < 0) {
                    // contact forces
                    Eigen::Vector3f fc = -normal.dot(f) * normal;
                    Eigen::Vector3f ff = -coefFriction * (-normal.dot(f)) * vt;
                    particle.addForce(fc + ff);
                }
            }
        }
    }
    vel = Eigen::Vector3f::Zero();
}
// BowlTerrain //

BowlTerrain::BowlTerrain() {
    modelMatrix = util::translate(position) * util::rotateDegree(-90, 0, 0) * util::scale(radius, radius, radius);
}

TerrainType BowlTerrain::getType() { return TerrainType::Bowl; }

void BowlTerrain::updatePose(const double time) { 
    position[1] = std::sin(time) + radius / sqrt(2);
    vel[1] = std::cos(time);
}

Eigen::Vector3f BowlTerrain::getPose() { return position; }

void BowlTerrain::handleCollision(const float delta_T, Jelly& jelly) {
    constexpr float eEPSILON = 0.01f;
    constexpr float coefResist = 0.8f;
    constexpr float coefFriction = 0.2f;
    // TODO#3-2: Handle collision when a particle collide with the sphere terrain.
    //   If collision happens:
    //      1. Directly update particles' velocity
    //      2. Apply contact force to particles when needed
    // Note:
    //   1. There are `jelly.getParticleNum()` particles.
    //   2. See TODOs in `Jelly::computeInternalForce` and functions in `particle.h` if you don't know how to access
    //   data. 
    // Hint:
    //   1. Review "particles.pptx¡¨ from p.14 - p.19
    //   1. Use a.norm() to get length of a.
    //   2. Use a.normalize() to normalize a inplace.
    //          a.normalized() will create a new vector.
    //   3. Use a.dot(b) to get dot product of a and b.
    
    for (int i = 0; i < jelly.getParticleNum(); i++) {
        Particle& particle = jelly.getParticle(i);
        Eigen::Vector3f x = particle.getPosition();
        float d = (x - position).norm();
        if (x[1] <= position[1] - radius / sqrt(2) && 
            sqrt(pow(x[0] - position[0], 2) + pow(x[2] - position[2], 2)) < radius / sqrt(2)) {

            Eigen::Vector3f v = particle.getVelocity() - vel;
            Eigen::Vector3f f = particle.getForce();
            Eigen::Vector3f normal = (position - x).normalized();
            Eigen::Vector3f vn = v.dot(normal) * normal;
            Eigen::Vector3f vt = v - vn;

            if (radius - d < eEPSILON) {
                if (normal.dot(v) < -eEPSILON) {
                    // collision
                    particle.setVelocity(-coefResist * vn + vt);
                }

                if (normal.dot(f) < 0) {
                    // contact forces
                    Eigen::Vector3f fc = -normal.dot(f) * normal;
                    Eigen::Vector3f ff = -coefFriction * (-normal.dot(f)) * vt;
                    particle.addForce(fc + ff);
                }
            }
        }
    }
    vel = Eigen::Vector3f::Zero();
}
}  // namespace simulation
