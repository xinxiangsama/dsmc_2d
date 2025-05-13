#include "Segment.h"
#include "random/Random.h"
#include <iostream>
#include <memory>
#include "Param.h"
extern std::unique_ptr<Random> randomgenerator;

bool Segment::isHit(const Particle::Coord &position) const
{   
    Particle::Coord m_point = 0.5 * (m_leftpoint->getPosition() + m_rightpoint->getPosition());
    // Particle::Coord m_point = m_rightpoint->getPosition();
    // Check if the particle is within a certain distance from the wall
    double distance = (position - m_point).dot(m_normal);
    return distance < 0.0;
    // return true;
}

void Segment::Reflect(Particle* particle, const double &dt) const
{
    Particle::Coord m_point = 0.5 * (m_leftpoint->getPosition() + m_rightpoint->getPosition());
    auto velocity = particle->getvelocity();
    auto position = particle->getposition();
    position -= velocity * dt; //back to previous position
    auto distance = (position - m_point).dot(m_normal);
    auto v_normal = velocity.dot(m_normal);
    if(v_normal > 0.0){
        return;
    }
    auto t_hit = abs( - distance / v_normal);
    t_hit = std::min(t_hit, dt);
    if(t_hit > dt || t_hit < 0.0){
        return;
    }
    auto hit_position = position + velocity * t_hit;
    /*=================Specular reflection===============*/
    // auto new_velocity = velocity  + 2 * abs(v_normal) * m_normal;
    // ================= 漫反射采样（Maxwell速度分布） =================
    // 构造局部坐标系：法线 m_normal + 两个正交切向量 tang1, tang2
    Particle::Coord tang1, tang2;
    // if (std::abs(m_normal[0]) > 0.01)
    //     tang1 = Particle::Coord(0.0, 1.0, 0.0).cross(m_normal).normalized();
    // else
    //     tang1 = Particle::Coord(1.0, 0.0, 0.0).cross(m_normal).normalized();

    tang1 = Eigen::Vector3d(-m_normal.y(), m_normal.x(), m_normal.z());
    tang2 = m_normal.cross(tang1).normalized();

    // auto tan1 {Eigen::Vector3d(-m_normal.y(), m_normal.x(), m_normal.z())};
    // 热速度采样
    double rand1 = randomgenerator->getrandom01();
    double V = sqrt(-log(rand1)) * sqrt(2.0 * boltz * T_wall / mass);

    double theta = 2.0 * M_PI * randomgenerator->getrandom01();
    double v = V * cos(theta);
    double w = V * sin(theta);
    double u = sqrt(-log(randomgenerator->getrandom01())) * sqrt(2.0 * boltz * T_wall / mass);  // 法向正方向

    // 在局部坐标中拼装新速度，然后映射回全局坐标
    Particle::Coord new_velocity = u * m_normal + v * tang1 + w * tang2;
    // Particle::Coord new_velocity = u * m_normal + v * tan1;
    particle->setvelocity(new_velocity);

    // auto newErot = -log(randomgenerator->getrandom01()) * boltz * T_wall;
    // particle->setRotationalEnergy(newErot);
    // Update the position of the particle
    auto new_position = hit_position + particle->getvelocity() * abs(dt - t_hit);
    particle->setposition(new_position);


    /*==============动量统计===============*/
    // auto tangent = Eigen::Vector3d(-m_normal.y(), m_normal.x(), m_normal.z());
    // auto delta_normal = mass * Fn * std::abs((velocity - new_velocity).dot(m_normal));
    auto delta_normal = mass * Fn * (u - velocity.dot(m_normal));
    // auto delta_tangent = mass * Fn * ((velocity - new_velocity) - ((velocity - new_velocity).dot(m_normal)) * m_normal).norm();
    // auto delta_tangent = mass * Fn * std::abs((velocity - new_velocity).dot(tangent));
    auto delta_tangent = mass * Fn * (v - velocity.dot(tang1));

    auto delta_horizontal = mass * Fn * (velocity.x() - new_velocity.x()) / tau;

    normal_momentum += delta_normal;
    tangent_momentum += delta_tangent;
    horizontal_momentum += delta_horizontal;
}

const double &Segment::getlength() const
{
    return m_length;
}

const LargrangianPoint::Coord &Segment::getsloop() const
{
    return m_sloop;
}

const LargrangianPoint::Coord &Segment::getnormal() const
{
    return m_normal;
}

LargrangianPoint *Segment::getleftpoint() const
{
    return m_leftpoint;
}

LargrangianPoint *Segment::getrightpoint() const
{
    return m_rightpoint;
}

const double &Segment::getNormalMomentum()
{
    return normal_momentum;
}

const double &Segment::getTangentMomemtum()
{
    return tangent_momentum;
}

const double &Segment::getHorizontalMometum()
{
    return horizontal_momentum;
}

void Segment::setlength(const double &length)
{
    m_length = length;
}

void Segment::setsloop(const Coord &sloop)
{
    m_sloop = sloop;
}

void Segment::setnormal(const Coord &normal)
{
    m_normal = normal;
}

void Segment::setleftpoint(LargrangianPoint *leftpoint)
{
    m_leftpoint = leftpoint;
}

void Segment::setrightpoint(LargrangianPoint *rightpoint)
{
    m_rightpoint = rightpoint;
}

void Segment::clearNormalMomentum()
{
    normal_momentum = 0;
}

void Segment::clearTangentMomemtum()
{
    tangent_momentum = 0;
}

void Segment::clearHorizontalMomentum()
{
    horizontal_momentum = 0;
}
