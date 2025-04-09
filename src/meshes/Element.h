#pragma once
#include "./Face.h"
#include <array>
#include <memory>
#include <Eigen/Dense>

class Element 
{
public:
    using Coord = Eigen::Vector2d;
    Element() = default;
    virtual ~Element() = default;

    // Accessors
    const Coord& getposition() const;
    const double& getL1() const;
    const double& getL2() const;
    const double& getvolume() const;
    const std::array<std::unique_ptr<Face>, 4>& getfaces() const;
    Face* getface(size_t index) const;
    // Note: The index should be in the range [0, 3] for a 2D element
    // Modifiers
    void setposition(const Coord& position);
    void setL1(const double& L1);
    void setL2(const double& L2);
    void setvolume(const double& volume);
    void setface(size_t index, std::unique_ptr<Face>&& face);
    void setface(size_t index, std::unique_ptr<Face>& face);

protected:
    Coord m_position; // Position of the element
    double m_L1; // Length in the first dimension
    double m_L2; // Length in the second dimension
    double m_volume;  // Volume of the element
    std::array<std::unique_ptr<Face>, 4> m_faces; // Array of faces associated with the element, order is left, right, bottom, top
    // Note: The size of the array should match the number of faces for the element type
    // For example, a 2D quadrilateral element has 4 faces
};