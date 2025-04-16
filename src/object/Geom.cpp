#include "Geom.h"
#include <iostream>
#include "Param.h"
#include <unordered_map>

void Geom::SortSegment2Element(Element *element)
{
    auto& vertices = element->getvertices();

    // distinguish wall or fluid vertice
    for(auto& vertice : vertices){
        auto p = vertice->getPosition();
        if(cn_PnPolyXpositive({p.x(), p.y()}) && cn_PnPolyXnegetive({p.x(), p.y()}) && cn_PnPolyYpositive({p.x(), p.y()}) && cn_PnPolyYnegetive({p.x(), p.y()})){
            vertice->setWall();
        }
    }
    // compute intersection point in cell dege

    

}

int Geom::cn_PnPolyXpositive(const Eigen::Vector2d &P)
{
    int cn = 0;
    // loop through all edges of polygon
    for(int i = 0; i < numLagrangianPoints; ++i){
        auto index1 = i;
        auto index2 = (i + 1) % numLagrangianPoints;
        auto point1 = m_points[index1];
        auto point2 = m_points[index2];
        auto V1 = Eigen::Vector2d{point1.getPosition()[0], point1.getPosition()[1]};
        auto V2 = Eigen::Vector2d{point2.getPosition()[0], point2.getPosition()[1]};

        if(((V1.y() <= P.y()) && (V2.y() > P.y())) // an upward crossing
            ||((V1.y() > P.y()) && (V2.y() <= P.y()))){ //an downward crossing
                double vt = static_cast<double>((P.y() - V1.y()) / (V2.y() - V1.y()));
                if (std::abs(V2.y() - V1.y()) < 1e-12) {
                    continue; // Avoid division by zero
                    std::cerr <<"infinity"<<std::endl;
                }
                double ixc = V1.x() + vt * (V2.x() - V1.x());
                if(P.x() < ixc){
                    cn = !cn;
                }
            }
    }

    return cn; // 0 if even (out), and 1 if odd (in)
}

int Geom::cn_PnPolyYpositive(const Eigen::Vector2d &P)
{
    int cn = 0;
    for (int i = 0; i < numLagrangianPoints; ++i) {
        int index1 = i;
        int index2 = (i + 1) % numLagrangianPoints;

        auto point1 = m_points[index1];
        auto point2 = m_points[index2];
        Eigen::Vector2d V1(point1.getPosition()[0], point1.getPosition()[1]);
        Eigen::Vector2d V2(point2.getPosition()[0], point2.getPosition()[1]);

        if (((V1.x() <= P.x()) && (V2.x() > P.x())) || ((V1.x() > P.x()) && (V2.x() <= P.x()))) {
            double vt = (P.x() - V1.x()) / (V2.x() - V1.x());
            if (std::abs(V2.x() - V1.x()) < 1e-12) {
                continue; // Avoid division by zero
                std::cerr <<"infinity"<<std::endl;
            }
            double iyc = V1.y() + vt * (V2.y() - V1.y());

            if (P.y() < iyc) {
                cn = !cn;
            }
        }
    }

    return cn;
}


int Geom::cn_PnPolyXnegetive(const Eigen::Vector2d &P)
{
    int cn = 0;
    // loop through all edges of polygon

    for(int i = 0; i < numLagrangianPoints; ++i){
        auto index1 = i;
        auto index2 = (i + 1) % numLagrangianPoints;
        auto point1 = m_points[index1];
        auto point2 = m_points[index2];
        auto V1 = Eigen::Vector2d{point1.getPosition()[0], point1.getPosition()[1]};
        auto V2 = Eigen::Vector2d{point2.getPosition()[0], point2.getPosition()[1]};

        if(((V1.y() <= P.y()) && (V2.y() > P.y())) // an upward crossing
            ||((V1.y() > P.y()) && (V2.y() <= P.y()))){ //an downward crossing
                double vt = static_cast<double>((P.y() - V1.y()) / (V2.y() - V1.y()));
                if (std::abs(V2.y() - V1.y()) < 1e-12) {
                    continue; // Avoid division by zero
                    std::cerr <<"infinity"<<std::endl;
                }
                double ixc = V1.x() + vt * (V2.x() - V1.x());
                if(P.x() > ixc){
                    cn = !cn;
                }
            }
    }

    return cn;
}

int Geom::cn_PnPolyYnegetive(const Eigen::Vector2d &P)
{   
    int cn = 0;
    for (int i = 0; i < numLagrangianPoints; ++i) {
        int index1 = i;
        int index2 = (i + 1) % numLagrangianPoints;

        auto point1 = m_points[index1];
        auto point2 = m_points[index2];
        Eigen::Vector2d V1(point1.getPosition()[0], point1.getPosition()[1]);
        Eigen::Vector2d V2(point2.getPosition()[0], point2.getPosition()[1]);

        if (((V1.x() <= P.x()) && (V2.x() > P.x())) || ((V1.x() > P.x()) && (V2.x() <= P.x()))) {
            double vt = (P.x() - V1.x()) / (V2.x() - V1.x());
            if (std::abs(V2.x() - V1.x()) < 1e-12) {
                continue; // Avoid division by zero
                std::cerr <<"infinity"<<std::endl;
            }
            double iyc = V1.y() + vt * (V2.y() - V1.y());

            if (P.y() > iyc) {
                cn = !cn;
            }
        }
    }

    return cn;
}
