#include "Geom.h"
#include <iostream>
#include <unordered_map>

void Geom::SortSegment2Element(Element *element)
{
    const auto& elemPos = element->getposition();
    const double Lx = element->getL1();
    const double Ly = element->getL2();

    // compute the box boundary of element
    double xmin = elemPos[0] - 0.5 * Lx;
    double xmax = elemPos[0] + 0.5 * Lx;
    double ymin = elemPos[1] - 0.5 * Ly;
    double ymax = elemPos[1] + 0.5 * Ly;

    // for 2d problem
    Eigen::Vector2d point1{xmin, ymin};
    Eigen::Vector2d point2{xmin, ymax};
    Eigen::Vector2d point3{xmax, ymin};
    Eigen::Vector2d point4{xmax, ymax};
    auto ElementV =  std::array<Eigen::Vector2d, 4>{point1, point2, point3, point4};

    int num {};
    for(auto& P : ElementV){
        if(cn_PnPolyX(P).get() && cn_PnPolyY(P).get()){
            ++num;
            std::cout << "intersection positon is"<<cn_PnPolyX(P)->transpose()<<std::endl;
        }
    }

    if(num != 0 && num != 4){
        std::cout << "num is" << num <<std::endl;
    }
}

std::unique_ptr<Eigen::Vector2d> Geom::cn_PnPolyX(const Eigen::Vector2d &P)
{
    std::vector<std::unique_ptr<Eigen::Vector2d>> intersectionP;

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
                double ixc = V1.x() + vt * (V2.x() - V1.x());
                if(P.x() < ixc){
                    intersectionP.emplace_back(std::make_unique<Eigen::Vector2d>(ixc, P.y()));
                }
            }
    }
    std::unique_ptr<Eigen::Vector2d> res;

    if(intersectionP.size() == 1){
        res = std::move(intersectionP[0]);
    }
    return res;
}

std::unique_ptr<Eigen::Vector2d> Geom::cn_PnPolyY(const Eigen::Vector2d &P)
{
    std::vector<std::unique_ptr<Eigen::Vector2d>> intersectionP;

    for (int i = 0; i < numLagrangianPoints; ++i) {
        int index1 = i;
        int index2 = (i + 1) % numLagrangianPoints;

        auto point1 = m_points[index1];
        auto point2 = m_points[index2];
        Eigen::Vector2d V1(point1.getPosition()[0], point1.getPosition()[1]);
        Eigen::Vector2d V2(point2.getPosition()[0], point2.getPosition()[1]);

        if (((V1.x() <= P.x()) && (V2.x() > P.x())) || ((V1.x() > P.x()) && (V2.x() <= P.x()))) {
            double vt = (P.x() - V1.x()) / (V2.x() - V1.x());
            double iyc = V1.y() + vt * (V2.y() - V1.y());

            if (P.y() < iyc) {
                intersectionP.emplace_back(std::make_unique<Eigen::Vector2d>(P.x(), iyc));
            }
        }
    }

    std::unique_ptr<Eigen::Vector2d> res;

    if(intersectionP.size() == 1){
        res = std::move(intersectionP[0]);
    }
    return res;
}
