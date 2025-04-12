#include "Geom.h"
#include <iostream>

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

    for (auto& segment : m_segments)
    {   
        const auto& p1 = segment.getleftpoint()->getPosition();
        const auto& p2 = segment.getrightpoint()->getPosition();

        // compute the cubic box of seement
        double seg_xmin = std::min(p1[0], p2[0]);
        double seg_xmax = std::max(p1[0], p2[0]);
        double seg_ymin = std::min(p1[1], p2[1]);
        double seg_ymax = std::max(p1[1], p2[1]);


        if (segmentIntersectsBox({p1[0], p1[1]}, {p2[0], p2[1]}, xmin, xmax, ymin, ymax)) {
            element->insertsegment(&segment);
        }
        // // 检查是否有重叠（快速剔除）
        // if (xmax < seg_xmin || xmin > seg_xmax || ymax < seg_ymin || ymin > seg_ymax)
        // continue;

        // // 粗略判断交集，可进一步精确判断相交
        // element->insertsegment(&segment);
    }
}

bool ccw(const std::array<double, 2>& A, const std::array<double, 2>& B, const std::array<double, 2>& C) {
    return (C[1]-A[1]) * (B[0]-A[0]) > (B[1]-A[1]) * (C[0]-A[0]);
}

bool segmentIntersectsSegment(const std::array<double, 2>& A, const std::array<double, 2>& B,
                               const std::array<double, 2>& C, const std::array<double, 2>& D) {
    return (ccw(A, C, D) != ccw(B, C, D)) && (ccw(A, B, C) != ccw(A, B, D));
}
bool Geom::segmentIntersectsBox(const std::array<double, 2>& p1, const std::array<double, 2>& p2,
    double xmin, double xmax, double ymin, double ymax) {
    // 包围盒快速剔除
    if (xmax < std::min(p1[0], p2[0]) || xmin > std::max(p1[0], p2[0]) ||
    ymax < std::min(p1[1], p2[1]) || ymin > std::max(p1[1], p2[1])) {
        return false;
    }

    // 定义矩形四个顶点
    std::array<double, 2> bl = {xmin, ymin};  // bottom-left
    std::array<double, 2> br = {xmax, ymin};  // bottom-right
    std::array<double, 2> tr = {xmax, ymax};  // top-right
    std::array<double, 2> tl = {xmin, ymax};  // top-left

    // 测试与矩形四边是否相交
    if (segmentIntersectsSegment(p1, p2, bl, br)) return true;
    if (segmentIntersectsSegment(p1, p2, br, tr)) return true;
    if (segmentIntersectsSegment(p1, p2, tr, tl)) return true;
    if (segmentIntersectsSegment(p1, p2, tl, bl)) return true;

    // 线段完全包含在矩形内（不与边界相交，但在内部）
    if (p1[0] >= xmin && p1[0] <= xmax && p1[1] >= ymin && p1[1] <= ymax &&
    p2[0] >= xmin && p2[0] <= xmax && p2[1] >= ymin && p2[1] <= ymax) {
        return true;
    }

    return false;
}
