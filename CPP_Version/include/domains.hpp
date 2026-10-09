#pragma once
#include <cstddef>
#include <iostream>
#include <utility>
#include <vector>

#include "utils.hpp"
#include "vec.hpp"




class Domain
{
public:
    virtual ~Domain() = default;
    [[nodiscard]] virtual bool isInside(const Vec2D& p) const =0;
};


class Cuboid:  public Domain
{
public:
    Cuboid(const Vec2D& lower, const Vec2D& upper)
        : lowerBounds(lower), upperBounds(upper) {}

    [[nodiscard]] bool isInside(const Vec2D& p) const override
    {
        for (uint i = 0; i < 2; ++i)
        {
            if (p[i] < lowerBounds[i] || p[i] > upperBounds[i])
            {
                return false;
            }
        }
        return true;
    }

private:
    Vec2D lowerBounds;
    Vec2D upperBounds;
};


class Circle: public Domain
{
public:
    Circle(const Vec2D& center, double radius)
        : center(center), radius(radius) {}

    [[nodiscard]] bool isInside(const Vec2D& p) const override
    {
        double distanceSquared = 0;
        for (uint i = 0; i < 2; ++i)
        {
            double diff = p[i] - center[i];
            distanceSquared += diff * diff;
        }
        return distanceSquared <= radius * radius;
    }

private:
    Vec2D center;
    double radius;
};





class IntersectionDomain: public Domain
{
    const std::vector<std::pair<Domain&, bool>>& regionsAndParities;
public:
    IntersectionDomain(std::vector<std::pair<Domain&, bool>>& _regionsAndParities): regionsAndParities(_regionsAndParities) {}
    [[nodiscard]] bool isInside(const Vec2D& p) const override
    {
        for (const auto& [region, parity] : regionsAndParities)
        {
            if (region.isInside(p) != parity)
            {
                return false;
            }
        }
        return true;
    }
};


class PointHull: public Domain
{
public:
    PointHull(std::vector<Vec2D> hullPoints)
        : points(std::move(hullPoints)) {}

    [[nodiscard]] bool isInside(const Vec2D& testPoint) const override
    {
        const std::size_t numPoints = points.size();
        if (numPoints < 3) return false;

        const Vec2D rayDirection{1.0, 0.0};
        unsigned int numIntersections = 0;
        for (std::size_t i = 0; i < numPoints; i++)
        {
            const Vec2D& boundaryPoint = points[i];
            const Vec2D segmentDirection = points[(i + 1) % numPoints] - points[i];
            const double denom = VecMath::cross2D(rayDirection, segmentDirection);
            if (std::abs(denom) <= EPSILON) continue;
            const double t = VecMath::cross2D(boundaryPoint - testPoint, segmentDirection) / denom;
            if (t < 0.0) continue;
            const double u = VecMath::cross2D(boundaryPoint - testPoint, rayDirection) / denom;
            if (u >= 0.0 && u < 1.0) numIntersections++;
        }
        return numIntersections % 2 != 0;
    }

    [[nodiscard]] std::size_t numPoints() const { return points.size(); }
    [[nodiscard]] const std::vector<Vec2D>& vertices() const { return points; }

    void printPoints() const
    {
        for (const Vec2D& point : points)
        {
            std::cout << point[0] << " " << point[1] << std::endl;
        }
    }

private:
    std::vector<Vec2D> points;
};
