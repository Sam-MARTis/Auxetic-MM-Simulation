#pragma once
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


class HyperSphere: public Domain
{
public:
    HyperSphere(const Vec2D& center, double radius)
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

//

