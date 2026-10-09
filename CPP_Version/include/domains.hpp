#pragma once
#include "vec.hpp"



template <typename T, uint N>
class Domain
{
public:
    virtual ~Domain() = default;
    virtual bool isInside(const Vec<T, N>& p) const =0;
};

template <typename T, uint N>
class HyperCuboid:  public Domain< T, N >
{
public:
    HyperCuboid(const Vec<T, N>& lower, const Vec<T, N>& upper)
        : lowerBounds(lower), upperBounds(upper) {}

    bool isInside(const Vec<T, N>& p) const override
    {
        for (uint i = 0; i < N; ++i)
        {
            if (p[i] < lowerBounds[i] || p[i] > upperBounds[i])
            {
                return false;
            }
        }
        return true;
    }

private:
    Vec<T, N> lowerBounds;
    Vec<T, N> upperBounds;
};

template<typename T, uint N>
class HyperSphere: public Domain< T, N >
{
public:
    HyperSphere(const Vec<T, N>& center, T radius)
        : center(center), radius(radius) {}

    bool isInside(const Vec<T, N>& p) const override
    {
        T distanceSquared = 0;
        for (uint i = 0; i < N; ++i)
        {
            T diff = p[i] - center[i];
            distanceSquared += diff * diff;
        }
        return distanceSquared <= radius * radius;
    }

private:
    Vec<T, N> center;
    T radius;
};





template <typename T, uint N>
class IntersectionDomain: public Domain< T, N >
{
    const std::vector<std::pair<Domain<T, N>&, bool>>& regionsAndParities;
public:
    IntersectionDomain(std::vector<std::pair<Domain<T, N>&, bool>>& _regionsAndParities): regionsAndParities(_regionsAndParities) {}
    bool isInside(const Vec<T, N>& p) const override
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

