#pragma once

#include <cmath>
#include <stdint.h>


inline double randValueMine(const double min, const double max)
{
    return min + (max-min)*(1.0*rand()/(1.0*RAND_MAX));
}



struct Vec2D
{
    double x;
    double y;
    Vec2D(double x=0, double y=0){ this->x=x; this->y=y; };
    Vec2D operator+(const Vec2D& other) const
    {
        return {x + other.x, y + other.y};
    }
    void operator+=(const Vec2D& other)
    {
        x += other.x;
        y += other.y;
    }
    void operator-=(const Vec2D& other)
    {
        x -= other.x;
        y -= other.y;
    }
    Vec2D operator*(const double multiplier) const
    {
        return {x*multiplier, y*multiplier};
    }
    void operator*=(const double multiplier)
    {
        x *= multiplier;
        y *= multiplier;
    }
    Vec2D operator/(const double divisor) const
    {
        const double multiplier = 1./divisor;
        return {x*multiplier, y*multiplier};
    }
    Vec2D operator-(const Vec2D& other) const
    {
        return {x - other.x, y - other.y};
    }
    [[nodiscard]] double norm() const
    {
        return sqrt(x*x + y*y);
    }
    [[nodiscard]] Vec2D normalized() const
    {
        const double magnitude = norm();
        const double multiplier = 1./(magnitude);
        return {x*multiplier, y*multiplier};
    }

};

inline double distanceSquared(const Vec2D& p1, const Vec2D& p2)
{
    return (p1.x - p2.x)*(p1.x - p2.x) + (p1.y - p2.y)*(p1.y - p2.y);
}

inline double crossVec2D(const Vec2D& p1, const Vec2D& p2)
{
    return p1.x*p2.y - p1.y*p2.x;
}
inline bool orientation(const Vec2D& p, const Vec2D& q, const Vec2D& r)
{
    double val = (q.y-p.y)*(r.x - q.x) - (q.x-p.x)*(r.y - q.y);
    return val <= 0;
}
inline bool segmentsIntersects2D(const Vec2D& p1, const Vec2D& p2, const Vec2D& p3, const Vec2D& p4)
{
    const bool o1 = orientation(p1, p2, p3);
    const bool o2 = orientation(p1, p2, p4);
    const bool o3 = orientation(p3, p4, p1);
    const bool o4 = orientation(p3, p4, p2);
    return o1!=o2 && o3!=o4;
}

