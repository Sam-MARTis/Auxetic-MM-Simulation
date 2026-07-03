#pragma once



struct Vec2D
{
    double x;
    double y;
    Vec2D(double x=0, double y=0){ this->x=x; this->y=y; };
    Vec2D operator+(const Vec2D& other) const
    {
        return {x + other.x, y + other.y};
    }
    Vec2D operator-(const Vec2D& other) const
    {
        return {x - other.x, y - other.y};
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

