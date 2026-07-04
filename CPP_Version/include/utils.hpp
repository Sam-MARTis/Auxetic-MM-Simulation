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
    void operator+=(const Vec2D& other)
    {
        x += other.x;
        y += other.y;
    }
    Vec2D operator*(const double multiplier) const
    {
        return {x*multiplier, y*multiplier};
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
        return std::sqrt(x*x + y*y);
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

z