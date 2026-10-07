#include <iostream>
#include <vector>
#include <cmath>

#include "loader.hpp"

void check(const std::string& name, const Boundary& b, const Vec2D& p, bool expected)
{
    bool res = b.isInsideBoundary(p);
    std::cout << name << ": point(" << p[0] << "," << p[1] << ") => " << (res?"INSIDE":"OUTSIDE")
              << " [expected: " << (expected?"INSIDE":"OUTSIDE") << "]"
              << (res==expected?" OK":" FAIL") << std::endl;
}

Boundary makeSquare(double x0, double y0, double size)
{
    std::vector<Vec2D> pts = {
        {x0, y0}, {x0+size, y0}, {x0+size, y0+size}, {x0, y0+size}
    };
    return Boundary(pts);
}

Boundary makeTriangle()
{
    std::vector<Vec2D> pts = { {0.0,0.0}, {1.0,0.0}, {0.0,1.0} };
    return Boundary(pts);
}

Boundary makeRegularHexagon(double cx, double cy, double r)
{
    std::vector<Vec2D> pts;
    for (int i=0;i<6;i++) {
        double ang = M_PI/3.0 * i;
        pts.push_back({cx + r * cos(ang), cy + r * sin(ang)});
    }
    return Boundary(pts);
}

int main()
{
    std::cout<<"Running Boundary::isInsideBoundary tests\n";

    // Square tests
    Boundary square = makeSquare(0.0,0.0,1.0);
    check("Square", square, {0.5,0.5}, true);
    check("Square", square, {1.5,0.5}, false);
    check("Square", square, {0.0,0.0}, true); // vertex
    check("Square", square, {1.0,0.5}, true); // edge

    // Triangle tests
    Boundary tri = makeTriangle();
    check("Triangle", tri, {0.2,0.2}, true);
    check("Triangle", tri, {0.8,0.8}, false);
    check("Triangle", tri, {0.0,0.0}, true); // vertex

    // Hexagon tests
    Boundary hex = makeRegularHexagon(0.0, 0.0, 1.0);
    check("Hexagon", hex, {0.0,0.0}, true);
    check("Hexagon", hex, {2.0,0.0}, false);
    check("Hexagon", hex, {0.9, 0.0}, true);

    return 0;
}
