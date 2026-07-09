
#include <iostream>
#include <nanoflann.hpp>
#include <vector>

struct Vec2D
{
    double x, y;
    Vec2D(double x=0.0, double y=0.0) : x(x), y(y) {}
};
struct PointCloud
{
    std::vector<Vec2D> pointsPos;
    size_t kdtree_get_point_count() const { return pointsPos.size(); };
    double kdtree_get_pt(const size_t idx, int dim) const { return dim==0?pointsPos[idx].x:pointsPos[idx].y; };
    template <class BBOX>
    bool kdtree_get_bbox(BBOX&) const { return false;}
};
using KDTree = nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<double, PointCloud>, PointCloud, 2>;
int main()
{
    std::vector<Vec2D> points;
    for (int i=0; i< 100; i++)
    {
        points.emplace_back(i*i*0.01,i);
    }
    PointCloud cloud{.pointsPos =  points};
    KDTree tree(2, cloud, nanoflann::KDTreeSingleIndexAdaptorParams(10));
    tree.buildIndex();
    

    Vec2D query(12.0, 35.0);

    const size_t K = 5;
    std::vector<size_t> indices(K);
    std::vector<double> distances(K);
    nanoflann::KNNResultSet<double> result(K);
    result.init(indices.data(), distances.data());
    tree.findNeighbors(result, &query.x);
    for (int i=0; i< K; i++)
    {
        const Vec2D& near = points[indices[i]];
        std::cout << near.x << " " << near.y << std::endl;
    }

    std::cout << "Hello world!" << std::endl;
}