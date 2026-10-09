Ok so the main data is carried by te `Mesh` class. This is used in the solver. 

It has an edges struct with various edge parameters cached in the struct. Each struct stores the properties of edge in one struct and the connectivity in another struct.
It also has a nodes struct that stores positions. This is main class for node information. 

I also have an edgeKey data structure that is used to uniquely identify edges in the mesh. This is useful for operations that require edge lookups or modifications.

`segmentDiscontinuities` is for internal cuts and  NOT boundaries.
Maybe I didn't use boundaries because it cut off the edges connecting to the boundaries


`findNearestNeighbours()` is pretty self explanatory. Uses nanoflann to find nearest neighbours of a point in the mesh.

`doesEdgeIntersectDiscontinuity()` checks if a given edge intersects with any discontinuities in the mesh. This is important for ensuring that the mesh remains valid and does not have any overlapping or intersecting edges.
Uses `segmentDiscontinuities` to check for intersections.


The `Nodes` class is position info. `nodePosFlat` is for Eigen solver. 

I need to figure out what Boundary and BoundaryEdgeInformation are

Ok so `BoundaryEdgeInformation` contains the boundary and the start and end indices of the boundaryNodes corresponding to a particular boundary id.
This might be a bit of a problem. We will need to reorganize boundayNodes indices so the start and end are sufficient to identify the boundary nodes. 
Basically all boundary nodes belonging to the same boundary should be contiguous in the boundaryNodes array.

Ok but where is the boundaryNodes array though. 

Dear past me... you idiot. 
Why are you storing boundary vec positions directly in the boundary struct. 
Why not refer by indices???
Okay okay, i don't think I store all boundary points, just ones needed to define the boundary shape. 



