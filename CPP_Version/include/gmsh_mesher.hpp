#pragma once

#include <string>
#include <utility>
#include <vector>

#include "constants.hpp"
#include "mesh.hpp"
#include "utils.hpp"

struct GmshMeshOptions
{
    // Target element size in the (normalised) design space.
    double targetMeshSize = DEFAULT_GMSH_MESH_SIZE;
    // 0 (default) meshes the cuts as exact slits: the mesh conforms to the cut and
    // the nodes along it are split into two coincident rows. A positive value
    // removes a slot of that width instead.
    double slotWidth = 0.0;
    // How far a slot runs past the outline, so the cut cleanly opens the notch.
    double slotExtension = 0.05;
    bool writeMeshFile = false;
    std::string meshFilePath = "gmsh_mesh.msh";
    int verbosity = 0;
};

// Meshes a closed boundary polygon with the Gmsh OpenCASCADE kernel and returns a
// spring network built from the resulting triangle edges. See
// DeveloperReference.md section 5.
class GmshMesher
{
public:
    static SpringMeshDescription build(const std::vector<Vec2D>& boundaryPolygon,
                                       const std::vector<std::pair<Vec2D, Vec2D>>& cuts,
                                       const GmshMeshOptions& options = {});
};
