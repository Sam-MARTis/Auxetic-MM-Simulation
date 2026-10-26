# CPP_Version

Spring-network solver for bistable/auxetic metamaterials: static energy
minimisation over a tagged boundary.

```bash
cd build
cmake --build . -j$(nproc)
./runThis                  # GUI demo (needs a display)
./boundary_tests           # headless checks
./gmsh_mesh_tests          # gmsh meshing checks (needs the gmsh SDK)
```

Local notes (gitignored): `AgentsReference.md` for the API and commands,
`DeveloperReference.md` for design rationale and traps.

Geometry comes from `resources/DXFs`: the DXF is parsed with ezdxf
(`src/subroutines/dxfParser.py`), meshed through the gmsh C++ API, and the mesh
read back as a spring network. Without gmsh the legacy node-sampling mesher is
used instead.

Needs SFML 3, Eigen3, OpenMP, UMFPACK, nanoflann and (optionally) gmsh.
