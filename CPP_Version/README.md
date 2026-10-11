# CPP_Version

Spring-network solver for bistable/auxetic metamaterials: static energy
minimisation over a tagged boundary.

```bash
cd build
cmake --build . -j$(nproc)
./runThis                  # GUI demo (needs a display)
./boundary_tests           # headless checks
```

Local notes (gitignored): `AgentsReference.md` for the API and commands,
`DeveloperReference.md` for design rationale and traps.

Needs SFML 3, Eigen3, OpenMP, UMFPACK and nanoflann. Geometry comes from
`resources/DXFs` via `src/subroutines/dxfParser.py`.
