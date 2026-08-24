# voxel-editor

An exact boundary-representation polyhedron engine (union / subtract /
intersect over vert-edge-face structures) with SDL voxel-editor frontends,
in C++ and Python.

## Building

Eigen is vendored in `eigen/` — the only external dependency is SDL2
(plus OpenMP for the parallel engine).

```
# macOS
brew install sdl2

# Debian / Ubuntu
sudo apt install libsdl2-dev
```

Then:

```
make
```

This builds every C++ tool:

| target | what it is | needs |
|---|---|---|
| `voxel_editor` | SDL voxel editor, first version | SDL2 |
| `voxel_editor3` | SDL voxel editor, current version | SDL2 |
| `intersect_polyhedron` | boolean kernel used as a subprocess by the Python engine | — |
| `polyhedron` | standalone C++ port of the engine | — |
| `polyhedron_parallel` | multithreaded engine | OpenMP¹ |

¹ On macOS: `brew install libomp` (the Makefile finds it automatically; the
target is skipped if OpenMP is unavailable). On Linux, GCC/Clang's built-in
`-fopenmp` is used.

`make editor` builds just the two SDL editors; `make clean` removes binaries.

## Python

The Python engine (`polyhedron.py`) shells out to `intersect_polyhedron`
for boolean operations — build that target first. `voxel_editor3.py` is the
current Python frontend; `camera.py`, `utility.py`, `view_intersects.py`
are supporting modules; `test_polyhedron.py` exercises the engine.

## Vendored dependencies

`eigen/` contains [Eigen 3.4.0](https://gitlab.com/libeigen/eigen)
(header-only, MPL2-licensed — see `eigen/COPYING.MPL2`).
