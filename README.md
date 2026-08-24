# voxel-editor

An exact boundary-representation polyhedron engine (union / subtract /
intersect over vert-edge-face structures) with SDL voxel-editor frontends,
in C++ and Python.

## Building

Every compile-time dependency is vendored, so all `.cpp` files compile
with no include paths at all (`c++ -std=c++17 file.cpp`). The only thing
to install is the SDL2 *library*, needed when linking the two editors:

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

| target | what it is | links against |
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

- `mini_eigen.h` — a self-contained implementation of the slice of Eigen
  this project uses (dynamic `MatrixXd`/`VectorXd` and column-pivoted
  Householder QR solves on small dense systems). Drop-in compatible with
  the Eigen spellings in the source; verified against Eigen 3.4 by fuzzing
  (200k systems, zero accept/reject divergences at engine tolerance) and
  by an end-to-end run of `intersect_polyhedron` producing geometrically
  identical output.
- `eigen/` — the full [Eigen 3.4.0](https://gitlab.com/libeigen/eigen)
  headers (MPL2, see `eigen/COPYING.MPL2`), kept available for anything
  that wants the real library; the .cpp files themselves build against
  `mini_eigen.h`.
- `SDL2/` — the [SDL2](https://libsdl.org) 2.32.8 public headers
  (zlib license, see `SDL2/LICENSE.txt`), so the editors compile without
  SDL installed. Linking still uses your system's SDL2 library. The
  GLES/EGL/opengl extension headers are omitted (not used here).
