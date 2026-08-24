# Builds every C++ tool in this repo. Eigen is vendored in eigen/ -- the only
# external dependency is SDL2 (and OpenMP for polyhedron_parallel).
#
#   macOS:  brew install sdl2
#   Debian: sudo apt install libsdl2-dev
#
# Then:  make          (builds everything)
#        make editor   (just the SDL editors)

CXX      ?= c++
CXXFLAGS ?= -std=c++17 -O2
EIGEN    := -I eigen

SDL_CFLAGS := $(shell sdl2-config --cflags 2>/dev/null || pkg-config --cflags sdl2 2>/dev/null)
SDL_LIBS   := $(shell sdl2-config --libs   2>/dev/null || pkg-config --libs   sdl2 2>/dev/null)

# OpenMP: GCC and Linux clang take -fopenmp directly; Apple clang needs libomp
# (brew install libomp). Auto-detected below.
UNAME := $(shell uname -s)
ifeq ($(UNAME),Darwin)
  OMP_PREFIX := $(shell brew --prefix libomp 2>/dev/null)
  ifneq ($(OMP_PREFIX),)
    OMP_CFLAGS := -Xpreprocessor -fopenmp -I$(OMP_PREFIX)/include
    OMP_LIBS   := -L$(OMP_PREFIX)/lib -lomp
  endif
else
  OMP_CFLAGS := -fopenmp
  OMP_LIBS   := -fopenmp
endif

ALL := voxel_editor voxel_editor3 intersect_polyhedron polyhedron
ifneq ($(OMP_CFLAGS),)
  ALL += polyhedron_parallel
endif

all: $(ALL)

editor: voxel_editor voxel_editor3

voxel_editor: voxel_editor.cpp
	$(CXX) $(CXXFLAGS) $(EIGEN) $(SDL_CFLAGS) -o $@ $< $(SDL_LIBS)

voxel_editor3: voxel_editor3.cpp
	$(CXX) $(CXXFLAGS) $(EIGEN) $(SDL_CFLAGS) -o $@ $< $(SDL_LIBS)

intersect_polyhedron: intersect_polyhedron.cpp
	$(CXX) $(CXXFLAGS) $(EIGEN) -o $@ $<

polyhedron: polyhedron.cpp
	$(CXX) $(CXXFLAGS) $(EIGEN) -o $@ $<

polyhedron_parallel: polyhedron_parallel.cpp
	$(CXX) $(CXXFLAGS) $(EIGEN) $(OMP_CFLAGS) -o $@ $< $(OMP_LIBS)

clean:
	rm -f voxel_editor voxel_editor3 intersect_polyhedron polyhedron polyhedron_parallel
	rm -rf *.dSYM

.PHONY: all editor clean
