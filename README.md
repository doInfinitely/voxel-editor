To compile `voxel_editor.cpp` use
```
g++ -std=c++14 -I eigen/ -o voxel_editor voxel_editor.cpp `sdl2-config --cflags --libs`
```
or
```
g++ -std=c++14 -I eigen/ -o voxel_editor voxel_editor.cpp `pkg-config --cflags --libs sdl2`
```

Eigen can be obtained from https://gitlab.com/libeigen/eigen

```
g++ -std=c++14 -I eigen/ -o intersect_polyhedron intersect_polyhedron.cpp `pkg-config --cflags --libs sdl2`
```

```
g++ -std=c++14 -I eigen/ -I or-tools/include/ -o polyhedron polyhedron.cpp  -L ~/Code/voxel-editor/or-tools/lib -lortools `pkg-config --cflags --libs sdl2`
```

```
g++ -std=c++17 -I eigen/  -I /opt/homebrew/include -L/opt/homebrew/lib -lortools -labsl_log_internal_message -labsl_log_internal_check_op -labsl_strings -labsl_base -labsl_raw_logging_internal -labsl_log_internal_conditions polyhedron.cpp -o polyhedron
```

```
g++ -std=c++17 -I eigen/ polyhedron.cpp -o polyhedron `pkg-config --cflags --libs sdl2`

```

```
g++ -std=c++17 -I /opt/homebrew/include -L/opt/homebrew/lib -lortools test_ortools.cpp -o test_ortools
```

```
clang++ -O3 -Xpreprocessor -fopenmp -std=c++17 \
  -I/opt/homebrew/include/eigen3 \
  -I/opt/homebrew/opt/libomp/include \
  -L/opt/homebrew/opt/libomp/lib -lomp \
  polyhedron_parallel.cpp -o polyhedron_parallel
```
