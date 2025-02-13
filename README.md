To compile `voxel_editor.cpp` use
```
g++ -std=c++14 -I eigen/ -o voxel_editor voxel_editor.cpp `sdl2-config --cflags --libs`
```
or
```
g++ -std=c++14 -I eigen/ -o voxel_editor voxel_editor.cpp `pkg-config --cflags --libs sdl2`
```

Eigen can be obtained from https://gitlab.com/libeigen/eigen
