# Soft B

## Demo

![ezgif-5-a8310b6e8a](https://user-images.githubusercontent.com/80531783/234442338-22096d73-cd47-4a99-86e6-4d7b239c00ef.gif)

## Build on Microsoft Windows with Visual Studio 2017/2019

### Instruction

- Open SoftSim.sln
- Build
- Executable will be in ./bin

## Build on other platforms and/or compilers

### Prerequisite

- [Cmake](https://cmake.org) (version >= 3.14)
- Compiler (e.g. GCC)

### Instruction

- Run:

```bash=
cmake -S . -B build
cmake --build build --config Release --target install --parallel 8
```
- Executable will be in ./bin

If you are building on Linux, you need one of these dependencies, usually `xorg-dev`

- `xorg-dev` (For X11)
- `libwayland-dev wayland-protocols extra-cmake-modules libxkbcommon-dev` (For Wayland)
- `libosmesa6-dev` (For OSMesa)
