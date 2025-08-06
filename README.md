# GDMathLib.jl: Game development mathematic kit

Ever wanted to build a game or a game engine but stumbled on tricky math concepts ?
You want to build your game but you don't have the time to learn about quaternions ?
Then this package is for you.
**GDMathLib.jl** takes away all the hassle and boilerplates of game maths and offer a high performance, ready to use toolkit.

## Installation

```julia
julia> ]add GDMathLib
```

For the development version

```julia
julia> ]add https://github.com/Gesee-y/GDMathLib.jl
```

## Features

### Math Foundations

* **Static Vectors** (`Vec2`, `Vec3`, `Vec4`, etc.): Allocated on the stack so no GC pressure.
* **Typed Hierarchy**: Unified `VectorType{T}` system including `Quaternion{T}`.
* **Optimized Performance**: Uses `@inbounds`, `@simd`, and loop unrolling.

### Utilities

* **Functions**: `lerp`, `inverse_lerp`, `projection_matrix`, etc.
* **Transformations**: The `Transformation{N}` struct lets you manage batched affine transforms efficiently, with manual matrix updates via `update_transform_matrix!`.

### Coordinate Systems

* Work in alternative spaces like `PolarCoord`, `SphericalCoord`, `CylindricalCoord`, etc.
* Convert easily with `ToCartesian(...)`.

### Colors

* Built-in color palette (`RED`, `GREEN`, `PURPLE`, etc.) ready to use with other libraries.

## Example

```julia
using GDMathLib

a = Vec3f(1,2,3)
b = Vec3f(7,5,9)

c = lerp(a, b, 0.5)
d = vnorm(a)
```

## License

This package is licensed under the MIT License.

## Bug report

If you encounter any bug or precision problem, don't hesitate to leave an issue.

It helps us improve the package.
