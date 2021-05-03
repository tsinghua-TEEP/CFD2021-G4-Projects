# CFD2021Projects API Documentation

```@meta
CurrentModule = CFD2021Projects
DocTestSetup = quote
    using Pkg
    Pkg.activate(normpath(joinpath(@__DIR__, "../..")))
    using CFD2021Projects
    include(normpath(joinpath(@__DIR__, "../../src/misc-util.jl")))
    include(normpath(joinpath(@__DIR__, "../../Project-1/src/transfinite-interpolate.jl")))
    using .transfinite_interpolate
end
```

## Functions

```@autodocs
Modules = [CFD2021Projects, CFD2021Projects.transfinite_interpolate]
Order   = [:function, :type]
```

## Index

```@index
```

```@meta
DocTestSetup = nothing
```
