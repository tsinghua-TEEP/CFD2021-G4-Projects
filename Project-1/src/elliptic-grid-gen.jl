# Copyright 2021 Gravifer
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#     http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
#=
  Author: Gravifer, X. H. Zhu, H. J. Yan, H. Y. Sun
  Date: 2021-04-10
  Description: TBD
=#

module elliptic_grid_gen

using LinearAlgebra, ProgressMeter

"""
    inverted_poisson_2d_jacobi_step(xs, ys[, P, Q])
    inverted_poisson_2d_jacobi_step!(output_xs, output_ys, xs, ys[, P, Q])

Take a Jacobi step for an array-pair represented 2d grid.

TODO: generalize to higher dimensions.

# Arguments
- `xs::Matrix{<: Number}`: array of `x`-coordinates.
- `ys::Matrix{<: Number}`: array of `y`-coordinates.
- `P::Function`: (optional) `ξ` related source term.
- `Q::Function`: (optional) `η` related source term.

# Output
- `output_xs::Matrix{<: Number}`: array of `x`-coordinates, at the next step.
- `output_ys::Matrix{<: Number}`: array of `y`-coordinates, at the next step.

# Examples
```jldoctest
```
"""
function inverted_poisson_2d_jacobi_step(
    xs::Matrix{T} , ys::Matrix{T} ,
    P::S = nothing, Q::S = nothing            ) where {T <: Number, S <: Union{Nothing, Function}}

    (axes(xs) == axes(ys)) ||
        throw(DimensionMismatch("Dimensions of x ($(axes(xs))) and y ($(axes(ys))) does not match."))
    output_xs = similar(xs)
    output_ys = similar(ys)
    inverted_poisson_2d_jacobi_step!(output_xs, output_ys, xs, ys, P, Q)
    return output_xs, output_ys
end
function inverted_poisson_2d_jacobi_step!(
    output_xs::Matrix{T}, output_ys::Matrix{T},
           xs::Matrix{T},        ys::Matrix{T},
    P::S = nothing      , Q::S = nothing      ) where {T <: Number, S <: Union{Nothing, Function}}

    (axes(output_xs) == axes(output_ys) == axes(xs) == axes(ys)) ||
        throw(DimensionMismatch("Dimensions of input ($((axes(xs), axes(ys))) and output ($((axes(output_xs), axes(output_ys)))) are not compatible"))
    M = [l-1 for l in size(xs)]
    "neighboring x's"; begin
        @inline x_l(u::CartesianIndex) = u[1] == axes(xs)[1][begin]                               ? zero(eltype(xs)) : xs[u+CartesianIndex(-1, 0)]
        @inline x_r(u::CartesianIndex) = u[1] == axes(xs)[1][ end ]                               ? zero(eltype(xs)) : xs[u+CartesianIndex( 1, 0)]
        @inline x_d(u::CartesianIndex) =                               u[2] == axes(xs)[2][begin] ? zero(eltype(xs)) : xs[u+CartesianIndex( 0,-1)]
        @inline x_u(u::CartesianIndex) =                               u[2] == axes(xs)[2][ end ] ? zero(eltype(xs)) : xs[u+CartesianIndex( 0, 1)]
        @inline xld(u::CartesianIndex) = u[1] == axes(xs)[1][begin] || u[2] == axes(xs)[2][begin] ? zero(eltype(xs)) : xs[u+CartesianIndex(-1,-1)]
        @inline xrd(u::CartesianIndex) = u[1] == axes(xs)[1][ end ] || u[2] == axes(xs)[2][begin] ? zero(eltype(xs)) : xs[u+CartesianIndex( 1,-1)]
        @inline xlu(u::CartesianIndex) = u[2] == axes(xs)[1][begin] || u[2] == axes(xs)[2][ end ] ? zero(eltype(xs)) : xs[u+CartesianIndex(-1, 1)]
        @inline xru(u::CartesianIndex) = u[2] == axes(xs)[1][ end ] || u[2] == axes(xs)[2][ end ] ? zero(eltype(xs)) : xs[u+CartesianIndex( 1, 1)]
    end
    "neighboring y's"; begin
        @inline y_l(u::CartesianIndex) = u[1] == axes(ys)[1][begin]                               ? zero(eltype(ys)) : ys[u+CartesianIndex(-1, 0)]
        @inline y_r(u::CartesianIndex) = u[1] == axes(ys)[1][ end ]                               ? zero(eltype(ys)) : ys[u+CartesianIndex( 1, 0)]
        @inline y_d(u::CartesianIndex) =                               u[2] == axes(ys)[2][begin] ? zero(eltype(ys)) : ys[u+CartesianIndex( 0,-1)]
        @inline y_u(u::CartesianIndex) =                               u[2] == axes(ys)[2][ end ] ? zero(eltype(ys)) : ys[u+CartesianIndex( 0, 1)]
        @inline yld(u::CartesianIndex) = u[1] == axes(ys)[1][begin] || u[2] == axes(ys)[2][begin] ? zero(eltype(ys)) : ys[u+CartesianIndex(-1,-1)]
        @inline yrd(u::CartesianIndex) = u[1] == axes(ys)[1][ end ] || u[2] == axes(ys)[2][begin] ? zero(eltype(ys)) : ys[u+CartesianIndex( 1,-1)]
        @inline ylu(u::CartesianIndex) = u[2] == axes(ys)[1][begin] || u[2] == axes(ys)[2][ end ] ? zero(eltype(ys)) : ys[u+CartesianIndex(-1, 1)]
        @inline yru(u::CartesianIndex) = u[2] == axes(ys)[1][ end ] || u[2] == axes(ys)[2][ end ] ? zero(eltype(ys)) : ys[u+CartesianIndex( 1, 1)]
    end
    "1-order diffs"; begin
        @inline ∂xξ(u::CartesianIndex) = M[1]/2 * (x_r(u) - x_l(u))
        @inline ∂xη(u::CartesianIndex) = M[2]/2 * (x_u(u) - x_d(u))
        @inline ∂yξ(u::CartesianIndex) = M[1]/2 * (y_r(u) - y_l(u))
        @inline ∂yη(u::CartesianIndex) = M[2]/2 * (y_u(u) - y_d(u))
    end
    "2-order diffs"; begin
        @inline α(  u::CartesianIndex) =  ∂xη(u)^2      + ∂yη(u)^2
        @inline β(  u::CartesianIndex) =  ∂xξ(u)*∂xη(u) + ∂yξ(u)*∂yη(u)
        @inline γ(  u::CartesianIndex) =  ∂xξ(u)^2      + ∂yξ(u)^2
        @inline bh( u::CartesianIndex) =  M[1]^2 * α(u)
        @inline bv( u::CartesianIndex) =  M[2]^2 * γ(u)
        @inline bp( u::CartesianIndex) = (bh(u) + bv(u)) * 2
        @inline cx( u::CartesianIndex) = -M[1]*M[2]/2 * β(u) * (xld(u) - xrd(u) - xlu(u) + xru(u))
        @inline cy( u::CartesianIndex) = -M[1]*M[2]/2 * β(u) * (yld(u) - yrd(u) - ylu(u) + yru(u))
    end
    "source related terms"; begin
        @inline sx( u::CartesianIndex) = (typeof(P) <: Nothing ? 0 : α(u)*P(u)*∂xξ(u)) + (typeof(Q) <: Nothing ? 0 : γ(u)*Q(u)*∂xη(u))
        @inline sy( u::CartesianIndex) = (typeof(P) <: Nothing ? 0 : α(u)*P(u)*∂yξ(u)) + (typeof(Q) <: Nothing ? 0 : γ(u)*Q(u)*∂yη(u))
    end
    "actual computation"
    for u in CartesianIndices(xs)[[begin end], :]
        output_xs[u] = xs[u]
        output_ys[u] = ys[u]
    end
    for u in CartesianIndices(xs)[:, [begin end]]
        output_xs[u] = xs[u]
        output_ys[u] = ys[u]
    end
    Threads.@threads for u in CartesianIndices(xs)[begin+1:end-1, begin+1:end-1] # this is good enough for a fixed dirichlet boundary
        output_xs[u] = (bh(u)*(x_l(u) + x_r(u)) + bv(u)*(x_d(u) + x_u(u)) + cx(u) + sx(u)) / bp(u)
        output_ys[u] = (bh(u)*(y_l(u) + y_r(u)) + bv(u)*(y_d(u) + y_u(u)) + cy(u) + sy(u)) / bp(u)
    end
    return output_xs, output_ys
end

"""
    inverted_poisson_2d_jacobi_iterate(xs, ys[, P, Q])
    inverted_poisson_2d_jacobi_iterate!(xs, ys[, P, Q[, cache_xs, cache_ys]])

Full Jacobi relaxation of an array-pair represented 2d grid.

TODO: generalize to higher dimensions.

# Arguments
- `xs::Matrix{<: Number}`: array of `x`-coordinates.
- `ys::Matrix{<: Number}`: array of `y`-coordinates.
- `ε::Float64=1e-10`: (optional) tolerance of residue.
- `P::Function`: (optional) `ξ` related source term.
- `Q::Function`: (optional) `η` related source term.
- `cache_xs::Matrix{<: Number}`: (optional) buffer array of `x`-coordinates.
- `cache_ys::Matrix{<: Number}`: (optional) buffer array of `y`-coordinates.

# Output
- `output_xs::Matrix{<: Number}`: relaxed array of `x`-coordinates.
- `output_ys::Matrix{<: Number}`: relaxed array of `y`-coordinates.

# Examples
```jldoctest
```
"""
function inverted_poisson_2d_jacobi_iterate(
    xs::Matrix{T} , ys::Matrix{T} , ε::Float64=1e-10,
    P::S = nothing, Q::S = nothing              ) where {T <: Number, S <: Union{Nothing, Function}}

    cache_xs = xs
    cache_ys = ys
    inverted_poisson_2d_jacobi_iterate!(cache_xs, cache_ys, ε, P, Q)
    return cache_xs, cache_ys
end
function inverted_poisson_2d_jacobi_iterate!(
    xs::Matrix{T} , ys::Matrix{T} , ε::Float64=1e-10,
    P::S = nothing, Q::S = nothing,
    cache_xs::R = nothing, cache_ys::R = nothing) where {T <: Number, S <: Union{Nothing, Function}, R <: Union{Nothing, Matrix{T}}}

    !(typeof(cache_xs) <: Nothing) || (cache_xs = similar(xs))
    !(typeof(cache_ys) <: Nothing) || (cache_ys = similar(xs))
    @inline residue() = max(norm(cache_xs - xs), norm(cache_ys - ys))
    @inline steppair()=(inverted_poisson_2d_jacobi_step!(cache_xs, cache_ys, xs, ys, P, Q);
                        inverted_poisson_2d_jacobi_step!(xs, ys, cache_xs, cache_ys, P, Q))
    steppair()
    prog = ProgressThresh(ε, "Elliptic relaxation:")
    while residue() > ε
        ProgressMeter.update!(prog, residue())
        sleep(0.1)
        steppair()
    end
    return xs, ys
end

"""
To be done.
"""
function dirichlet_boundary_conform end

"""
To be done.
"""
function nuemann_boundary_conform end

end # module elliptic_grid_gen
