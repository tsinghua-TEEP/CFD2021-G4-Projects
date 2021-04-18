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

using LinearAlgebra

"""
Work in progress.
"""
function inverted_poisson_2d_jacobi_step(xs::Matrix{T}, ys::Matrix{T}) where {T <: Number}
    (axes(xs) == axes(ys)) ||
        throw(DimensionMismatch("Dimensions of x ($(axes(xs))) and y ($(axes(ys))) does not match."))
    output_xs = similar(xs)
    output_ys = similar(xs)
    inverted_poisson_2d_jacobi_step!(output_xs, output_ys, xs, ys)
    return output_xs, output_ys
end
function inverted_poisson_2d_jacobi_step!(output_xs::Matrix{T}, output_ys::Matrix{T}, xs::Matrix{T}, ys::Matrix{T}) where {T <: Number}
    (axes(output_xs) == axes(output_ys) == axes(xs) == axes(ys)) ||
        throw(DimensionMismatch("Dimensions of input ($((axes(xs), axes(ys))) and output ($((axes(output_xs), axes(output_ys)))) are not compatible"))
    M = [l-1 for l in size(xs)]

    begin # neighboring x's
        x_l(u::CartesianIndex) = u[1] == axes(xs)[1][begin] ? zero(eltype(xs)) : xs[u-CartesianIndex(1,0)]
        x_r(u::CartesianIndex) = u[1] == axes(xs)[1][ end ] ? zero(eltype(xs)) : xs[u+CartesianIndex(1,0)]
        x_d(u::CartesianIndex) = u[2] == axes(xs)[2][begin] ? zero(eltype(xs)) : xs[u-CartesianIndex(0,1)]
        x_u(u::CartesianIndex) = u[2] == axes(xs)[2][ end ] ? zero(eltype(xs)) : xs[u+CartesianIndex(0,1)]
        xld(u::CartesianIndex) = u[1] == axes(xs)[1][begin] || u[2] == axes(xs)[2][begin] ? zero(eltype(xs)) : xs[u+CartesianIndex(-1,-1)]
        xrd(u::CartesianIndex) = u[1] == axes(xs)[1][ end ] || u[2] == axes(xs)[2][begin] ? zero(eltype(xs)) : xs[u+CartesianIndex( 1,-1)]
        xlu(u::CartesianIndex) = u[2] == axes(xs)[1][begin] || u[2] == axes(xs)[2][ end ] ? zero(eltype(xs)) : xs[u+CartesianIndex(-1, 1)]
        xru(u::CartesianIndex) = u[2] == axes(xs)[1][ end ] || u[2] == axes(xs)[2][ end ] ? zero(eltype(xs)) : xs[u+CartesianIndex( 1, 1)]
    end
    begin # neiboring y's
        y_l(u::CartesianIndex) = u[1] == axes(ys)[1][begin] ? zero(eltype(ys)) : ys[u-CartesianIndex(1,0)]
        y_r(u::CartesianIndex) = u[1] == axes(ys)[1][ end ] ? zero(eltype(ys)) : ys[u+CartesianIndex(1,0)]
        y_d(u::CartesianIndex) = u[2] == axes(ys)[2][begin] ? zero(eltype(ys)) : ys[u-CartesianIndex(0,1)]
        y_u(u::CartesianIndex) = u[2] == axes(ys)[2][ end ] ? zero(eltype(ys)) : ys[u+CartesianIndex(0,1)]
        yld(u::CartesianIndex) = u[1] == axes(ys)[1][begin] || u[2] == axes(ys)[2][begin] ? zero(eltype(ys)) : ys[u+CartesianIndex(-1,-1)]
        yrd(u::CartesianIndex) = u[1] == axes(ys)[1][ end ] || u[2] == axes(ys)[2][begin] ? zero(eltype(ys)) : ys[u+CartesianIndex( 1,-1)]
        ylu(u::CartesianIndex) = u[2] == axes(ys)[1][begin] || u[2] == axes(ys)[2][ end ] ? zero(eltype(ys)) : ys[u+CartesianIndex(-1, 1)]
        yru(u::CartesianIndex) = u[2] == axes(ys)[1][ end ] || u[2] == axes(ys)[2][ end ] ? zero(eltype(ys)) : ys[u+CartesianIndex( 1, 1)]
    end
    begin # 1-order diffs
        ∂xξ(u::CartesianIndex) = M[1]/2 * (x_r(u) - x_l(u))
        ∂xη(u::CartesianIndex) = M[2]/2 * (x_u(u) - x_d(u))
        ∂yξ(u::CartesianIndex) = M[1]/2 * (y_r(u) - y_l(u))
        ∂yη(u::CartesianIndex) = M[2]/2 * (y_u(u) - y_d(u))
    end
    begin # 2-order diffs
        α(u::CartesianIndex) = ∂xη(u)^2 + ∂yη(u)^2
        β(u::CartesianIndex) = ∂xξ(u)*∂xη(u) + ∂yξ(u)*∂yη(u)
        γ(u::CartesianIndex) = ∂xξ(u)^2 + ∂yξ(u)^2
        bh(u::CartesianIndex) = M[1]^2 * α(u)
        bv(u::CartesianIndex) = M[2]^2 * γ(u)
        bp(u::CartesianIndex) = 2 * (bh(u) + bv(u))
        cx(u::CartesianIndex) = - M[1]*M[2]/2 * β(u) * (xld(u) - xrd(u) - xlu(u) + xru(u))
        cy(u::CartesianIndex) = - M[1]*M[2]/2 * β(u) * (yld(u) - yrd(u) - ylu(u) + yru(u))
    end
    for u in CartesianIndices(xs)
        output_xs[u] = (bh(u)*(x_l(u) + x_r(u)) + bv(u)*(x_d(u) + x_u(u)) + cx(u)) / bp(u)
        output_ys[u] = (bh(u)*(y_l(u) + y_r(u)) + bv(u)*(y_d(u) + y_u(u)) + cy(u)) / bp(u)
    end
    return output_xs, output_ys
end

"""
Work in progress.
"""
function inverted_poisson_2d_jacobi_iterate(xs::Matrix{T}, ys::Matrix{T}, ε=1e-10)
    cache_xs = xs
    cache_ys = ys
    inverted_poisson_2d_jacobi_iterate!(cache_xs, cache_ys)
    return cache_xs, cache_ys
end
function inverted_poisson_2d_jacobi_iterate!(xs::Matrix{T}, ys::Matrix{T}, ε=1e-10)
    cache_xs = similar(xs)
    cache_ys = similar(xs)
    residue() = max(norm(cache_xs - xs), norm(cache_ys - ys))
    while residue() > ε
        inverted_poisson_2d_jacobi_step!(cache_xs, cache_ys, xs, ys)
        inverted_poisson_2d_jacobi_step!(xs, ys, cache_xs, cache_ys)
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
