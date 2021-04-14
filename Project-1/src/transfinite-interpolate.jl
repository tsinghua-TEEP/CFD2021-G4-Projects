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

using StaticArrays
using OffsetArrays
include("../../src/misc-util.jl")
using .__CFD2021__misc_util__: tuplejoin

"""
    transfinite_interpolate_2d(f_lo, f_hi[, v_cr])

Return a interpolating function defined on ``[0, 1]²``.

**Notice**: does not check or warn when extrapolating.

See also: https://en.wikipedia.org/wiki/Transfinite_interpolation

TODO: generalize to higher dimensions.

# Arguments
- `f_lo::Tuple{Function, Function}`: a pair of functions.
- `f_hi::Tuple{Function, Function}`: a pair of functions.
- `v_cr::SMatrix{2, 2} = nothing`: (optional) values to take at the corners;
    use average of edge functions if not provided.
    Still, in practical use you should make sure the ends of the edges meet.

# Output
- `interpolating_function(u::NTuple{2, T} where T <: Number):: Number`:
    Take in a length-2 static vector of numbers and return a number.

# Examples
```jldoctest
julia> interp_func = transfinite_interpolate_2d( (x->2x, y->3y), (x->4x+3, y->2+5y) )
(::var"#interpolating_function#1"{Tuple{var"#2#6", var"#3#7"}, Tuple{var"#4#8", var"#5#9"}}) (generic function with 1 method)

julia> interp_func((0.5,0.5))
3.0
```
"""
function transfinite_interpolate_2d(
    f_lo::Tuple{Function, Function},
    f_hi::Tuple{Function, Function},
    v_cr::Union{Nothing, SMatrix{2, 2}} = nothing)::Function
    #= handling of optionals suggested by
       https://stackoverflow.com/questions/42499528/julia-convention-for-optional-arguments;
       though Base.Void is deprecated and now called Nothing instead, as pointed out by
       https://discourse.julialang.org/t/void-not-defined-version-problem/14183 =#
    if typeof(v_cr) <: Nothing
        v_cr = SA[(f_lo[1](  0  )+f_lo[2](  0  ))/2 (f_lo[1](  1  )+f_hi[2](  0  ))/2
                  (f_hi[1](  0  )+f_lo[2](  1  ))/2 (f_hi[1](  1  )+f_hi[2](  1  ))/2]
    end
    
    function interpolating_function(u::NTuple{2, T} where T <: Number):: Number
     #= linear-part from edges   bilinear-part from vertices =#
        (1-u[2])*f_lo[1](u[1]) - (1-u[1])*(1-u[2])*v_cr[1,1] +
        (  u[2])*f_hi[1](u[1]) - (  u[1])*(1-u[2])*v_cr[1,2] +
        (1-u[1])*f_lo[2](u[2]) - (1-u[1])*(  u[2])*v_cr[2,1] +
        (  u[1])*f_hi[2](u[2]) - (  u[1])*(  u[2])*v_cr[2,2]
    end
    return interpolating_function
end

"""
    transfinite_interpolate_2d(v_lo, v_hi[, v_cr])
    transfinite_interpolate_2d!(interpolated_array, v_lo, v_hi[, v_cr])

Return a interpolated array defined on a discretization of ``[0, 1]²``.

**Notice**: the in-place version does not check dimension-compatibility of its inputs.

TODO: generalize to higher dimensions.

# Arguments
- `v_lo::Tuple{Vector{<: Real}, Vector{<: Real}}`: a pair of vectors.
- `v_hi::Tuple{Vector{<: Real}, Vector{<: Real}}`: a pair of vectors.
- `v_cr::SMatrix{2, 2, <: Real} = nothing`: (optional) values to take at the corners;
    use average of edge values if not provided.
    Still, in practical use you should make sure the ends of the edges meet.

# Output
- `interpolated_array::Array{Real, 2}`

# Examples
```jldoctest
julia> transfinite_interpolate_2d(([0,0.5,1], [0,0.5,1]), ([1,1.5,2], [1,1.5,2]))
3×3 Matrix{Real}:
 0.0  0.5  1.0
 0.5  1.0  1.5
 1.0  1.5  2.0
```
"""
function transfinite_interpolate_2d(
    v_lo::Tuple{Vector{T}, Vector{S}} where {T <: Real, S <: Real},
    v_hi::Tuple{Vector{T}, Vector{S}} where {T <: Real, S <: Real},
    v_cr::Union{Nothing, SMatrix{2, 2, T} where T <: Real} = nothing)::Array{Real, 2}

    #= check dimension-compatibility before allocating anything =#
    length(v_lo) == length(v_hi) || throw(
        DimensionMismatch(" incompatible grid specification ($(length(v_lo))-d at one end but $(length(v_hi))-d at the other)."))
    N = length(v_lo); errmsg = "" #= In this function, dimensions are always indexed from 1, but not so for each axis of them =#
    for dim in eachindex(v_lo)  if ! (axes(v_lo[dim]) == axes(v_hi[dim]))
        errmsg *= "grid mismatch in dimension $(dim) ($(axes(v_lo[dim])) at one end but $(axes(v_hi[dim])) at the other)."
    end end; errmsg == "" || throw(DimensionMismatch(errmsg))

    interpolated_array = Array{Float64}(undef, (length(v) for v in v_lo)...)
    interpolated_array = OffsetArray(interpolated_array, tuplejoin((axes(v) for v in v_lo)...))
    transfinite_interpolate_2d!(interpolated_array, v_lo, v_hi, v_cr)
    return interpolated_array
end
function transfinite_interpolate_2d!(
    interpolated_array::OffsetArray,
    v_lo::Tuple{Vector{T}, Vector{S}} where {T <: Real, S <: Real},
    v_hi::Tuple{Vector{T}, Vector{S}} where {T <: Real, S <: Real},
    v_cr::Union{Nothing, SMatrix{2, 2, T} where T <: Real} = nothing)::Array{Real, 2}

    if typeof(v_cr) <: Nothing # e.g, the following lines still assumes range ``1:2``.
        v_cr = SA[(v_lo[1][begin]+v_lo[2][begin])/2 (v_lo[1][ end ]+v_hi[2][begin])/2
                  (v_hi[1][begin]+v_lo[2][ end ])/2 (v_hi[1][ end ]+v_hi[2][ end ])/2]
    end; M = [length(v)-1 for v in v_lo] # this is only for scaling, and never for iterating

    # loop Version
    for u in CartesianIndices(interpolated_array) #= #= Deprecated =# Iterators.product(IndexCartesian(), [eachindex(v) for v in v_lo]...) =#
      interpolated_array[u] =
      #= linear-part   from   edges      bilinear-part   from   vertices    =#
        (1-u[2]/M[2])*v_lo[1][u[1]] - (1-u[1]/M[1])*(1-u[2]/M[2])*v_cr[1,1] +
        (  u[2]/M[2])*v_hi[1][u[1]] - (  u[1]/M[1])*(1-u[2]/M[2])*v_cr[1,2] +
        (1-u[1]/M[1])*v_lo[2][u[2]] - (1-u[1]/M[1])*(  u[2]/M[2])*v_cr[2,1] +
        (  u[1]/M[1])*v_hi[2][u[2]] - (  u[1]/M[1])*(  u[2]/M[2])*v_cr[2,2]
    end
    return interpolated_array
end
