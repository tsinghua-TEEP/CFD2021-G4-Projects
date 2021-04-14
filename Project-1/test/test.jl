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

using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "../..")))
using Plots, OffsetArrays
using BenchmarkTools

"""
    yt0012(ξ)

Compute ``yₜ`` (i.e., the half-thickness relative to the chord)
of NACA0012 airfoils from the normalized chord position ``ξ=x/c``.

Additionally defined upon ``[1, 1.005]`` for a 2-order smooth closure.

See also: https://en.wikipedia.org/wiki/NACA_airfoil

TODO: generalize for 00xx airfoils; furthermore, for all 4-digit airfoils.

# Examples
```jldoctest
julia> yt0012(0.13)
0.051193437032115506
```
"""
function yt0012(ξ)
  ε = 0.0050544107511712
  return 5 * 0.12 *  (
  if    0 <= ξ <= 1
     0.2969 * sqrt(ξ)  -
     0.1260 *      ξ   -
     0.3516 *      ξ^2 +
     0.2843 *      ξ^3 -
     0.1015 *      ξ^4
  elseif 1 < ξ <= 1+ε # 2-order smooth joining of the tail edge.
     0.03423616658680458 * sqrt(1+ε - ξ)   +
     0.12491972188290168 *     (1+ε - ξ)   +
    11.579398976610927   *     (1+ε - ξ)^2 +
    12.21202728162769    *     (1+ε - ξ)^3
  else
     0
  end)
end


include("../../src/misc-util.jl")
using .__CFD2021__misc_util__: tuplejoin
# temp_v = ([0,0.5,1], [0,0.5,1])
# temp_t = (axes(v) for v in temp_v)
# @btime tuplejoin(temp_t...)

include("../src/transfinite-interpolate.jl")
import .transfinite_interpolate: transfinite_interpolate_2d!, transfinite_interpolate_2d
# v_lo = ([LinRange(0, 1, 10^4)...], [LinRange(0, 1, 10^4)...]); v_hi = ([LinRange(1, 2, 10^4)...], [LinRange(1, 2, 10^4)...])
# interpolated_array = @btime Array{Float64}(undef, (length(v) for v in v_lo)...) # 20.500 μs (7 allocations: 762.94 MiB)
# @btime transfinite_interpolate_2d!(interpolated_array, v_lo, v_hi) # 1.535 s (1 allocation: 96 bytes)
# res = @btime transfinite_interpolate_2d(v_lo, v_hi) # 7.581 s (100000005 allocations: 2.98 GiB)
# @show res