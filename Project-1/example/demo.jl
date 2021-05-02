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
# using Revise
using BenchmarkTools
using Plots
using OffsetArrays

myblue = RGB(0.368417, 0.506779, 0.709798) # <==> 94, 129, 181 (the default style in Wolfram Mathematica)

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
Base.@pure function yt0012(ξ::T where {T<:Real})
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
        11.579398976610927  *     (1+ε - ξ)^2 +
        12.21202728162769   *     (1+ε - ξ)^3
        else
            0
        end)
end

const Mx, My = 100+1, 20+1

include(normpath(joinpath(@__DIR__, "../../src/misc-util.jl")))
using .__CFD2021__misc_util__: tuplejoin
"""
- usage of ``tuplejoin``:
```jldoctest
  julia> temp_v = ([0,0.5,1], [0,0.5,1])
  julia> temp_t = (axes(v) for v in temp_v)
  julia> @btime tuplejoin(temp_t...)
```
"""

include(normpath(joinpath(@__DIR__, "../src/transfinite-interpolate.jl")))
import .transfinite_interpolate: transfinite_interpolate_2d!, transfinite_interpolate_2d

"""
- Try to use the in-place version when applicable;
  there is a tremendous performance boost.
```jldoctest
  julia> v_lo = ([LinRange(0, 1, 10^4)...], [LinRange(0, 1, 10^4)...]); v_hi = ([LinRange(1, 2, 10^4)...], [LinRange(1, 2, 10^4)...])
  
  julia> interpolated_array = @btime Array{Float64}(undef, (length(v) for v in v_lo)...)
  20.500 μs (7 allocations: 762.94 MiB)
  
  julia> @btime transfinite_interpolate_2d!(interpolated_array, v_lo, v_hi)
  1.535 s (1 allocation: 96 bytes)
  
  julia> res = @btime transfinite_interpolate_2d(v_lo, v_hi)
  7.581 s (100000005 allocations: 2.98 GiB)
```
"""

"""
- example of a cartesian grid;
```jldoctest
  x_lo = ([LinRange(0, 1, Mx)...], [LinRange(0, 0, My)...]); x_hi = ([LinRange(0, 1, Mx)...], [LinRange(1, 1, My)...]);
  y_lo = ([LinRange(0, 0, Mx)...], [LinRange(0, 1, My)...]); y_hi = ([LinRange(1, 1, Mx)...], [LinRange(0, 1, My)...]);
```
"""

"""
- example of a half-O grid aroundNACA0012;
"""
x_lo = ([sinpi(ξ)^2 for ξ in LinRange(0, 0.5, Mx)], -exp.(LinRange(0, log(6), My)) .+ 1);
y_lo = (yt0012.(x_lo[1]), [LinRange(0, 0, My)...]);
x_hi = (5cospi.(LinRange(1, 0, Mx)), +exp.(LinRange(0, log(5), My)));
y_hi = (5sinpi.(LinRange(1, 0, Mx)), [LinRange(0, 0, My)...]);

"""
- example of a whole-O grid aroundNACA0012;
```jldoctest
  x_lo = ([         sinpi.(LinRange(-.5,0, Mx)).^2;          sinpi.(LinRange(0, .5, Mx)).^2],
             exp.(LinRange(0, log(5), My)                  )   );
  y_lo = ([-yt0012.(sinpi.(LinRange(-.5,0, Mx)).^2); yt0012.(sinpi.(LinRange(0, .5, Mx)).^2)],
          [LinRange(0, 0, My)...]);
  x_hi = (5cospi.(LinRange( 2, 0, 2Mx)),
             exp.(LinRange(0, log(5), My)                  )   );
  y_hi = (5sinpi.(LinRange( 2, 0, 2Mx)),
          [LinRange(0, 0, My)...]);
```
"""

p = plot(; xlims=(-5, 5), ylims=(-5, 5), aspect_ratio=:equal)
# p = plot(xlims = (-.6, .6), ylims = (0, .5), aspect_ratio = :equal)
interpolated_x = Array{Float64}(undef, (length(v) for v in x_lo)...);
interpolated_y = Array{Float64}(undef, (length(v) for v in y_lo)...);
@btime transfinite_interpolate_2d!(interpolated_x, x_lo, x_hi)
@btime transfinite_interpolate_2d!(interpolated_y, y_lo, y_hi)
for u in axes(interpolated_x)[begin]  plot!(p,interpolated_x[u,:], interpolated_y[u,:], label = nothing, color = myblue); end
for u in axes(interpolated_x)[ end ]  plot!(p,interpolated_x[:,u], interpolated_y[:,u], label = nothing, color = myblue); end
# for u in axes(interpolated_x)[begin]  plot!(p,interpolated_x[u,:],-interpolated_y[u,:], label = nothing, color = myblue); end
# for u in axes(interpolated_x)[ end ]  plot!(p,interpolated_x[:,u],-interpolated_y[:,u], label = nothing, color = myblue); end
plot!(p, yt0012, 0, 1; label=nothing, color=:red)
display(p) # savefig(normpath(joinpath(@__DIR__, "img/tf-hO.svg")))
