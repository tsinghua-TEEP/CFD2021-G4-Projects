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

using Plots

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

myblue = RGB(0.368417, 0.506779, 0.709798) # <==> 94, 129, 181 (the default style in Wolfram Mathematica)

function grid_display(
    xs::Matrix{T} , ys::Matrix{T};
    xlims=(-5, 5), ylims=(-5, 5), aspect_ratio=:equal, color=myblue, title=nothing, mirror::Bool=false) where {T <: Number}

    (axes(xs) == axes(ys)) ||
        throw(DimensionMismatch("Dimensions of x ($(axes(xs))) and y ($(axes(ys))) does not match."))
    p = plot(; xlims=xlims, ylims=ylims, aspect_ratio=aspect_ratio)
    p = grid_display!(p, xs, ys; xlims=xlims, ylims=ylims, aspect_ratio=aspect_ratio, color=color, title=title, mirror=mirror)
    return p
end

function grid_display!(plot_ref,
    xs::Matrix{T} , ys::Matrix{T};
    xlims=(-5, 5), ylims=(-5, 5), aspect_ratio=:equal, color=myblue, title=nothing, mirror::Bool=false) where {T <: Number}

    (axes(xs) == axes(ys)) ||
        throw(DimensionMismatch("Dimensions of x ($(axes(xs))) and y ($(axes(ys))) does not match."))
    p = (typeof(plot_ref) <: Nothing) ? plot(; xlims=xlims, ylims=ylims, aspect_ratio=aspect_ratio) : plot_ref
    "normal"; begin
        for u in axes(xs)[begin]  plot!(p, interpolated_x[u,:], interpolated_y[u,:], label=nothing, color=color); end
        for u in axes(ys)[ end ]  plot!(p, interpolated_x[:,u], interpolated_y[:,u], label=nothing, color=color); end
    end
    if mirror
        for u in axes(xs)[begin]  plot!(p, interpolated_x[u,:],-interpolated_y[u,:], label=nothing, color=color); end
        for u in axes(ys)[ end ]  plot!(p, interpolated_x[:,u],-interpolated_y[:,u], label=nothing, color=color); end
    end
    plot!(p,        yt0012   , 0, 1; label=nothing, color=:black)
    plot!(p, (ξ)-> -yt0012(ξ), 0, 1; label=nothing, color=:black)
    title!(p, title)
    display(p)
    return p
end
