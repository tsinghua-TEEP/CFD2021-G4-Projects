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

using Test
using BenchmarkTools
using OffsetArrays

const Mx, My = 3, 3

include(normpath(joinpath(@__DIR__, "../../src/misc-util.jl")))
using .__CFD2021__misc_util__: tuplejoin

include("../src/transfinite-interpolate.jl")
import .transfinite_interpolate: transfinite_interpolate_2d!, transfinite_interpolate_2d

# normal group
x_lo = ([LinRange(0, 1, Mx)...], [LinRange(0, 0, My)...]); x_hi = ([LinRange(0, 1, Mx)...], [LinRange(1, 1, My)...]);
y_lo = ([LinRange(0, 0, Mx)...], [LinRange(0, 1, My)...]); y_hi = ([LinRange(1, 1, Mx)...], [LinRange(0, 1, My)...]);
in_situ_interp_x = Array{Float64}(undef, (length(v) for v in x_lo)...);
in_situ_interp_y = Array{Float64}(undef, (length(v) for v in y_lo)...);
transfinite_interpolate_2d!(in_situ_interp_x, x_lo, x_hi)
transfinite_interpolate_2d!(in_situ_interp_y, y_lo, y_hi)
ab_situ_interp_x = transfinite_interpolate_2d(x_lo, x_hi)
ab_situ_interp_y = transfinite_interpolate_2d(y_lo, y_hi)

@test in_situ_interp_x == ab_situ_interp_x == [0.0 0.0 0.0; 0.5 0.5 0.5; 1.0 1.0 1.0]
@test in_situ_interp_y == ab_situ_interp_y == [0.0 0.5 1.0; 0.0 0.5 1.0; 0.0 0.5 1.0]

x_lo = map(x->OffsetArray(x, OffsetArrays.Origin(0)), x_lo)
y_lo = map(x->OffsetArray(x, OffsetArrays.Origin(0)), y_lo)
x_hi = map(x->OffsetArray(x, OffsetArrays.Origin(0)), x_hi)
y_hi = map(x->OffsetArray(x, OffsetArrays.Origin(0)), y_hi)
in_situ_interp_x = OffsetArray(in_situ_interp_x, tuplejoin((axes(v) for v in x_lo)...))
in_situ_interp_y = OffsetArray(in_situ_interp_y, tuplejoin((axes(v) for v in y_lo)...))
transfinite_interpolate_2d!(in_situ_interp_x, x_lo, x_hi)
transfinite_interpolate_2d!(in_situ_interp_y, y_lo, y_hi)
ab_situ_interp_x = transfinite_interpolate_2d(x_lo, x_hi)
ab_situ_interp_y = transfinite_interpolate_2d(y_lo, y_hi)

@test in_situ_interp_x == ab_situ_interp_x != [0.0 0.0 0.0; 0.5 0.5 0.5; 1.0 1.0 1.0]
@test in_situ_interp_y == ab_situ_interp_y != [0.0 0.5 1.0; 0.0 0.5 1.0; 0.0 0.5 1.0]
@test in_situ_interp_x == ab_situ_interp_x == OffsetArray([0.0 0.0 0.0; 0.5 0.5 0.5; 1.0 1.0 1.0], OffsetArrays.Origin(0))
@test in_situ_interp_y == ab_situ_interp_y == OffsetArray([0.0 0.5 1.0; 0.0 0.5 1.0; 0.0 0.5 1.0], OffsetArrays.Origin(0))

x_lo = ([LinRange(0, 1, Mx)...], [LinRange(0, 0, My)...]); x_hi = ([LinRange(0, 1, 2Mx+1)...], [LinRange(1, 1, 2My+1)...]);
@test_throws DimensionMismatch transfinite_interpolate_2d(x_lo, x_hi)
