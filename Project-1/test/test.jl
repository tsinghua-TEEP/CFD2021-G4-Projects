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
using BenchmarkTools

# include("../../src/misc-util.jl")
# using .__CFD2021__misc_util__: tuplejoin
# temp_v = ([0,0.5,1], [0,0.5,1])
# temp_t = (axes(v) for v in temp_v)
# @btime tuplejoin(temp_t...)

include("../src/transfinite-interpolate.jl")
res = @btime transfinite_interpolate_2d(([0,0.5,1], [0,0.5,1]), ([1,1.5,2], [1,1.5,2]))
@show res