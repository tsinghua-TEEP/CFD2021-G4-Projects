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

include(normpath(joinpath(@__DIR__, "../src/misc-util.jl")))
include(normpath(joinpath(@__DIR__, "../src/CFD2021Projects.jl")))
using .__CFD2021__misc_util__: tuplejoin
@test tuplejoin((1, 2)) == (1, 2)
temp_v = ([0,0.5,1], [0,0.5,1], [0,0.5,1], [0,0.5,1])
temp_t = (axes(v) for v in temp_v)
@test tuplejoin(temp_t...) == (Base.OneTo(3), Base.OneTo(3), Base.OneTo(3), Base.OneTo(3))

include(normpath(joinpath(@__DIR__, "../Project-1/test/runtests.jl")))
