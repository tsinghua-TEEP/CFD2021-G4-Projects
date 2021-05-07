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
Pkg.activate(normpath(joinpath(@__DIR__, "..")))
using Documenter
using CFD2021Projects
doctest__setup__ = quote
    include(normpath(joinpath(@__DIR__, "../src/misc-util.jl")))
    include(normpath(joinpath(@__DIR__, "../Project-1/src/transfinite-interpolate.jl")))
    import .transfinite_interpolate: transfinite_interpolate_2d!, transfinite_interpolate_2d
end
eval(doctest__setup__)

DocMeta.setdocmeta!(CFD2021Projects, :DocTestSetup, doctest__setup__; recursive=true)

makedocs(sitename = "CFD2021Projects",
         modules  = [CFD2021Projects])
