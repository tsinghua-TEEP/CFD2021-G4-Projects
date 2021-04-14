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

"""
To be done.
"""
function inverted_poisson_jacobi_iterate(arr, ε=1e-10)
end

"""
To be done.
"""
function dirichlet_boundary_conform end

"""
To be done.
"""
function nuemann_boundary_conform end
