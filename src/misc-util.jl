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
module __CFD2021__misc_util__
# as suggested by https://discourse.julialang.org/t/efficient-tuple-concatenation/5398
@inline tuplejoin(x::Tuple, y::Tuple, z...) = (x..., tuplejoin(y, z...)...)
@inline tuplejoin(x::Tuple, y::Tuple) = (x..., y...)
@inline tuplejoin(t::Tuple) = t
end # module __CFD2021__misc_util__
