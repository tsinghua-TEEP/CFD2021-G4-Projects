# CFD2021-G4-Projects

<!-- [![CMake](https://github.com/tsinghua-TEEP/CFD2021-G4-Projects/actions/workflows/cmake.yml/badge.svg)](https://github.com/tsinghua-TEEP/CFD2021-G4-Projects/actions/workflows/cmake.yml) !-->
<!-- [![Python package with Conda](https://github.com/tsinghua-TEEP/CFD2021-G4-Projects/actions/workflows/python-package-conda.yml/badge.svg)](https://github.com/tsinghua-TEEP/CFD2021-G4-Projects/actions/workflows/python-package-conda.yml) !-->
[![Codacy grade](https://img.shields.io/codacy/grade/8ddf95075915482d8708388554f16386?label=quality&logo=Codacy)](https://www.codacy.com?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=tsinghua-TEEP/CFD2021-G4-Projects&amp;utm_campaign=Badge_Grade)<!-- ![codacy](https://app.codacy.com/project/badge/Grade/8ddf95075915482d8708388554f16386?label=) !-->
[![CircleCI build status](https://img.shields.io/circleci/build/gh/tsinghua-TEEP/CFD2021-G4-Projects?label=build&logo=CircleCI&token=9b51e15e5ced695a347386f06bdc605e23e7d8e5)](https://github.com/tsinghua-TEEP/CFD2021-G4-Projects/actions)<!-- ![circleci](https://circleci.com/gh/tsinghua-TEEP/CFD2021-G4-Projects.svg?style=shield&label=CircleCI&logo=CircleCI&circle-token=9b51e15e5ced695a347386f06bdc605e23e7d8e5) !-->
[![CodeCov coverage](https://img.shields.io/codecov/c/gh/tsinghua-TEEP/CFD2021-G4-Projects?logo=CodeCov&logoColor=white&token=9R7SWYU9W5)](https://codecov.io/gh/tsinghua-TEEP/CFD2021-G4-Projects)<!-- ![codecov](https://codecov.io/gh/tsinghua-TEEP/CFD2021-G4-Projects/branch/main/graph/badge.svg?token=9R7SWYU9W5&logoColor=white&style=flat) !-->
[![License: Apache-2.0](https://img.shields.io/badge/license-APL2-blue.svg?logo=Apache&style=flat)](https://github.com/tsinghua-TEEP/CFD2021-G4-Projects/blob/main/LICENSE)
[![Wakatime timing](https://wakatime.com/badge/github/tsinghua-TEEP/CFD2021-G4-Projects.svg?logo=WakaTime&style=flat)](https://wakatime.com/badge/github/tsinghua-TEEP/CFD2021-G4-Projects)  
[![ColPrac: Contributor's Guide](https://img.shields.io/badge/ColPrac-Contributor's%20Guide-blueviolet?logo=Julia&logoColor=white&style=flat)](https://github.com/SciML/ColPrac)
[![GitHub Repository Discussions](https://img.shields.io/badge/GitHub-Discussions-blueviolet?logo=github&style=flat)](https://github.com/tsinghua-TEEP/CFD2021-G4-Projects/discussions)
[![Zoom meeting](https://img.shields.io/static/v1?logo=LiveChat&logoColor=white&label=meeting&message=Zoom&style=flat&color=2D8CFF)](https://us02web.zoom.us/j/88643726401?pwd=V3BNdTV4TWlvZmZkd2VoSHZ0Y2Q0Zz09)
<!-- [![Dependabot Status](https://api.dependabot.com/badges/status?host=github&repo=tsinghua-TEEP/CFD2021-G4-Projects)](https://dependabot.com) !-->

This repository dedicates to the team projects
related to the [computational fluid dynamics course](http://reserves.lib.tsinghua.edu.cn/Courses/CourseDetail?courseId=cb3f2412-7ba3-465f-a8b8-c24061b136d9) by [Prof. Ren](http://www.hy.tsinghua.edu.cn/info/1154/1826.htm), 2021 spring, at SAS, THU.  
Efforts are contributed by a joined 4-member group.  
For practice and performance, we attempt to use the
[![Julia](https://img.shields.io/static/v1?logo=Julia&logoColor=white&label=&message=Julia&color=9558B2)](https://julialang.org)
programming language in these projects.

## Configuration

The project is mainly written in
Julia<!-- [![Julia](https://img.shields.io/static/v1?logo=Julia&logoColor=white&label=&message=Julia&color=9558B2)](https://julialang.org) !-->,
albeit
[![FORTRAN](https://img.shields.io/static/v1?logo=Fortran&label=&message=FORTRAN&color=4D41B1)](https://fortran-lang.org)
[![C/C++](https://img.shields.io/static/v1?logo=Coursera&label=&message=C/C%2B%2B&color=00599C)](https://isocpp.org)
[![Python](https://img.shields.io/static/v1?logo=Python&logoColor=white&label=&message=Python&color=3776AB)](https://www.python.org)
[![Wolfram](https://img.shields.io/static/v1?logo=Wolfram-Language&logoColor=white&label=&message=Wolfram&color=DD1100)](https://www.wolfram.com)
[![MATLAB](https://img.shields.io/static/v1?logo=MathWorks&logoColor=white&label=&message=MATLAB&color=0076A8)](https://www.mathworks.com)
may also marginally used.
As [CFD](https://en.wikipedia.org/wiki/Computational_fluid_dynamics) is a mature community we utilized various robust tools built by the community, which require user/contributors to do some preparation before setting the project up.

### Prerequisites

To use or contribute to the project, these tools should be previously installed:

- [![Julia](https://img.shields.io/static/v1?logo=Julia&logoColor=white&label=Julia&message=1.6.0+and+above&color=9558B2)](https://julialang.org)
- [![Python](https://img.shields.io/static/v1?logo=Python&logoColor=white&label=Python&message=3.8+++and+above&color=3776AB)](https://www.python.org)
- [![Conda](https://img.shields.io/static/v1?logo=Anaconda&logoColor=white&label=Conda&message=4.5+++and+above&color=44A833)](https://docs.conda.io)

### Setup

- Clone the repository with [![Git](https://img.shields.io/static/v1?logo=Git&logoColor=white&label=&message=Git&color=F05032)](https://git-scm.com):
```shell
git clone git@github.com:tsinghua-TEEP/CFD2021-G4-Projects.git
```
- Instantiate the dependencies.
  - *(to be implemented)* In the root directory of the repository, run:
    ```shell
    python ./setup.py
    ```
  - For now, you need to instantiate the dependencies yourself.
    - **Julia**: type ``julia`` in the OS shell to enter the Julia REPL, then type ``]`` to enter the Pkg interface;
      from there:
    ```jldoctest
    (@v1.6) pkg> activate .
    (CFD2021Projects) pkg> instantiate
    ```
      in case you need to update the dependencies after that, do
    ```jldoctest
    (CFD2021Projects) pkg> update; precompile
    ```
    - **Conda**: in the repo root, from the OS shell execute
    ```shell
    conda env create --prefix ./.conda/env/CFD2021-G4-Projects --file ./conda-environment.yml
    ```
    - **Git submodules**: some of the dependencies are git submodules. In the repo, from the OS shell execute
    ```shell
    git submodule update --init --recursive
    ```

## Contents

Projects are listed as below:

### Project 1

Analytical and elliptic grid generation
for geometrical and general setups.

#### TODOs:
- [ ] General purpose solver for the inverted Poisson equation. (see ``gauss-seidel.jl`` from ``CFD-Julia``)
- [ ] Processing general [NACA airfoils](https://en.wikipedia.org/wiki/NACA_airfoil)
      (formula generation, discretization, splining, etc.)
- [ ] Port the output to standard
      [![Gmsh](https://img.shields.io/static/v1?logo=Vercel&logoColor=white&label=&message=Gmsh&color=000000)](https://gmsh.info)
      format.
