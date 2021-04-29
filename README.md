# CFD2021-G4-Projects

<!-- [![CMake](https://github.com/tsinghua-TEEP/CFD2021-G4-Projects/actions/workflows/cmake.yml/badge.svg)](https://github.com/tsinghua-TEEP/CFD2021-G4-Projects/actions/workflows/cmake.yml) !-->
<!-- [![Python package with Conda](https://github.com/tsinghua-TEEP/CFD2021-G4-Projects/actions/workflows/python-package-conda.yml/badge.svg)](https://github.com/tsinghua-TEEP/CFD2021-G4-Projects/actions/workflows/python-package-conda.yml) !-->
[![codacy](https://app.codacy.com/project/badge/Grade/8ddf95075915482d8708388554f16386?label=)](https://www.codacy.com?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=tsinghua-TEEP/CFD2021-G4-Projects&amp;utm_campaign=Badge_Grade)
[![circleci](https://circleci.com/gh/tsinghua-TEEP/CFD2021-G4-Projects.svg?style=shield&label=CircleCI&logo=CircleCI&circle-token=9b51e15e5ced695a347386f06bdc605e23e7d8e5)](https://github.com/tsinghua-TEEP/CFD2021-G4-Projects/actions)
[![codecov](https://codecov.io/gh/tsinghua-TEEP/CFD2021-G4-Projects/branch/main/graph/badge.svg?token=9R7SWYU9W5&logoColor=white&style=flat)](https://codecov.io/gh/tsinghua-TEEP/CFD2021-G4-Projects)
[![License: Apache-2.0](https://img.shields.io/badge/License-Apache%202.0-blue.svg?logo=Apache&style=flat)](https://github.com/tsinghua-TEEP/CFD2021-G4-Projects/blob/main/LICENSE)
[![wakatime](https://wakatime.com/badge/github/tsinghua-TEEP/CFD2021-G4-Projects.svg?logo=WakaTime&style=flat)](https://wakatime.com/badge/github/tsinghua-TEEP/CFD2021-G4-Projects)  
[![ColPrac](https://img.shields.io/badge/ColPrac-Contributor's%20Guide-blueviolet?logo=Julia&logoColor=white&style=flat)](https://github.com/SciML/ColPrac)
[![Discussions](https://img.shields.io/badge/GitHub-Discussions-blueviolet?logo=github&style=flat)](https://github.com/tsinghua-TEEP/CFD2021-G4-Projects/discussions)
<!-- [![Dependabot Status](https://api.dependabot.com/badges/status?host=github&repo=tsinghua-TEEP/CFD2021-G4-Projects)](https://dependabot.com) !-->

This repository dedicates to the team projects
related to the CFD course, 2021 spring, at SAS, THU.  
Efforts are contributed by a joined 4-member group.  
For practice and performance, we attempt to use the
[Julia](https://julialang.org) programming language in these projects.

## Configuration

The project is mainly written in [Julia](https://julialang.org), albeit FORTRAN, C/C++, Python, Wolfram, MATLAB may also marginally used. As CFD is a mature community we utilized various robust tools built by the community, which require user/contributors to do some preparation before setting the project up.

### Prerequisites

To use or contribute to the project, these tools should be previously installed:

- [Julia](https://julialang.org)   (1.6.0 and above)
- [Python](https://www.python.org) (3.8   and above)
- [Conda](https://docs.conda.io)   (4.5   and above)

### Setup

- Clone the repository with git:
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
- [ ] Port the output to standard [Gmsh](http://gmsh.info) format.
