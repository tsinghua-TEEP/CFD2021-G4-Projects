# The first team project of the CFD course.

对下列两种几何形状, 通过求解椭圆型方程生成贴体网格.

- [ ] 在直径为15的圆内,有一个直径为2的圆, 二者中心重合.
- [ ] NACA0012翼型(尖尾缘). 远场形状为圆形或者半圆与矩形组合, 直径是机翼弦长的5倍.

## Implementation details

对于要求的两种场景, 首先进行解析网格生成, 然后再进行椭圆网格生成.
要注意, 椭圆生成中要求解的是经反变换的 Poisson 方程;
这是一个非线性程度很高的椭圆型方程,
因此最好首先经验地得到不错的解析网格, 再使用椭圆方法.

为了实现一个通用的网格生成器, 总是认为计算平面的计算域是
$\left(\xi, \eta\right) \in \left[0, 1\right]^2$;
从而, 当我们的网格大小是 $M_x \times M_y$,
我们求解的点总是
$\xi_i = \frac{i}{M_x}, \eta_j = \frac{j}{M_y}$;
对于非周期边界, 一般使下标从 1 开始;
对于周期边界和有物理意义的边界,
通过 [OffsetArray](https://github.com/JuliaArrays/OffsetArrays.jl) 使用从 0 开始的下标
(有一定的 overhead, 但对于一份教学代码是完全可以接受的);
**注意: ``OffsetArray`` 对应「场」的几何图像; 在类似数组拼接的任务上, 隐含使用了「列表」的几何图像, 应当转换为普通数组操作, 再化为 ``OffsetArray``.**  
由于我们只使用 1 阶, 2 阶有限差分格式,
使用 [FiniteDiff.jl](https://github.com/JuliaDiff/FiniteDiff.jl)
比使用 [FiniteDifferences.jl](https://github.com/JuliaDiff/FiniteDifferences.jl)
要更合适.

### 圆环

很容易写出此时的亚纯映射:

$$\begin{aligned}&
\begin{cases}
  \xi  = \frac{\theta}{2\pi} = \frac{1}{2\pi}\arctan{\left(x,y\right)}\\
  \eta = \frac{\ln{{\;  r  \;}/R_\text{i}  }}{   \ln{    R_\text{e}    /R_\text{i}}}
       = \frac{\ln{{(x^2+y^2)}/R_\text{i}^2}}{2\;\ln{\;\;R_\text{e}\;\;/R_\text{i}}}
\end{cases} \quad
\begin{cases}
  \theta = 2\pi \xi\\
  r      = \left(\frac{R_\text{e}}{R_\text{i}}\right)^\eta\!R_\text{i}
\end{cases} \\&
\begin{cases}
  x = r \cos{\theta}
    = \left(\frac{R_\text{e}}{R_\text{i}}\right)^\eta\!R_\text{i} \cos{2\pi \xi}\\
  y = r \sin{\theta}
    = \left(\frac{R_\text{e}}{R_\text{i}}\right)^\eta\!R_\text{i} \sin{2\pi \xi}
\end{cases}
\end{aligned}$$

### NACA0012

如果生成 O 型网格, 则应当选择圆形远场. 应当平滑地封闭尾缘.  
如果生成 C 型网格, 则应当选择矩形远场. 可以任意地封闭尾缘, 包括使用尖尾缘.  
此外存在翼形的对齐问题:
公开资料中多将尾缘对齐到圆心/原点, 但实际上这一选择看来是相当任意的.
我们约定在不进行实际流场计算时, 总是将前缘对齐到圆心/原点.

此外, 考虑到我们本次计算的是 00xx 对称翼形, 只对上半平面进行计算.

## References
- [Guide](http://www.people.virginia.edu/~rjr/mae672/projects/GridGeneration.pdf) from a similar course at Virginia U., intended for a FORTRAN implementation.
- P.R. Eiseman. [*Grid Generation for Fluid Mechanics Computations*](https://doi.org/10.1146/annurev.fl.17.010185.002415). Ann. Rev. Fluid Mech. 1985. 17: 487-522.
- S.P. Spekreijse. [*Elliptic Grid Generation Based on Laplace Equations and Algebraic Transformations*](https://doi.org/10.1006/jcph.1995.1078). Journal of Computational Physics 118, 38-61.
- S.P. Spekreijse, J.W. Boerstoel. [*Multiblock grid generation, Part I: Elliptic grid generation methods for structured grids*](https://core.ac.uk/download/pdf/80112194.pdf). National Aerospace Laboratory NLR. Amsterdam, Belgium, March 1996.
- M.A. Sadybekov, A.A. Dukenbayeva. [*Direct and inverse problems for the Poisson equation with equality of flows on a part of the boundary*](https://doi.org/10.1080/17476933.2018.1517340). Complex Variables and Elliptic Equations Vol. 64, 2019, 777-791.
- Toby Driscoll. [Schwarz-Christoffel mapping](http://www.math.udel.edu/~driscoll/research/conformal.html).
- SchwarzChristoffel.jl 包文档. [说明](https://jdeldre.github.io/SchwarzChristoffel.jl).
- Gmsh 4.8.3 文档. [Gmsh 文件格式](http://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format).
- tlanyan. [Gmsh 网格文件格式说明](https://tlanyan.me/gmsh%E7%BD%91%E6%A0%BC%E6%96%87%E4%BB%B6%E6%A0%BC%E5%BC%8F%E8%AF%B4%E6%98%8E/). 
