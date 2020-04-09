SpecSolve
====================

A MATLAB implementation for computing spectral measures of self-adjoint operators [1] written by Matt Colbrook and Andrew Horning. The code supports

(1) Differential operators on the real line, with variable coefficients,

(2) Integral operators on [-1,1], with smooth kernels,

(3) Infinite matrices, with finite bandwidth or rapid off-diagonal decay [2].

See `Example_*.m`  files for several examples that appear in [1]. **SpecSolve** functions **diffMeas()**, **intMeas()**, and **rseMeas()** make use of the **Chebfun** software package for computing with functions, available for download at https://www.chebfun.org/.


diffMeas()
====================
Computes a smoothed approximation to the spectral measure of an ordinary, linear differential operator `L` acting on functions on the real line. `L` has the form `L = a_0 + a_1 D_1 + ... + a_p D_p` where `D_k` is the kth derivative operator and `a_k=a_k(x)` is a smooth variable coefficient.


intMeas()
====================
Computes a smoothed approximation to the spectral measure of a linear integral operator `L` acting on functions on `[-1,1]`. `L` has the form `[Lu](x) = a(x)u(x) + \int K(x,y)u(y) dy` where `a(x)` and `K(x,y)` are smooth functions on `[-1,1]` and `[-1,1]^2`, respectively.


infmatMeas()
===================
Computes a smoothed approximation to the spectral measure of a lattice operator `A` acting on square summable sequences.


rseMeas()
===================
Computes a smoothed approximation to the spectral measure of radial Schrodinger operators. See **radialSchrodinger_example.m** for a worked example.

References
===================

[1] M. J. Colbrook, A. Horning, and A. Townsend. "Computing spectral measures of self-adjoint operators." arXiv preprint. 

[2] M. J. Colbrook. "Computing spectral measures and spectral types." arXiv preprint arXiv:1908.06721v2, 2019.

