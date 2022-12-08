SpecSolve
====================

A MATLAB implementation for computing spectral measures of self-adjoint operators [1] written by Matt Colbrook and Andrew Horning. The code supports

(1) Differential operators on the real line, with variable coefficients,

(2) Integral operators on [-1,1], with smooth kernels,

(3) Infinite matrices, with finite bandwidth or rapid off-diagonal decay [2].

(4) General function handles for computing the resolvent and inner products.

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
Computes a smoothed approximation to the spectral measure of radial Schrodinger operators. See `Example_radialSchrodinger.m` for a worked example.

genMeas()
===================
Computes a smoothed approximation to the spectral measure given function handles for the resolvent and inner products.

References
===================

[1] M. J. Colbrook, A. Horning, and A. Townsend. "Computing spectral measures of self-adjoint operators." SIAM Review 63.3 (2021): 489-524.

[2] M. J. Colbrook. "Computing spectral measures and spectral types." Communications in Mathematical Physics 384.1 (2021): 433-501.

[3] M. J. Colbrook, A. Horning. "SpecSolve: Spectral methods for spectral measures." ICOSAHOM (2022).

