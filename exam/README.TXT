Adaptive 1D Integrator with Random Nodes
=========================================

DESCRIPTION
-----------
This project implements an adaptive one-dimensional integrator where the
quadrature abscissas at each recursion level are chosen randomly (stratified
sampling), rather than at fixed predetermined points.

The algorithm works as follows:
  - At each subinterval, N points are drawn uniformly at random and their
    mean and variance are computed.
  - The error is estimated as |mean_new - mean_inherited| * (b-a), reusing
    the mean passed down from the parent interval. On the very first call
    (no parent), the sample variance is used instead: 3*sqrt(Var/N)*(b-a).
  - If the error exceeds the tolerance, the interval is halved and the
    algorithm recurses, halving the absolute tolerance by sqrt(2) at each
    level to keep the global error budget consistent.
  - A depth cap of 52 prevents infinite recursion near true singularities.

Files:
  quad_random.cc / quad_random.h
      Core implementation: integrate_random (plain random abscissas) and
      integrate_random_cc (with Clenshaw-Curtis variable transformation).

  sfuns.cc / sfuns.h
      Predefined-node adaptive integrators (open 4-point rule and
      Clenshaw-Curtis), used as the baseline for comparison.

  main.cc
      Drives all experiments and writes data files to data/.

  comparison.gnu / sweep.gnu
      Gnuplot scripts that produce SVG plots in plots/.

CLENSHAW-CURTIS TRANSFORMATION
-------------------------------
Plain random abscissas converge poorly (or not at all) on integrals with
endpoint singularities such as 1/sqrt(x) and ln(x)/sqrt(x) on [0,1].
integrate_random_cc applies the substitution x = cos(pi*t), which maps
[0,1] in t to [0,1] in x and folds the Jacobian sin(pi*t) into the
integrand, damping the singularity. The t-interval is opened by eps_t=1e-9
on each side to avoid evaluating the integrand exactly at x=0.

OPTIMAL N
---------
The sweep experiment (write_sweep in main.cc, plotted by sweep.gnu) runs
integrate_random_cc for N in {2, 4, 8, 16, 32, 64, 128, 256} on three
test integrals, averaging over 8 runs each (acc = eps = 1e-4):

  - integral of sqrt(x)        on [0,1]  = 2/3
  - integral of 1/sqrt(x)      on [0,1]  = 2
  - integral of ln(x)/sqrt(x)  on [0,1]  = -4

Result: N=4 to N=8 gives the best trade-off between function evaluations
and accuracy. Small N leads to many recursive splits (high ncalls); large N
wastes evaluations per level without proportional error reduction. N=4 is
used as the default throughout.

COMPARISON WITH PREDEFINED-NODE INTEGRATORS
--------------------------------------------
Four methods are compared on four integrals (including sqrt(1-x^2) = pi/4):

  fixed         open 4-point adaptive rule (sfuns)
  fixed_cc      Clenshaw-Curtis adaptive rule (sfuns)
  random_N*     adaptive random integrator, N points per level
  random_cc_N*  adaptive random + Clenshaw-Curtis transform

The comparison data (ncalls, actual error, estimated error, value) is
written to data/<integral>.data and plotted by comparison.gnu.

Key findings:
  - For smooth integrands, the fixed-node integrators are more efficient
    (fewer function evaluations for the same accuracy).
  - For singular integrands, random+CC is competitive and more robust,
    reliably hitting the tolerance where plain random abscissas struggle.
  - The fixed CC rule is generally the most efficient overall; random+CC
    is a close second and handles edge cases (e.g. depth-capped intervals)
    more gracefully via its non-finite sample guard.

BUILD AND RUN
-------------
  make        — compiles, runs, and produces all plots
  make run    — compile and run only (no plots)
  make clean  — deletes the generated data and plots

Requires: g++ with C++23 support, gnuplot (for plots only).

SELF-EVALUATION
---------------
7 / 10

The core algorithm, point reuse, CC transform, N sweep, and comparison
are all implemented and working. The main limitation is that the random
integrator is inherently slower than the fixed-node rules on smooth
functions, which is expected but worth noting. The depth cap and non-finite
guard add robustness but are somewhat ad hoc. A more principled approach
(e.g. importance sampling or a better variance estimator) could improve
convergence on the hardest singularities.
