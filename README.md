# Fast Algorithms for Optimal Damping in Mechanical Systems

MATLAB implementations of the **BB-step projected residual algorithm (BBPR)** and the **spectral projected gradient method (SPG)** described by Qingna Li and Françoise Tisseur. The algorithms compute nonnegative coefficients for viscous dampers in a mechanical system. The objective is the trace of the solution of a Lyapunov equation associated with the selected vibration modes.

The driver `do_table2.m` compares BBPR and SPG in the experiments of Section 5.1 (Table 2). The Julia comparison of SPG and FODA in Section 5.2 is a separate experiment.

## Running the Table 2 experiments

Open MATLAB and set the current folder to the directory containing the seven `.m` files. The driver generates the test matrices, including the beam examples. The code uses MATLAB's `normalize` function and implicit array expansion.

To run all nine Table 2 configurations with both BBPR and SPG, enter the following command in MATLAB:

```matlab
do_table2
```

To run a single configuration, change the loop header in `do_table2.m` from `for pb = (1:9)` to the desired case index. For example, use `for pb = 1` for `damp1-a` (order 4).

| Case | Problem | System order | Dampers | Selected modes | Initial coefficients |
|---:|---|---:|---:|---:|---|
| 1 | damp1-a | 4 | 1 | 4 | e |
| 2 | damp1-b | 20 | 1 | 20 | e |
| 3 | damp1-c | 20 | 2 | 20 | e |
| 4 | damp1-c | 20 | 2 | 20 | 10e |
| 5 | damp2-a | 801 | 3 | 27 | e |
| 6 | damp2-b | 1601 | 3 | 27 | e |
| 7 | damp2-c | 2001 | 3 | 20 | e |
| 8 | beam-a | 200 | 3 | 40 | e |
| 9 | beam-b | 1000 | 5 | 80 | e |

Here e is the vector of ones. Cases 3 and 4 use the same system and damper locations with different initial coefficients.

## Results

The driver prints the iteration count (`#iter`), the number of SPG iterations requiring backtracking (`#ls`), the number of complete eigendecompositions during optimization (`#eig`), and the final objective value when available. A single SPG iteration can require several backtracking steps. Each additional trial point contributes to `#eig`. The initial modal decomposition of the mass and stiffness matrices is common to both methods and is excluded from this count.

The driver sets the tolerance for the projected residual to `1e-8` and the tolerance for the relative change in the damping vector to `1e-5`. Both tests must hold. Each method permits at most 1000 iterations.

After a run, the damping vectors and solver information for the last case are available in the MATLAB workspace:

```matlab
nu_BBPR
nu_SPG
info_BBPR.reshist     % norms of the projected residual for BBPR
info_SPG.reshist      % norms of the projected residual for SPG
info_SPG.iter_ls      % iterations requiring backtracking
```

## MATLAB files

| File | Purpose |
|---|---|
| `do_table2.m` | Constructs the test problems, computes the modal data, and calls both solvers. |
| `BBPR.m` | Implements BBPR. Call with `[nu,info] = BBPR(opts)`. |
| `SPG.m` | Implements SPG with a nonmonotone line search. Call with `[nu,info] = SPG(opts)`. |
| `EA_eig.m` | Computes the eigendecomposition for the quadratic eigenvalue problem with diagonal internal damping and low-rank external damping in modal coordinates. |
| `SEA.m` | Implements the structured Ehrlich–Aberth iteration used by `EA_eig.m`. |
| `apply_lyapunov_correction.m` | Applies a Lyapunov correction to a block of vectors. |
| `damped_beam.m` | Generates the finite element beam matrices used to construct cases 8 and 9. |

## Reference

Qingna Li and Françoise Tisseur, *Fast Algorithms for Optimal Damping in Mechanical Systems*. [arXiv:2601.05404](https://arxiv.org/abs/2601.05404).
