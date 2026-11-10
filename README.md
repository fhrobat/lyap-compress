# LYAP-COMPRESS

Contains the MATLAB code for the memory-efficient computation of the solution of the Lyapunov equation

<p align="center"> $AX + XA = cc^T$ </p>
	
where $A$ is a symmetric positive definite real matrix and $c$ is a real vector. The algorithm (Lanczos with compression, `compress`) is presented in [1].

The main function is ``` lyap_compress/lyap_compress.m ```, with additional helper functions contained in the ```utils``` folder. Minimal usage:

```matlab
setup_paths;
options.tol = 1e-6;       % relative tolerance on the residual norm
options.maxmem = 120;     % maximum number of vectors of length N kept in memory
% optional: options.eigmin, options.eigmax (estimates of the extreme eigenvalues of A)
[U, Y, result] = lyap_compress(A, c, options);   % X is approximated by U*Y*U'
```

If `options.eigmin` and `options.eigmax` are not provided, they are estimated from the first `maxmem - 1` Lanczos iterations (performed with full orthogonalization).

## Numerical experiments

This repository additionally contains the MATLAB code to reproduce the numerical experiments in [1]:

* ```test_2DLap.m```: 2D Laplacian (Section 5.1)
* ```test_2DLap_sensitivity.m```: sensitivity to the approximation of the smallest eigenvalue (Section 5.1, sensitivity figures)
* ```test_3DLap.m```: 3D Laplacian (Section 5.2)
* ```test_fenics_rail.m```: FEniCS Rail example (heat diffusion in 2D, Section 5.3)
* ```test_heat.m```: heat diffusion in a steel bracket (heat diffusion in 3D, Section 5.4)

The script ```run_all.sh``` runs all experiments sequentially, one MATLAB process per experiment and size:

```
MATLAB=matlab NTHREADS=4 MEMLIMIT=12G ./run_all.sh
```

(```MEMLIMIT``` is optional and requires ```systemd-run```; a run exceeding it is stopped and the next one starts.) The results are saved in ```results/<machine>_<date>/```:

* ```summary.md```: the tables of the paper (compress, two-pass, restart), the results of all the methods and the exit status of every run;
* ```<test>.csv```: the results of every experiment (all sizes);
* ```raw/```: the results saved by the test scripts after every method;
* ```logs/```, ```jobs.csv```, ```machine.txt```: MATLAB output, exit status and duration of every run, machine and MATLAB version.

A single size of an experiment can be run from MATLAB with, e.g., ```setup_paths; size_idx = 2; test_2DLap```.
For two-pass, the reported number of matrix-vector products includes both passes.

The methods we compare against (two-pass Lanczos, rational Krylov methods with recycled/adaptive poles and with multishift CG) are contained in the ```comparison``` folder; these methods rely on helper functions contained in the ```utils``` folder.

[1] Angelo A. Casulli, Francesco Hrobat, Daniel Kressner, Lanczos with compression for symmetric matrix Lyapunov equations.

## Dependencies

* ```lyap_compress/lyap_compress.m``` and the methods used for comparison require the Control System Toolbox (function ```lyap```). If the extreme eigenvalues satisfy $\lambda_{\max}/\lambda_{\min} > 10^7$, the Zolotarev poles are computed with asymptotic expansions, which require the Symbolic Math Toolbox.
* ```test_heat.m``` requires the Partial Differential Equation Toolbox.
* The methods ```restart``` and ```extended``` (```test_2DLap.m```, ```test_3DLap.m```) require the code of https://gitlab.com/katlund/compress-and-restart-KSM, to be placed in ```external/compress-and-restart-KSM```:
  ```
  git clone https://gitlab.com/katlund/compress-and-restart-KSM.git external/compress-and-restart-KSM
  ```
* ```test_fenics_rail.m``` requires the FEniCS Rail matrices, available at https://zenodo.org/records/5113560. Download and unzip ```fenics_rail_79841.zip```, ```fenics_rail_317377.zip```, ```fenics_rail_1265537.zip``` into a folder ```<src>``` and run
  ```matlab
  setup_paths; build_rail_matrices('<src>', 'matrices');
  ```
  which creates ```matrices/rail_N.mat```.
