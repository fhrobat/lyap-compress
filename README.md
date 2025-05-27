# LYAP-COMPRESS

Contains the MATLAB code for the memory-efficient computation of the solution of the Lyapunov equation

$AX + XA = cc^T$
	
where $A$ is a symmetric positive definite real matrix and $c$ is a real vector. The algorithm is presented in [1].

The main function is ``` lyap_compress/lyap_compress.m ```, with additional helper functions contained in the utils folder.

This repository additonally contains the MATLAB code to reproduce the numerical experiments in [1]:

```test_4DLap.m```: first experiment, 4D Laplacian

```test_Fenics_Rail.m```: Model order reduction: example 1 (FEniCS_Rail)

```test_MOR.m```: Model order reduction: example 2

Some functions required to run the two-pass Lanczos method that we compare against our method are contained in the two_pass_lanc folder; the two-pass method rely on helper functions contained in the utils folder.

[1] Angelo A. Casulli, Francesco Hrobat, Daniel Kressner, Lanczos with compression for symmetric matrix Lyapunov equations.

## Dependencies
There are no external dependencies for the main function ``` lyap_compress/lyap_compress.m ```. The test scripts have the following dependencies, required for the comparison against the restart method:

* https://gitlab.com/katlund/compress-and-restart-KSM
