# MATH-453: Computational Linear Algebra
## Alessandro Tavazzi, Mariia Soroka

This is our project for Computational Linear Algebra course at EPFL. 
We have implemented different randomized iterative linear solvers described in [[1]](#1) and [[2]](#2).

### About the project
Given a real matrix $A \in \mathbb{R}^{m×n}$ and a real vector $b \in \mathbb{R}^m$ the objective is to find
a solution of the linear system $Ax = b$.

Here we assume that the system is consistent, i.e. there is such $x^{\*}$ that $Ax^{\*} = b$.

Common approaches to tackling this problem fall into two categories: direct methods and iterative methods. 
Direct methods can be prohibitively
memory-consuming for large systems and that is the reason why iterative solvers are usually used
for large scale problems. 

The goal of this project was to study so called randomized
iterative methods, to rederive main convergence results and to perform some
numerical experiments.

### Reproducing results
To run the code and to reproduce our results it is necessary to install [Lightspeed Toolbox](https://github.com/tminka/lightspeed) for Matlab 
to get access to $\texttt{flops}$ function that is missing in current version of Matlab.



### References
<a id="1">[1]</a> 
Thomas Strohmer and Roman Vershynin. A randomized Kaczmarz algorithm with exponential convergence. J. Fourier Anal. Appl., 15(2):262--278, 2009.

<a id="1">[2]</a> 
Robert M. Gower and Peter Richtárik. Randomized iterative methods for linear systems. SIAM J. Matrix Anal. Appl., 36(4):1660--1690, 2015
