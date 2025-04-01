# GAP code accompanying the paper “Asymmetric graphs with quantum symmetry”

This repository contains GAP code to reproduce the computations from [DSR23+].
There are 5 main scripts and 3 files containing helper functions.
The main scripts are:

  * `A7_computation.g`: verify that the solution group $\Gamma(M_{A_7}, 0)$ is isomorphic to the triple cover of $A_7$ (equivalently, the unique perfect group of order $7560$).
  * `A8_computation.g`: verify that the solution group $\Gamma(M_{A_8}, 0)$ is isomorphic to $A_8$.
  * `Mathieu.g`: verify that, if $H$ is one of the five sporadic Mathieu groups $M_{11}$, $M_{12}$, $M_{22}$, $M_{23}$ or $M_{24}$, then $\Gamma(M_H, 0)$ is perfect.
  * `PSL2q.g`: verify that, if $H = \textrm{PSL}(2,q)$ with $q$ odd and $19 \leq q \leq 349$, then $\Gamma(M_H, 0)$ is perfect.
  * `small_perfect_groups.g`: find all perfect groups $H$ with $|H| \leq 10080$ such that $\Gamma(M_H, 0)$ is perfect. The code can easily be modified to extend the search range.

The helper files are:

  * `build_system_from_group.g`: a GAP function that takes as input a finite group $H$ and returns as output the matrix $M_H$ (as defined in [DSR23+, Definition 6.1]).
  * `build_sparse_F2_system_from_group.g`: a GAP function that takes as input a finite group $H$ and returns as output the matrix $M_H$, but this time using the [`SparseMatrix`](https://docs.gap-system.org/pkg/gauss/doc/chap3.html) data structure from the GAP package [`gauss`](https://homalg-project.github.io/pkg/Gauss). This script is more suitable when $H$ is large, but it requires the extra package `gauss` and the output of this function is not compatible with `solution_group.g`.
  * `solution_group.g`: a GAP function that takes as input a matrix $M \in \mathbb{F}_2^{m \times n}$ and a vector $b \in \mathbb{F}_2^m$ and returns as output the solution group $\Gamma(M,b)$ of the linear system $Mx = b$ (as a finitely presented group).


## References

[DRS23+].  Josse van Dobben de Bruyn, David E. Roberson, and Simon Schmidt.
“Asymmetric graphs with quantum symmetry”.
Preprint, 2023.  https://arxiv.org/abs/2311.04889
