# implicit-func-relaxations

Suppose an *implicit function* $\mathbf{x}$ is defined in terms of a
known *residual function* $\mathbf{f}$ as satisfying:

$$
\mathbf{f}(\mathbf{x}(\mathbf{p}),\mathbf{p}) = \mathbf{0},
\qquad\forall\mathbf{p}\in P.
$$

Suppose we can also generate convex/concave relaxations $\mathbf{f}^{\text{cv}}/\mathbf{f}^{\text{cc}}$ for $\mathbf{f}$,
perhaps in the sense of
[McCormick.jl](https://github.com/PSORLab/McCormick.jl), and that we
know a set $X$ containing the range of $\mathbf{x}$. This
repository illustrates our new approach for constructing
convex/concave relaxations for the unknown implicit function
$\mathbf{x}$ in terms of known information, and contains our Julia code for all numerical examples in the
corresponding manuscript.

This implementation was developed by Huiyi Cao in Julia. This repository is
tied to the accompanying manuscript, and will not be updated except for bug
fixes. If you make use of this code, please cite the accompanying manuscript.

This work was supported by the McMaster Advanced Control Consortium
(MACC), and by the Natural Sciences and Engineering Research Council of Canada (NSERC) under Grant RGPIN-2017-05944.

## Dependencies

- McCormick.jl v0.13.1
- DataFrames.jl v1.3.4
- IntervalArithmetic.jl v0.20.6
- JuMP.jl v0.21.5
- Ipopt.jl v0.6.5
- NLsolve.jl v4.5.1
- Plots.jl v1.31.3
- LaTeXStrings v1.3.0

## Method outline

With the setup described above, we describe convex/concave relaxations
of the implicit function $\mathbf{x}$ as follows, for each component
$i$ and each $\mathbf{p}\in P$:

$$
\begin{align*}
x_i^{\mathrm{cv}}(\mathbf{p})&:=\inf_{\mathbf{z}\in X} z_i
\quad\text{subject
to}\quad\mathbf{f}^{\text{cv}}(\mathbf{z},\mathbf{p})\leq\mathbf{0}\leq\mathbf{f}^{\text{cc}}(\mathbf{z},\mathbf{p}),
\\
x_i^{\mathrm{cc}}(\mathbf{p})&:=\sup_{\mathbf{z}\in X} z_i
\quad\text{subject
to}\quad\mathbf{f}^{\text{cv}}(\mathbf{z},\mathbf{p})\leq\mathbf{0}\leq\mathbf{f}^{\text{cc}}(\mathbf{z},\mathbf{p}). 
\end{align*}
$$

When monotonicity can be exploited, these relaxations are evaluated
particularly easily.
This formulation does not actually require existence or uniqueness of
$\mathbf{x}$ on the presumed domain. Inverse functions and
constraint-satisfaction problems may be relaxed analogously. For more details, please refer to
the accompanying manuscript.
