### A Pluto.jl notebook ###
# v0.12.4

using Markdown
using InteractiveUtils

# ╔═╡ d615baaa-4c22-4df5-94b4-832f2a29063f
md"""
# Perturbation Theory, Backward Error and  Stable Algorithms

## Perturbation theory

__Question:__

> _How much does the result change with respect to the change in input data?_

For some function $f(x)$ and input datum $x$, we want to estimate
__absolute error__ w.r.t.  __change of__ input datum for $\delta x$,

$$
\| f(x+\delta x)-f(x)\|\leq \kappa \|\delta x\|,$$

and __relative error__ w.r.t. __relative change of__ input datum for $\displaystyle\frac{\| \delta x\|}{\|x\|}$,

$$
\frac{\| f(x+\delta x)-f(x)\|}{\| f(x) \|}\leq \kappa_R \frac{\| \delta x \|}{ \|x\|}.$$


We have 

$$
\| f(x+\delta x)-f(x)\| = \frac{\| f(x+\delta x)-f(x)\|}{\| \delta x \|} \|\delta x\| \equiv \kappa \|\delta x\|.$$

The quantity $\kappa$ is the __condition number__ or __condition__. 
It reminds us of a formula for derivative, and it tells us how much is the perurbation of input data magnified, at most.

Similarly, in the expression

$$
\frac{\| f(x+\delta x)-f(x)\|}{\| f(x) \|}= \frac{\| f(x+\delta x)-f(x)\|\cdot  \|x\| }{\|\delta x\| \cdot\| f(x)\|}
\cdot \frac{\|\delta x\|}{\|x\|} \equiv \kappa_R \frac{\|\delta x\|}{\|x\|},$$

so $\kappa_R$ tells us how much is the relative perurbation of input data relatively magnified, at most.
"""

# ╔═╡ 77a77817-06b7-4cfe-8326-840a8ef11997
md"""
## Backward error

Let the value $f(x)$ be computed by some algorithm $\mathrm{alg(x)}$.

The __algorithm error__ is 

$$
\|\mathrm{alg(x)}-f(x)\|,$$

and the __relative algorithm error__ is 

$$
\frac{\| \mathrm{alg}(x)-f(x)\|}{\| f(x) \|}.$$

These error are, in general, hard or impossible to estimate directly, 
and we study the  __backward error__,

$$
\mathrm{alg}(x)=f(x+\delta x),$$

instead. In other words, 

> the computed value of a function $f$ for input $x$ is equal to the exact value of $f$ for the input pertubed with some (unknown) perturbation.
"""

# ╔═╡ daef9a23-10bc-4a2e-af92-dccb1657b1aa
md"""
## Stable algorithms

An algoritam is __stable__ if for every input $x$ 

$$
\mathrm{alg}(x)=f(x+\delta x)$$

for some small $\delta x$.
"""

# ╔═╡ 56a6197c-c597-484c-aafd-f0e247cdfe90


# ╔═╡ Cell order:
# ╟─d615baaa-4c22-4df5-94b4-832f2a29063f
# ╟─77a77817-06b7-4cfe-8326-840a8ef11997
# ╟─daef9a23-10bc-4a2e-af92-dccb1657b1aa
# ╠═56a6197c-c597-484c-aafd-f0e247cdfe90
