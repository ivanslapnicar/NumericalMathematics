### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# ╔═╡ cefa7a08-9c72-4629-9b80-5b9ac9167fd1
begin
	import Pkg
	Pkg.activate(mktempdir())
	Pkg.add("PlutoUI")
end

# ╔═╡ 6a60299f-8efb-4210-bd51-5efeb88e450d
using PlutoUI, Random, LinearAlgebra

# ╔═╡ d3bf091a-28ba-4e5f-b09d-e5205c5abfa7
TableOfContents(title="📚 Table of Contents", aside=true)

# ╔═╡ b523e28f-4ed0-4b24-ba02-1806e38fe81e
md"""
# Norms


__Norm__ on a vector space $X$ is any function $\| \phantom{x} \| : X\to \mathbb{R}$ with  the following properties:

1.  $\| x\|=0\| \Leftrightarrow x=0$
2.  $\| \lambda x\|=|\lambda| \|x\|$
3.  $\| x+y\| \leq \|x\|+\|y\|$ (triangle inequality)
"""

# ╔═╡ d2fdc7a1-a28d-4daa-aee7-cfda4a353e8e
md"""
## Vector norms

For $X=\mathbb{R}^n$ we have

$$\|x\|_p=\big(\sum_{i=1}^n |x_i|^p\big)^{1/p}$$

Specially:

*  $\|x\|_1=\sum_{i=1}^n |x_i|\qquad$  (Manhattan norm or Taxicab norm)
*  $\|x\|_2=\sqrt{\sum_{i=1}^n x_i^2}= \sqrt{x\cdot x}\qquad$ (Euclidean norm)
*  $\|x\|_\infty = \max\limits_{i=1,\ldots,n} |x_i|\qquad$ (Maximum norm)
"""

# ╔═╡ dfa1aa39-0ea5-472d-af9e-a3a631438a15
begin
	Random.seed!(1244)
	x=rand(-7:7,5)
end

# ╔═╡ 562147d2-3146-44cb-9d54-16331d26828c
norm(x,1), norm(x), norm(x,Inf)

# ╔═╡ 6180e573-f082-4273-87a1-0699c167daea
md"""
## Matrix norms

From every vector norm we can derive a matrix norm (__induced norms__):

$$\|A\| = \max\limits_{x\neq 0} \frac{\|Ax\|}{\|x\|}=\max\limits_{\|x\|=1} \|Ax\|$$

Specially:

*  $\|A\|_1=\max\limits_{j=1:n} \sum_{i=1}^n |a_{ij}|\qquad$  (largest 1-norm of a column)
*  $\|A\|_{\infty}=\max\limits_{i=1:n} \sum_{j=1}^n |a_{ij}|\qquad$  (largest 1-norm of a row)
*  $\|A\|_2\qquad$  (largest singular value of matrix $A$)

_Frobenius_ or _Euclidean_ norm

$$\|A\|_F =\sqrt{\sum_{i,j=1}^n a_{ij}^2}$$

is not an induced norm.

Matrix norms also have the property

$$
\|A\cdot B\|\leq \|A\| \cdot \| B\|.$$
"""

# ╔═╡ 220d4a5d-4cbe-467a-b345-12a57f82aa12
A=rand(-4:4,5,5)

# ╔═╡ cc19cc86-b88a-486f-a5de-b74d2252999b
norm(A,1), norm(A), norm(A,2), norm(A,Inf), opnorm(A),maximum(svdvals(A))

# ╔═╡ e705726b-842f-4079-a5d9-079872749470
md"""
## Scalar (dot)  product, norm and orthogonality 

__Scalar product__ or __dot product__ on a vector space $X$ is every map 
$\cdot : X\times X \to \mathbb{R}$ with the following properties:

1.  $x\cdot x\geq 0$
1.  $x\cdot x=0 \Leftrightarrow x=0$
2.  $x\cdot y=y\cdot x$
3.  $(\alpha x)\cdot y =\alpha (x\cdot y)$
3.  $(x+y)\cdot z=x\cdot z+y \cdot z$

If scalar product is defined on a vector space, we can define norm as

$$\|x\|=\sqrt{x\cdot x}.$$

Also, if $x \cdot y=0$ we say that the vectors $x$ and $y$ are __mutually orthogonal (perpendicular)__.  

For example, the standard vector norm 

$$\|x\|_2=\sqrt{\sum_{i=1}^n x_i^2}= \sqrt{x\cdot x}$$

is defined by the dot product of vectors, 

$$x\cdot y=\sum_{i=1}^n  x_i y_i,$$

and vectors $x$ and $y$ are orthogonal, $x\perp y$, if
$x\cdot y=0$.

The scalar product of functions is defined via (definite) integral:

$$f\cdot g = \int_a^b f(x)g(x) \, dx.$$

The other definitions remain the same:

$$\| f\|_2= \sqrt{f\cdot f} = \sqrt{\int_a^b [f(x)]^2 \, dx},$$

$$f\perp g \Longleftrightarrow f\cdot g =0.$$
"""

# ╔═╡ Cell order:
# ╠═cefa7a08-9c72-4629-9b80-5b9ac9167fd1
# ╠═6a60299f-8efb-4210-bd51-5efeb88e450d
# ╠═d3bf091a-28ba-4e5f-b09d-e5205c5abfa7
# ╟─b523e28f-4ed0-4b24-ba02-1806e38fe81e
# ╟─d2fdc7a1-a28d-4daa-aee7-cfda4a353e8e
# ╠═dfa1aa39-0ea5-472d-af9e-a3a631438a15
# ╠═562147d2-3146-44cb-9d54-16331d26828c
# ╟─6180e573-f082-4273-87a1-0699c167daea
# ╠═220d4a5d-4cbe-467a-b345-12a57f82aa12
# ╠═cc19cc86-b88a-486f-a5de-b74d2252999b
# ╟─e705726b-842f-4079-a5d9-079872749470
