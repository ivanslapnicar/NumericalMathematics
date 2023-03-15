### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

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
# Scalar (dot)  product, norm and orthogonality

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

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[compat]
PlutoUI = "~0.7.9"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.1+0"

[[Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "8076680b162ada2a031f707ac7b4953e30667a37"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.2"

[[Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.20+0"

[[Parsers]]
deps = ["Dates"]
git-tree-sha1 = "438d35d2d95ae2c5e8780b330592b6de8494e779"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.0.3"

[[PlutoUI]]
deps = ["Base64", "Dates", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "Suppressor"]
git-tree-sha1 = "44e225d5837e2a2345e69a1d1e01ac2443ff9fcb"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.9"

[[Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[Suppressor]]
git-tree-sha1 = "a819d77f31f83e5792a76081eee1ea6342ab8787"
uuid = "fd094767-a336-5f1f-9728-57cf17d0bbfb"
version = "0.2.0"

[[Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.1.1+0"
"""

# ╔═╡ Cell order:
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
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
