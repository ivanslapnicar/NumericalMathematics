### A Pluto.jl notebook ###
# v0.19.20

using Markdown
using InteractiveUtils

# ╔═╡ c83e3f47-b1ea-443c-a235-457ed28bcb18
using PlutoUI, Random, LinearAlgebra

# ╔═╡ 5f87bb8e-34f9-4734-841b-066e906a28b3
TableOfContents(title="📚 Table of Contents", aside=true)

# ╔═╡ 1a406352-1739-4709-85bc-6ca3ecb19253
md"""
# Iterative Methods


For large systems of linear equations, and in particular for sparse systems (few non-zero elements), and if the matrix is _strictly diagonally dominant_ , the solution can be computed using __iterative methods__
(see [Numerička matematika, section 3.8](http://www.mathos.unios.hr/pim/Materijali/Num.pdf)):

## Contraction

__Definition.__ Function $F:\mathbb{R}^n\to \mathbb{R}^n$ is a __contraction__ if there is $q<1$ such that

$$\| F(x)-F(y)\| < q\|x-y\|\qquad \forall x,y.$$

## Fixed point theorem

__Banach Fixed-Point Theorem.__
If $F$ is a contraction, then the sequence defined by

$$x_{k+1}=F(x_k)$$

converges towards the unique vector $\tilde x$ for which

$$\tilde x = F(\tilde x).$$

Vector $\tilde x$ is the  __fixed point__ of the function $F$. Error in the  $k$-th step satisfies the bounds

$$\|x_k- \tilde x\| \leq \frac{q}{1-q} \|x_k-x_{k-1}\|$$

and

$$\|x_k- \tilde x\| \leq \frac{q^k}{1-q} \|x_1-x_{0}\|,$$

where the second bound is better. The speed of convergence is __linear__ ,

$$\|x_{k+1}-\tilde x\| \leq q\| x_k-\tilde x\|.$$

"""

# ╔═╡ 610ef7a4-f0a6-42c8-a2cc-1a03cb155a22
md"""

# Jacobi and  Gauss-Seidel methods

Let

$$F(x)=Bx+c,$$

where $B$ is a square matrix. Then

$$\| F(x)-F(y)\|=\| Bx+c-(By+c)\|=\|B(x-y)\| \leq \|B\| \|x-y\|,$$

so $F$ is a contraction if

$$\|B\|=q<1.$$

Consider the system  $Ax=b$. Decompose $A$ as

$$A=D\,(L+I+U)$$

where $D$ is a diagonal matrix, $L$ is strictly lower triangular matrix and $U$ is strictly upper triangular matrix.

## Jacobi method

Let

$$B=-(L+U), \quad c=D^{-1}b.$$


If the matrix $A$ is __strictly diagonally dominant__ ,

$$\| B\|_{\infty} = \max_i \sum_{{j=1} \atop {j\neq i}}^n \frac{|a_{ij}|}{|a_{ii}|}<1,$$

then the map $F$ is a contraction (here it is possible to take other norma, as well) so the sequence

$$x_{k+1}=-(L+U)x_k+c$$

converges towards the solution of the system $x$.

## Gauss-Seidel method

Let

$$
B=-(I+L)^{-1}U, \quad c=(I+L)^{-1}\, D^{-1}b.$$

Without proof we state the following: if the matrix $A$ is strictly diagonally dominant, then the map $F$ is a contraction, and the sequence

$$
x_{k+1}=-(I+L)^{-1}Ux_k+(I+L)^{-1}D^{-1}b,$$

or

$$
x_{k+1}=-Lx_{k+1}-Ux_k+D^{-1}b,$$

converges towards the solution $x$.
"""

# ╔═╡ 1ea94870-1389-11eb-1619-47636ac1d230
md"
Let us see the factorization $A=D(L+I+U)$:
"

# ╔═╡ 583e17b2-e189-49ae-8c80-0014d53c40c2
begin
	Random.seed!(123)
	n=8
	A=rand(n,n)
	# Make the matrix diagonally dominant
	A=A+n*I
	b=rand(n)
end

# ╔═╡ 4471d452-1389-11eb-39dc-67039a7b60e9
A

# ╔═╡ 47576810-1389-11eb-3e05-97dc3c3af6c8
D=Diagonal(A)

# ╔═╡ 4c917320-1389-11eb-2c06-9dd64ebf12cd
inv(D)*A

# ╔═╡ 52b234b0-1389-11eb-3d14-717752260ec5
L=inv(D)*tril(A,-1)

# ╔═╡ 5abe2470-1389-11eb-1a3e-43f6103379d2
U=inv(D)*triu(A,1)

# ╔═╡ 5aba7e24-0424-45f2-9716-3b32a71fc610
function Jacobi(A::Array,b::Array,x::Array)
    D=Diagonal(A)
    L=inv(D)*tril(A,-1)
    U=inv(D)*triu(A,1)
    tol=1000*eps()
    d=1.0
    B=-(L+U)
    c=inv(D)*b
    q=norm(B,Inf)
    # @show q
    while d>tol
        y=B*x+c
        d=norm(x-y,Inf)
        # @show d
        x=y
    end
    x,d
end

# ╔═╡ 213d2b7b-b742-4274-9bb0-e029aec6f892
# Starting vector
x₀=rand(n)

# ╔═╡ 91b52c67-de20-4bfc-9da6-3e04ed73b990
# x is the solution, d is norm of the difference of two final iterations
x,d=Jacobi(A,b,x₀)

# ╔═╡ 0d07f057-9012-42ad-bf52-31a0f14614df
# Residual
r=A*x-b

# ╔═╡ 2cb0db8c-11e6-49e3-baf7-9fdda352a26a
# Norm of the relative residual
norm(r)/(norm(A)*norm(x))

# ╔═╡ adbd72cb-4dcf-490b-bbcb-1d681358c455
function GaussSeidel(A::Array,b::Array,x::Array)
    D=Diagonal(A)
    L=inv(D)*tril(A,-1)
    U=inv(D)*triu(A,1)
    tol=1000*eps()
    d=1.0
    # B=-inv(I+L)*U
    B=-(I+L)\U
    c=(I+L)\(inv(D)*b)
    # @show norm(U,Inf)
    y=Vector{Float64}(undef,n)
    while d>tol
        y=B*x+c
        d=norm(x-y)
        x=y
    end
    x,d
end

# ╔═╡ 9ae4a166-f3ee-4850-89ff-4c0a41bae48c
xᵧ,dᵧ=GaussSeidel(A,b,x₀)

# ╔═╡ 4a832be5-26c0-46fa-bcbb-3afb4eb0cf2e
# Residual
A*xᵧ-b

# ╔═╡ 92d8992a-f55e-4576-a1aa-c8b0de9e5806
md"""
Let us measure the speed for larger matrices:
"""

# ╔═╡ 54e5e9c4-5ad2-45f9-a3ed-0aaced61663d
begin
	n₁=1024
	A₁=rand(n₁,n₁)+n₁*I
	b₁=rand(n₁)
	x₁=rand(n₁)
end

# ╔═╡ 8794b613-ddab-4fa9-82e6-eec4192705dd
@time GaussSeidel(A₁,b₁,x₁);

# ╔═╡ 02d632a3-3ab8-4b02-ba94-709880df6313
@time A\b;

# ╔═╡ aba2f7f6-8690-4ede-8282-ad925e4aae8d
md"
__Problem.__ Try to rewrite our functions to allocate less memory.
"

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
# ╠═c83e3f47-b1ea-443c-a235-457ed28bcb18
# ╠═5f87bb8e-34f9-4734-841b-066e906a28b3
# ╟─1a406352-1739-4709-85bc-6ca3ecb19253
# ╟─610ef7a4-f0a6-42c8-a2cc-1a03cb155a22
# ╟─1ea94870-1389-11eb-1619-47636ac1d230
# ╠═583e17b2-e189-49ae-8c80-0014d53c40c2
# ╠═4471d452-1389-11eb-39dc-67039a7b60e9
# ╠═47576810-1389-11eb-3e05-97dc3c3af6c8
# ╠═4c917320-1389-11eb-2c06-9dd64ebf12cd
# ╠═52b234b0-1389-11eb-3d14-717752260ec5
# ╠═5abe2470-1389-11eb-1a3e-43f6103379d2
# ╠═5aba7e24-0424-45f2-9716-3b32a71fc610
# ╠═213d2b7b-b742-4274-9bb0-e029aec6f892
# ╠═91b52c67-de20-4bfc-9da6-3e04ed73b990
# ╠═0d07f057-9012-42ad-bf52-31a0f14614df
# ╠═2cb0db8c-11e6-49e3-baf7-9fdda352a26a
# ╠═adbd72cb-4dcf-490b-bbcb-1d681358c455
# ╠═9ae4a166-f3ee-4850-89ff-4c0a41bae48c
# ╠═4a832be5-26c0-46fa-bcbb-3afb4eb0cf2e
# ╟─92d8992a-f55e-4576-a1aa-c8b0de9e5806
# ╠═54e5e9c4-5ad2-45f9-a3ed-0aaced61663d
# ╠═8794b613-ddab-4fa9-82e6-eec4192705dd
# ╠═02d632a3-3ab8-4b02-ba94-709880df6313
# ╟─aba2f7f6-8690-4ede-8282-ad925e4aae8d
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
