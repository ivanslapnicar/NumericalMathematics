### A Pluto.jl notebook ###
# v0.19.20

using Markdown
using InteractiveUtils

# ╔═╡ bc0768e0-ca87-46da-9260-7b2caabcdf7a
using PlutoUI, SparseArrays, LinearAlgebra, DelimitedFiles

# ╔═╡ e1fea405-5099-4c2c-8fd4-a02469ead9f6
TableOfContents(title="📚 Table of Contents", aside=true)

# ╔═╡ c1dc1294-a3a7-431f-8ecf-ea3c8448c44a
md"""
# PageRank
"""

# ╔═╡ 85992dca-d40f-461a-8984-ac57ddfff970
md"""
## Age of Search

google (and others)


* [55 billion pages](http://www.worldwidewebsize.com/), [4 billion searches daily](http://www.internetlivestats.com/google-search-statistics/)
* __PageRank__
* History, context - cookies, storing data (about You), [200+ parameters](http://backlinko.com/google-ranking-factors)

The age of search was followed by the __age of recommendation__, and, since 2022, we have the __age of chat__ (?). 

# Transition matrix (and graph)

* Graph Theory and Linear Algebra
* [C. Moler, Google PageRank](https://www.mathworks.com/moler/exm/chapters/pagerank.pdf)

Some programs:

* https://github.com/purzelrakete/Pagerank.jl

* https://gist.github.com/domluna/2b9358ccc89fee7d5e26

We try the example from  Moler's paper.
"""

# ╔═╡ d761f999-c402-49c4-bb4f-b52c23475db1
begin
	i = vec([ 2 6 3 4 4 5 6 1 1])
	j = vec([ 1 1 2 2 3 3 3 4 6])
end

# ╔═╡ daf1a60f-9827-418b-adbe-06d2046b1c96
G=sparse(i,j,1.0)

# ╔═╡ f874ba0c-9d9e-437b-9709-723ce34756c3
Matrix(G)

# ╔═╡ f9abf219-8880-4d45-a556-ba20ba114a56
begin
	G₁=similar(G)
	c=sum(G,dims=1)
	n=size(G,1)
	for j=1:n
	    if c[j]>0
	        G₁[:,j]=G[:,j]/c[j]
	    end
	end
end

# ╔═╡ 726e69a4-1606-466f-a927-2effba6bcaf2
Matrix(G₁)

# ╔═╡ 7c96e1d5-c076-4638-8c95-fe5bc8ab8936
md"""
* probability to follow some link is $p$
* probability to visit random page is $1-p$
* google uses $p=0.85$ ?
"""

# ╔═╡ ae67506f-65b2-4b67-af2e-3fc07f5044e9
begin
	p=0.85
	δ=(1-p)/n
end

# ╔═╡ ef68750f-be24-4d53-9d60-9715d4064ecf
z = ((1-p)*(c.!=0) + (c.==0))/n

# ╔═╡ 420af1cd-becc-4005-8534-96a6aaddbde5
A=p*G₁+ones(n)*z

# ╔═╡ 9d3b7528-4035-4304-bee9-a9407bced36f
sum(A,dims=1)

# ╔═╡ 90edb7b5-c882-41a6-a48b-ba15373f2283
md"""
# Random walk

Let us start a random walk from the vector $x_0=\begin{bmatrix} 1/n \\ 1/n \\ \vdots \\ 1/n \end{bmatrix}$ (or some other!).

The subsequent vectors are

$$\begin{aligned}
x_1&=A\cdot x_0 \cr
x_2&=A\cdot x_1 \cr
x_3&=A\cdot x_2\cr
& \ \vdots
\end{aligned}$$

Map $A(x)=Ax$ is not a contraction in the sence of the Banach Fixed Point Theorem since $\|A\|_1=1$, but is can be shown that it has a fixed point. Also, id $x\geq 0$ (componentwise), then $\|Ax\|_1=\|x\|_1$.

When the iterations __stabilize__:

$$
A\cdot x\approx x,$$

then $x[i]$ is the __rank__ of the page $i$.
"""

# ╔═╡ fece3020-0f09-11eb-0f69-237286bd58af
function PageRank(G₁::SparseMatrixCSC{Float64,Int64},steps::Int)
	G=copy(G₁)
	p=0.85
	c=sum(G,dims=1)/p
	n=size(G,1)
	for i=1:n
	    G.nzval[G.colptr[i]:G.colptr[i+1]-1]./=c[i]
	end
	e=ones(n)
	x=e/n
	z = vec(((1-p)*(c.!=0) + (c.==0))/n)
	for j=1:steps
	    x=G*x.+(z⋅x)
	end
	return x
end

# ╔═╡ a62cea43-948e-43d7-9df2-82ff096d04ba
md"
__We need to understand and use the CSC format.__
"

# ╔═╡ 739c238c-03db-4ee6-9fb7-f8e5b93282f8
fieldnames(typeof(G))

# ╔═╡ 6c6a8ce2-5483-45ed-b5c8-61e924b3eb1c
G

# ╔═╡ cb04da5e-0f08-11eb-21b9-8fdaea539145
Matrix(G)

# ╔═╡ 870c8ccd-7e2c-489c-957f-fc34651bb65f
G.colptr

# ╔═╡ 738462b0-62a9-4aed-8ca7-687fb51d52e2
G.nzval

# ╔═╡ 2e73e3d1-cc59-4977-a89e-bb9d1c2eb89f
G.rowval

# ╔═╡ 023e6d22-6505-4211-8fec-55ae732405bc
# Starting vector
x=ones(n)/n

# ╔═╡ fae1bfcb-ef52-4cd8-a066-cf138c8697f8
PageRank(G,15)

# ╔═╡ 5a02f8ad-3f97-4201-b903-9ed789721f81
md"""
## [Stanford web graph](http://snap.stanford.edu/data/web-Stanford.html)

A bigger test problem.
"""

# ╔═╡ eb9f4ac4-53fb-4a92-b33e-1d0074d83edb
# Create directory
if !isdir("files")
	mkdir("files")
end

# ╔═╡ b46547c4-6a2d-45b5-bad3-75ebf8b4451f
# Download the test file
download("https://ivanslapnicar.github.io/NumericalMathematics/files/web-Stanford.txt","./files/web-Stanford.txt")

# ╔═╡ 47394960-0f02-11eb-1ddf-cb6b81f096b4
W=readdlm("./files/web-Stanford.txt",Int,comments=true)

# ╔═╡ fb14f5c1-2b83-4534-95a3-4647d8e50738
#?sparse

# ╔═╡ 573a625a-ad1c-4133-bae3-342a7501b492
S=sparse(W[:,2],W[:,1],1.0)

# ╔═╡ 83fdc63f-4aac-45c0-a226-87a4830f697e
@time x100=PageRank(S,100);

# ╔═╡ b473be54-b1a2-4f34-80d8-742386c9535b
x101=PageRank(S,101);

# ╔═╡ d94eacf3-0008-4e8c-a5cd-a92fa1fd76d4
maximum(abs,(x101-x100)./x101)

# ╔═╡ 12eb4b3f-6fe1-4bba-9556-4378eab6e191
# Ranks
sort(x100,rev=true)

# ╔═╡ a445d5e2-adef-4d13-b0b3-0af37f7039d6
# Pages
sortperm(x100,rev=true)

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
DelimitedFiles = "8bb1440f-4735-579b-a4ab-409b98df4dab"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

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

[[DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

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

[[SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

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
# ╠═bc0768e0-ca87-46da-9260-7b2caabcdf7a
# ╠═e1fea405-5099-4c2c-8fd4-a02469ead9f6
# ╟─c1dc1294-a3a7-431f-8ecf-ea3c8448c44a
# ╟─85992dca-d40f-461a-8984-ac57ddfff970
# ╠═d761f999-c402-49c4-bb4f-b52c23475db1
# ╠═daf1a60f-9827-418b-adbe-06d2046b1c96
# ╠═f874ba0c-9d9e-437b-9709-723ce34756c3
# ╠═f9abf219-8880-4d45-a556-ba20ba114a56
# ╠═726e69a4-1606-466f-a927-2effba6bcaf2
# ╟─7c96e1d5-c076-4638-8c95-fe5bc8ab8936
# ╠═ae67506f-65b2-4b67-af2e-3fc07f5044e9
# ╠═ef68750f-be24-4d53-9d60-9715d4064ecf
# ╠═420af1cd-becc-4005-8534-96a6aaddbde5
# ╠═9d3b7528-4035-4304-bee9-a9407bced36f
# ╟─90edb7b5-c882-41a6-a48b-ba15373f2283
# ╠═fece3020-0f09-11eb-0f69-237286bd58af
# ╟─a62cea43-948e-43d7-9df2-82ff096d04ba
# ╠═739c238c-03db-4ee6-9fb7-f8e5b93282f8
# ╠═6c6a8ce2-5483-45ed-b5c8-61e924b3eb1c
# ╠═cb04da5e-0f08-11eb-21b9-8fdaea539145
# ╠═870c8ccd-7e2c-489c-957f-fc34651bb65f
# ╠═738462b0-62a9-4aed-8ca7-687fb51d52e2
# ╠═2e73e3d1-cc59-4977-a89e-bb9d1c2eb89f
# ╠═023e6d22-6505-4211-8fec-55ae732405bc
# ╠═fae1bfcb-ef52-4cd8-a066-cf138c8697f8
# ╟─5a02f8ad-3f97-4201-b903-9ed789721f81
# ╠═eb9f4ac4-53fb-4a92-b33e-1d0074d83edb
# ╠═b46547c4-6a2d-45b5-bad3-75ebf8b4451f
# ╠═47394960-0f02-11eb-1ddf-cb6b81f096b4
# ╠═fb14f5c1-2b83-4534-95a3-4647d8e50738
# ╠═573a625a-ad1c-4133-bae3-342a7501b492
# ╠═83fdc63f-4aac-45c0-a226-87a4830f697e
# ╠═b473be54-b1a2-4f34-80d8-742386c9535b
# ╠═d94eacf3-0008-4e8c-a5cd-a92fa1fd76d4
# ╠═12eb4b3f-6fe1-4bba-9556-4378eab6e191
# ╠═a445d5e2-adef-4d13-b0b3-0af37f7039d6
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
