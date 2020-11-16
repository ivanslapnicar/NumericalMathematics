### A Pluto.jl notebook ###
# v0.12.10

using Markdown
using InteractiveUtils

# ╔═╡ aa7480bd-e15a-4f07-bcea-f3f2c9f3f5fb
using LinearAlgebra

# ╔═╡ 6075a4cb-5931-49d1-987c-ffc0f40ebb12
using DelimitedFiles

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
* history, context - cookies, storing data (about You), [200+ parameters](http://backlinko.com/google-ranking-factors)

## PageRank

* Graph Theory and Linearna Algebra
* [C. Moler, Google PageRank](https://www.mathworks.com/moler/exm/chapters/pagerank.pdf)

Some programs:

* https://github.com/purzelrakete/Pagerank.jl

* https://gist.github.com/domluna/2b9358ccc89fee7d5e26

We try example from  Moler's paper.
"""

# ╔═╡ d761f999-c402-49c4-bb4f-b52c23475db1
begin
	i = vec([ 2 6 3 4 4 5 6 1 1])
	j = vec([ 1 1 2 2 3 3 3 4 6])
end

# ╔═╡ daf1a60f-9827-418b-adbe-06d2046b1c96
begin
	using SparseArrays
	G=sparse(i,j,1.0)
end

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
## Idea

Let us start a random walk from the vector $x_0=\begin{bmatrix} 1/n \\ 1/n \\ \vdots \\ 1/n \end{bmatrix}$ (or some other!).

The subsequent vectors are

$$\begin{aligned}
x_1&=A\cdot x_0 \cr
x_2&=A\cdot x_1 \cr
x_3&=A\cdot x_2\cr
& \ \vdots
\end{aligned}$$

Map $A(x)=Ax$ is not a contraction in the sence of the Banach Fixed Point Theorem since $\|A\|_1=1$, but is can be shown that it has a fixed point. Also, id $x\geq 0$ (componentwise), then $\|Ax\|_1=\|x\|_1$.

When the iterations _stabilize_ :

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

# ╔═╡ 47394960-0f02-11eb-1ddf-cb6b81f096b4
W=readdlm("../files/web-Stanford.txt",Int,comments=true)

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

# ╔═╡ 2188e21d-f3bb-41f3-8a9c-c425ba4f0887


# ╔═╡ Cell order:
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
# ╠═aa7480bd-e15a-4f07-bcea-f3f2c9f3f5fb
# ╠═fece3020-0f09-11eb-0f69-237286bd58af
# ╠═739c238c-03db-4ee6-9fb7-f8e5b93282f8
# ╠═6c6a8ce2-5483-45ed-b5c8-61e924b3eb1c
# ╠═cb04da5e-0f08-11eb-21b9-8fdaea539145
# ╠═870c8ccd-7e2c-489c-957f-fc34651bb65f
# ╠═738462b0-62a9-4aed-8ca7-687fb51d52e2
# ╠═2e73e3d1-cc59-4977-a89e-bb9d1c2eb89f
# ╠═023e6d22-6505-4211-8fec-55ae732405bc
# ╠═fae1bfcb-ef52-4cd8-a066-cf138c8697f8
# ╟─5a02f8ad-3f97-4201-b903-9ed789721f81
# ╠═6075a4cb-5931-49d1-987c-ffc0f40ebb12
# ╠═47394960-0f02-11eb-1ddf-cb6b81f096b4
# ╠═fb14f5c1-2b83-4534-95a3-4647d8e50738
# ╠═573a625a-ad1c-4133-bae3-342a7501b492
# ╠═83fdc63f-4aac-45c0-a226-87a4830f697e
# ╠═b473be54-b1a2-4f34-80d8-742386c9535b
# ╠═d94eacf3-0008-4e8c-a5cd-a92fa1fd76d4
# ╠═12eb4b3f-6fe1-4bba-9556-4378eab6e191
# ╠═a445d5e2-adef-4d13-b0b3-0af37f7039d6
# ╠═2188e21d-f3bb-41f3-8a9c-c425ba4f0887
