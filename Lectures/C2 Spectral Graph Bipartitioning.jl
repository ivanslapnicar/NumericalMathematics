### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# ╔═╡ 77c3ba33-a12d-4e48-bbeb-727166771e63
begin
	import Pkg
    Pkg.activate(mktempdir())
    Pkg.add([	
		Pkg.PackageSpec(name="PlutoUI"),
		Pkg.PackageSpec(name="LightGraphs"),
		Pkg.PackageSpec(name="GraphPlot"),
		Pkg.PackageSpec(name="SparseArrays"),
		Pkg.PackageSpec(name="Plots"),
		Pkg.PackageSpec(name="Distances"),
		Pkg.PackageSpec(name="Arpack"),
		Pkg.PackageSpec(name="Images")
    ])
end

# ╔═╡ 0cd92d78-1da6-435b-ad8e-60af29502e1d
begin
	# Necessary packages
	using PlutoUI, Random, LinearAlgebra, LightGraphs, GraphPlot
	using SparseArrays, Plots, Distances, Arpack, Images
end

# ╔═╡ 7ef4b0b3-9d0a-4356-bd16-14487f14e0ab
TableOfContents(title="📚 Table of Contents", aside=true)

# ╔═╡ cb61f761-b2ee-4818-9149-f85d95b76a1e
md"""
# Spectral Graph Bipartitioning

Many data clustering problems can be interpreted as clustering of vertices of graphs. __Graph bipartitioning problem__ is to partition vertices into subsets such that the connections within subsets are stronger than the connections between different subsets.

Partition of the vertices into two subsetts is done according to signs of the eigenvectors of the second smallest eigenvalue of the Laplacian matrix. 

__Prerequisites__

The reader should be familiar with the basic graph theory, linear algebra and, in particular,  eigenvalues and eigenvectors.
 
__Competences__

The reader should be able to apply graph spectral bipartitioning and recursive bipartitioning to data clustering problems.

__Credits.__ The notebook was initially derived from M.Sc. Thesis of Ivančica Mirošević.
"""

# ╔═╡ 5926c79a-6038-4dbc-b49a-5b05d02250a0
md"""
# Graphs

For more details, see [W. H. Haemers, Matrices and Graphs, in L. Hogben, Ed., 'Handbook of Linear Algebra', pp. 39.1-39.14, CRC Press, Boca Raton, 2014.](https://www.routledge.com/Handbook-of-Linear-Algebra/Hogben/p/book/9781138199897) and [S. Butler and F. Chung, Spectral Graph Theory, ibid., pp. 47.1-47.6](https://www.routledge.com/Handbook-of-Linear-Algebra/Hogben/p/book/9781138199897) and the references therein.

## Definitions

A __weighted graph__ is an ordered triplet $G=(V,E,\omega)$, where $V=\{1,2,3,...,n\}$ is the set of __vertices__ , $E=\{(i,j)\}$ is a set of __edges__ connecting vertices, and $\omega$ is a set of __weights__ of edges. We assume $G$ is undirected.
"""

# ╔═╡ a66a6e80-1f80-11eb-15f8-71ffaa5050e7
md"
__Adjacency matrix__ of graph $G$ is the matrix $A$ defined as

$$A_{ij}=\begin{cases} 1 \quad \textrm{if}\ (i,j)\in E, \\
0\quad  \textrm{otherwise} \end{cases}.$$

__Weight matrix__ of graph $G$ is the matrix $W$ defined as

$W_{ij}=\begin{cases} \omega(e) \quad \textrm{if}\ e=(i,j)\in E, \\
0\quad  \textrm{otherwise} \end{cases}.$

__Laplacian matrix__ of graph $G$ is the matrix 

$L=D-W,$ 

where 
$D=\mathop{\mathrm{diag}}(d_1,d_2,\ldots,d_n)$ with $d_i=\sum_{k=1}^n W_{ik}$ for $i=1,\ldots,n$.

__Normalized Laplacian matrix__ is the matrix

$L_n=D^{-1/2} L D^{-1/2}\equiv D^{-1/2} (D-W) D^{-1/2}$ 

(__diagonally scaled $L$__).

__Incidence matrix__ of graph $G$ is the $|V|\times |E|$ matrix $I_G$. Each row of $I_G$ corresponds to a vertex of $G$ and each column corresponds to an edge of $G$.
In the column corresponding to en edge $e=(i,j)$, all elements are zero except the ones in the $i$-th and $j$-th row, which are equal to $\sqrt{\omega(e)}$ and $-\sqrt{\omega(e)}$, respectively.
"

# ╔═╡ d4cf12ce-12fb-4562-bbea-e9b85cba94c2
md"""
### Examples

Graph types and algorithms are implemented in the package [LightGraphs.jl](https://github.com/JuliaGraphs/LightGraphs.jl). 

Plotting graphs is done by the packages
[GraphPlot.jl](https://github.com/JuliaGraphs/GraphPlot.jl).

As a small inconvenience, we can only plot unweighted graphs and plot weights as edge labels.
"""

# ╔═╡ 87cca1d3-b944-4d74-bf4c-1562f34cf91a
begin
	# Sources, targets and weights
	n=7
	sn=[1,1,1,2,2,3,3,3,5,5,6]
	tn=[2,3,4,4,5,4,6,7,6,7,7]
	wn=[2,3,4,7,1,3,2,1,7,3,5]
	[sn tn wn]
end

# ╔═╡ 3d7b4a01-9a7d-4111-a45d-1d31f1ef6278
begin
	# Create the graph
	G=Graph(n)
	for i=1:length(sn)
	    add_edge!(G,sn[i],tn[i])
	end
	G
end

# ╔═╡ 3b5ed1fa-75b0-43fb-a918-4c2e8f52d039
# What is the optimal bipartition?
gplot(G, nodelabel=1:n, edgelabel=wn)

# ╔═╡ 8f7d1f6c-ee4c-407f-80d9-d2abc2af949e
begin
	# We define some functions
	function WeightMatrix(src::Array,dst::Array,weights::Array)
	    n=nv(G)
	    sparse([src;dst],[dst;src],[weights;weights],n,n)
	end
	
	Laplacian(W::AbstractMatrix)=spdiagm(0=>vec(sum(W,dims=2)))-W
	
	function NormalizedLaplacian(L::AbstractMatrix)
	    D=1.0./sqrt.(diag(L))
		Diagonal(D)*L*Diagonal(D)
	end
end

# ╔═╡ b7f127d8-e8cd-4771-8091-c88f9ffe6b08
W=WeightMatrix(sn,tn,wn)

# ╔═╡ 3f5447bb-95cd-4520-b0b9-e85f33b323bd
Matrix(W)

# ╔═╡ ecb2a46e-f9dc-4e57-812f-0a2788afb202
begin
	L=Laplacian(W)
	Matrix(L)
end

# ╔═╡ 86e5f00f-13ee-4867-b0cd-187f2c782a46
Lₙ=NormalizedLaplacian(L)

# ╔═╡ cf26dbe9-3c4f-4c95-88ef-f197f360d759
Matrix(Lₙ)

# ╔═╡ 8fcd2440-4013-4b37-a747-2c0c7f6d4aeb
issymmetric(Lₙ)

# ╔═╡ 0c8d3572-3436-446f-952a-b19d0370ab38
# Let us compute the incidence matrix
function IncidenceMatrix(G::Graph, weights::Array)
    A=zeros(nv(G),ne(G))
    k=1
    for a in edges(G)
        A[a.dst,k]=sqrt.(weights[k])
        A[a.src,k]=-sqrt(weights[k])
        k+=1
    end
    A
end

# ╔═╡ 1437aadc-1c26-4eff-bb22-9cf9c15e1f5f
Iᵧ=IncidenceMatrix(G,wn)

# ╔═╡ 6529324b-49f7-4f01-b5cd-33668979b56b
md"""
## Facts

1.  $L=I_{G}I_{G}^{T}$.

2.  $L$ is symmetric PSD matrix.

3.  $L\mathbf{1}=0$ for $\mathbf{1}=[1,...,1]^{T}$, thus $0$ is an eigenvalue of $L$  and $\mathbf{1}$ is the corresponding eigenvector.

4. If $G$ has $c$ connected components, then $L$ has $c$ eigenvalues equal to $0$.

5. For every $x\in \mathbb{R}^{n}$, it holds
$x^{T}L x=\sum\limits_{i<j}W_{ij}(x_{i}-x_{j})^{2}$.

6. For every $x\in\mathbb{R}^{n}$ and $\alpha,\beta\in\mathbb{R}$, it holds
$(\alpha x+\beta \mathbf{1})^{T} L (\alpha x+\beta \mathbf{1}) 
=\alpha^{2} x^{T}L x$.

7. Assume that the eigenvalues of $L$ are increasingly ordered. Then,

$$
0=\lambda_1(L)\leq \lambda_2(L)\leq \cdots \leq\lambda_{n}(L)\leq 
2\max\limits_{i=1,\cdots ,n}d_{i}.$$

8.  $\sigma(L_n) \subseteq [0,2]$.
"""

# ╔═╡ a9ba64a6-bd17-4b86-acb4-e7c7d8f87dc1
md"""
### Examples
"""

# ╔═╡ 6a6f5db5-7dd4-4a26-ac1a-ef3ca3a5d6be
# Fact 1
norm(L-Iᵧ*Iᵧ')

# ╔═╡ e11bc23d-03e2-4780-9488-93f0adfc2fcb
# Facts 2 and 7
issymmetric(L), eigs(L)[1], 2*maximum(diag(L))

# ╔═╡ 7192f7d0-28b9-11eb-28dd-8db0019504ca
eigen(Matrix(L))

# ╔═╡ e5890dca-368f-437e-bc45-6a395714dae2
# Fact 3
L*ones(n)

# ╔═╡ 61a978a7-b2d2-43e7-b3f7-95bf4c4f54b9
begin
	# Fact 5
	x=rand(n)
	x'*L*x, sum([W[i,j]*(x[i]-x[j])^2 for i=1:n, j=1:n])/2
end

# ╔═╡ da7f85ed-73b0-44b9-a97f-4a8ed5e9a1d1
begin
	# Fact 6
	α,β=rand(),rand()
	(α*x+β*ones(n))'*L*(α*x+β*ones(n)), α^2*x'*L*x
end

# ╔═╡ fed7f592-be03-4baf-8ace-ff9ee78fee71
# Fact 8
eigvals(Matrix(Lₙ))

# ╔═╡ 9f31928e-c933-43bb-9a5c-bf55f29737cd
md"""
# Bipartitioning

## Definitions

Let $\pi=\{V_{1},V_{2}\}$ be a partition of $V$ with $V_1,V_2\neq \emptyset$.

__Cut__ of partition $\pi$ is the sum of weights of all 
edges between $V_1$ and $V_2$, 

$$\mathop{\mathrm{cut}}(\pi)\equiv \mathop{\mathrm{cut}}(V_1,V_2)=\sum\limits_{{\displaystyle i\in V_{1} \atop \displaystyle j\in V_{2}}}W_{ij}.$$

__Weight__ of vertex $i\in V$ is the sum of the weights of all egdges emanating from $i$,
$\omega(i)=\sum\limits_{j=1}^{n}W_{ij}$.

__Weight__ of a subset $\bar V\subset V$ is the sum of the weights of all vertices in $\bar V$, 

$\omega(\bar V)=\sum\limits_{\displaystyle i\in\bar V} \omega(i)$.

__Proportional cut__ of partition $\pi$ is

$$
\mathop{\mathrm{pcut}}(\pi)=\displaystyle\frac{\mathop{\mathrm{cut}}(\pi)}{|V_{1}|}+\frac{\mathop{\mathrm{cut}}(\pi)}{|V_{2}|}.$$

__Normalized cut__ of partition $\pi$ is

$$
\mathop{\mathrm{ncut}}(\pi)=\displaystyle\frac{\mathop{\mathrm{cut}}(\pi)}{\omega(V_{1})}+\frac{\mathop{\mathrm{cut}}(\pi)}{\omega(V_{2})}.$$
"""

# ╔═╡ 1e6d4c72-f25d-4bf5-99d6-ae0446be25dd
md"""
### Example

Consider the following partitions (all edges have unit weights):

 $(load(\"./files/cut2.png\"))

Left partition is $\pi$, right partition is $\pi'$.

|     Cut \ Partition  | $\pi$            |  $\pi'$     |
| ------- | ---------------- | ------------|
| $\mathop{\mathrm{cut}}$  |  $2$        |     $3$        |
| $\mathop{\mathrm{pcut}}$ | $\frac{2}{1}+\frac{2}{11}=2.18$|    $\frac{3}{6}+\frac{3}{6}=1$      |
| $\mathop{\mathrm{ncut}}$ | $\frac{2}{2}+\frac{2}{50}=1.04$ |  $\frac{3}{27}+\frac{3}{25}=0.23$|
"""

# ╔═╡ 02035eed-6d49-4253-a755-68a999e7e90e
md"""
## Facts

1. The informal description of the bipartitioning problem can be formulated as two problems, 

$$\mathop{\textrm{arg min}}\limits_{\pi} \mathop{\mathrm{pcut}}(\pi) \quad \textrm{or} \quad \mathop{\textrm{arg min}}\limits_{\pi} \mathop{\mathrm{ncut}}(\pi).$$ 

The first problem favors partitions into subsets with similar numbers of vertices, while the second problem favors partitions into subsets with similar weights.

2. Both problems are NP-hard.

3. __Approximate solutions can be computed by suitable relaxations in $O(n^2)$ operations.__

4. The partition $\pi$ is defined by the vector $y$ such that

$$
y_{i}=
\begin{cases}
\frac{1}{2} & \text{for } i\in V_1 \\
-\frac{1}{2} & \text{for } i\in V_2
\end{cases}$$

The proportional cut problem can be formulated as the  __discrete proportional cut__ problem

$$
\underset{\displaystyle \big|\mathbf{y}^{T}\mathbf{1} \big|\leq \beta}
{\min\limits_{\displaystyle y_{i}\in \{-\frac{1}{2},\frac{1}{2}\}}}
\frac{1}{2}\sum_{i,j}(y_{i}-y_{j})^{2}W_{ij}.$$

Parameter $\beta$ controls the number of vertices in each subset.

5. The normalized cut problem can be formulated as the __discrete normalized cut__ problem

$$
\underset{\displaystyle \big|y^{T}D\mathbf{1} \big|\leq \beta}
{\min\limits_{\displaystyle y_{i}\in \{-\frac{1}{2},\frac{1}{2}\}}}
\frac{1}{2}\sum_{i,j}(y_{i}-y_{j})^{2}W_{ij}.$$

Parameter $\beta$ controls the weights of each subset.

6. Using the Fact 5 above, the discrete proportional cut problem can be formulated as the __relaxed proportional cut__ problem

$$
\underset{\displaystyle y^{T}y=1}{\underset{\displaystyle \big| y^{T}\mathbf{1} \big|
\leq 2\frac{\beta}{\sqrt{n}}}
{\min\limits_{\displaystyle y\in \mathbb{R}^{n}}}} y^{T}L y.$$

Similarly, the discrete normalized cut problem can be formulated as the __relaxed normalized cut__ problem

$$
\underset{\displaystyle y^{T}Dy=1}{\underset{\displaystyle \big| y^{T}D\mathbf{1}\big|
\leq \displaystyle \frac{\beta}{\sqrt{\theta n}}}{\min\limits_{\displaystyle y\in
\mathbb{R}^{n}}}}y^{T}L_n y.$$

7. __The Main Theorem.__ Let $A\in \mathbb{R}^{n\times n}$ be a symmetric matrix with eigenvalues $\lambda _{1}<\lambda _{2}<\lambda_{3}\leq \cdots \leq \lambda _{n}$ and let $v^{[1]},v^{[2]},\ldots,v^{[n]}$ be the corresponding eigenvectors. For the fixed $0\leq \alpha <1$, the solution of the problem

$$
\underset{\displaystyle y^{T}y=1}{\underset{\displaystyle \left|y^{T}v^{[1]}\right|\leq \alpha}
{\min\limits_{\displaystyle y\in \mathbb{R}^{n}}}} y^{T}Ay$$

is $y=\pm \alpha v^{[1]}\pm \sqrt{1-\alpha^{2}}v^{[2]}$. 

For the proof see [D. J. Higham and M. Kibble, A Unified View of Spectral Clustering, Theorem 3.1, p. 7](http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.115.1591&rep=rep1&type=pdf).

8. For $0\leq \beta <\frac{n}{2}$, the solution of the relaxed proportional cut problem is

$$
y=\pm \frac{2\beta}{n\sqrt{n}}\mathbf{1}\pm
\sqrt{1-4\frac{\beta ^{2}}{n^{2}}}v^{[2]},$$

where $v^{[2]}$ is an eigenvector corresponding to $\lambda_2(L)$. $v^{[2]}$ the __Fiedler vector__. Since the first summand carries no information, $V$ is partitioned according to the signs of the components of $v^{[2]}$:

$$
V_{1}=\{i:v^{[2]}_i <0\}, \quad V_{2}=\{i:v^{[2]}_i \geq 0\}.$$

_Notice that the value of $\beta$ is irrelevant for the solution._

9. For $0\leq \beta <\sqrt{\theta n}\left\Vert D^{\frac{1}{2}}\mathbf{1} \right\Vert _{2},$ the solution of the relaxed normalized cut problem is

$$
y=\pm \frac{\beta }{\sqrt{\theta n}\left\Vert
D^{\frac{1}{2}} \mathbf{1}\right\Vert _{2}^{2}}\mathbf{1}\pm
\sqrt{1-\frac{\beta ^{2}}{ \theta n\left\Vert
D^{\frac{1}{2}}\mathbf{1}\right\Vert _{2}^{2}}}D^{-\frac{1
}{2}} v_n^{[2]},$$

where $v_n^{[2]}$ is an eigenvector corresponding to $\lambda_2(L_n)$. $V_n$ is partitioned according to the signs of the components of $v_n^{[2]}$, as above.

10. Neither of the relaxed algorithms is guaranteed to solve exactly the true proportional / normalized cut problem. However, the computed solutions are in the right direction. Whether to use proportional or normalized cut formulation, depends upon the specific problem.  
"""

# ╔═╡ c188b7cf-03e7-40b3-bf90-773978c7b410
# Voila!
eigs(L,nev=2,which=:SM, v0=ones(n))

# ╔═╡ f60b51b6-79c3-4a3f-8bfb-8f198cdb92a8
# For the normalized cut
eigs(Lₙ,nev=2,which=:SM, v0=ones(n))

# ╔═╡ bfac7a55-0d82-4549-b999-7db3272eaa20
md"""
### Concentric rings

A __complete graph__ has edges connecting each pair of vertices.

To a set of points $X=\{x_{1},x_{2},\cdots ,x_{m}\}$ , where $x_{i}\in\mathbb{R}^{n}$, we assign a weighted complete graph $G=(V,E)$ with $m$ vertices, where the vertex $j\in V$ corresponds to the point $x_j\in X$.

The main idea is to assign weight of an edge $e=(i,j)$ which reflects the distance between $x_i$ and $x_j$, something like $\omega(e)=\displaystyle\frac{1}{\mathop{\mathrm{dist}}(x_i,x_j)}$.

However, this has to be implemented with care. For example, using simple Euclidean distance yield the same results as the function `kmeans()`. In this example we use Gaussian kernel, that is

$$
\omega(e)=e^{\displaystyle -\|x_i-x_j\|_2^2/\sigma^2},$$

where the choice of $\sigma$ is based on experience.

The computation of various distances is implemented in the package [Distances.jl](https://github.com/JuliaStats/Distances.jl).

We will construct the Laplace matrix directly.
"""

# ╔═╡ 80ed91ad-6f9a-4f12-a1ee-bc8dc6b58168
begin
	# Two concentric circles
	k=2
	# Center
	Random.seed!(541)
	# center=[rand(-5:5),rand(-5:5)]
	center=[0,0]
	# Radii
	radii=randperm(10)[1:k]
	# Number of points in circles
	sizes=rand(1000:2000,k)
	center,radii,sizes
end

# ╔═╡ 2d8689c9-1b9e-4360-865c-23c2c8c5dcb9
begin
	# Generate points
	X=Array{Float64}(undef,2,sum(sizes))
	csizes=cumsum(sizes)
	# Random angles
	ϕ=2*π*rand(sum(sizes))
	for i=1:csizes[1]
		X[:,i]=center+radii[1]*[cos(ϕ[i]);sin(ϕ[i])] + (rand(2).-0.5)/50
	end
	for j=2:k
		for i=csizes[j-1]+1:csizes[j]
			X[:,i]=center+radii[j]*[cos(ϕ[i]);sin(ϕ[i])] + (rand(2).-0.5)/50
		end
	end
	scatter(X[1,:],X[2,:],title="Concentric rings", aspect_ratio=1,label="Points")
end

# ╔═╡ fa4ae713-a0b2-46ce-b5f7-558b26d82019
# Weight matrix
W₁=1 ./pairwise(SqEuclidean(),X)

# ╔═╡ 760e8ca2-596f-435f-befb-140b3d204f7f
begin
	# Laplacian matrix
	m=csizes[end]
	for i=1:m
	    W₁[i,i]=0
	end
	L₁=Diagonal(vec(sum(W₁,dims=2)))-W₁
	# Check Fact 3
	norm(L₁*ones(m))
end

# ╔═╡ 75956488-76ab-4262-a32f-0c53e1b3dc17
# Notice λ₁=0
E=eigs(L₁,nev=2,which=:SM, v0=ones(m))

# ╔═╡ a55260dd-f05d-4a65-9fdf-cb3fd4638857
begin
	# Define clusters
	C=ones(Int64,m)
	C[findall(E[2][:,2].>0)].=2
	C
end

# ╔═╡ ba164ff3-5ef4-4480-8af2-9dd9a99b1e8d
# Yet another plotting function
function plotKpartresult(C::Vector,X::Array)
	scatter(aspect_ratio=1)
    k=maximum(C)
    for j=1:k
        scatter!(X[1,findall(C.==j)],X[2,findall(C.==j)],label="Cluster $j")
    end
	scatter!(aspect_ratio=1)
end

# ╔═╡ 46de28bc-73ed-495c-a2e0-f565c1ce5651
plotKpartresult(C,X)

# ╔═╡ e004e4bf-ca04-4b41-947b-9f1ecd058abb
md"""
This is the same partitioning as obtained by `kmeans()`. Let us try Gaussian kernel. A rule of thumb is: if rings are close, use $\sigma<1$, if rings are apart, use $\sigma>1$.
"""

# ╔═╡ bb9d4d98-3aff-439e-b8a3-8347cf345839
begin
	σ=0.7 # 0.1
	W₂=exp.(-pairwise(SqEuclidean(),X)/σ^2)-I
	L₂=Diagonal(vec(sum(W₂,dims=2)))-W₂
	E₂=eigs(L₂,nev=2,which=:SM, v0=ones(m))
	C₂=ones(Int64,m)
	C₂[findall(E₂[2][:,2].>0)].=2
	plotKpartresult(C₂,X)
end

# ╔═╡ 012fa54a-5171-4980-a0ba-474a22c25981
md"""
# Recursive bipartitioning

## Definitions

Let $G=(V,E)$ be a weighted graph with weights $\omega$.

Let $\pi_k =\{V_{1},V_{2},...,V_{k}\}$ be a $k$-partition of $V$, with $V_i\neq \emptyset$ for $i=1,\ldots,k$.

The previous definition of $cut(\pi)\equiv cut(\pi_2)$ extends naturally to $k$-partition.

A __cut__ of a partition $\pi_k$ is 

$$
\mathop{\mathrm{cut}}(\pi_k)=\sum\limits_{\displaystyle i<j} \mathop{\mathrm{cut}}(V_{i},V_{j}),$$

where $\mathop{\mathrm{cut}}(V_{i},V_{j})$ is interpreted as a cut of the bipartition of the subgraph of $G$ with vertices $V_1\cup V_2$.

__Proportional cut__ of a partition $\pi_k$ is

$$
\mathop{\mathrm{pcut}}(\pi_k)=\underset{i<j}{\sum\limits_{i,j=1}^{k}} \left(
\frac{\mathop{\mathrm{cut}}(V_{i},V_{j})}{|V_{i}|}+\frac{\mathop{\mathrm{cut}}(V_{i},V_{j})}{|V_{j}|}\right) =
\sum_{i=1}^{k}\frac{\mathop{\mathrm{cut}}(V_{i},V\backslash V_{i})}{|V_{i}|}.$$

__Normalized cut__ of a partition $\pi_k$ is

$$
\mathop{\mathrm{ncut}}(\pi_k)=\underset{i<j}{\sum\limits_{i,j=1}^{k}} \left(
\frac{\mathop{\mathrm{cut}}(V_{i},V_{j})}{\omega(V_{i})}+\frac{\mathop{\mathrm{cut}}(V_{i},V_{j})}{\omega(V_{j})}\right) =
\sum_{i=1}^{k}\frac{\mathop{\mathrm{cut}}(V_{i},V\backslash V_{i})}{ \omega(V_{i})}.$$

## Facts

If we want to cluster vertices of graph $G=(V,E)$ into $k$ clusters, we can apply the following recursive algorithm:

1. __Initialization.__ Compute the bipartition $\pi=\{V_{1},V_{2}\}$ of $V$. Set the counter $c=2$.

2. __Recursion.__ While $c<k$ repeat:

    1. Compute the bipartition of each subset of $V$.
    
    2. Among all $(c+1)$-partitions, choose the one with the smallest $\mathop{\mathrm{pcut}}(\pi_{c+1})$ or $\mathop{\mathrm{ncut}}(\pi_{c+1})$, respectively.
    
    3. Set $c=c+1$.

3. __Stop.__

There is no guarantee for optimality of this algorithm. Clearly, the optimal $k$-partiton may be a subpartition of one of the discarded partitions.
"""

# ╔═╡ f3eb8ad4-55a8-4e76-a581-07f9a7cc75f7


# ╔═╡ Cell order:
# ╠═77c3ba33-a12d-4e48-bbeb-727166771e63
# ╠═0cd92d78-1da6-435b-ad8e-60af29502e1d
# ╠═7ef4b0b3-9d0a-4356-bd16-14487f14e0ab
# ╟─cb61f761-b2ee-4818-9149-f85d95b76a1e
# ╟─5926c79a-6038-4dbc-b49a-5b05d02250a0
# ╟─a66a6e80-1f80-11eb-15f8-71ffaa5050e7
# ╟─d4cf12ce-12fb-4562-bbea-e9b85cba94c2
# ╠═87cca1d3-b944-4d74-bf4c-1562f34cf91a
# ╠═3d7b4a01-9a7d-4111-a45d-1d31f1ef6278
# ╠═3b5ed1fa-75b0-43fb-a918-4c2e8f52d039
# ╠═8f7d1f6c-ee4c-407f-80d9-d2abc2af949e
# ╠═b7f127d8-e8cd-4771-8091-c88f9ffe6b08
# ╠═3f5447bb-95cd-4520-b0b9-e85f33b323bd
# ╠═ecb2a46e-f9dc-4e57-812f-0a2788afb202
# ╠═86e5f00f-13ee-4867-b0cd-187f2c782a46
# ╠═cf26dbe9-3c4f-4c95-88ef-f197f360d759
# ╠═8fcd2440-4013-4b37-a747-2c0c7f6d4aeb
# ╠═0c8d3572-3436-446f-952a-b19d0370ab38
# ╠═1437aadc-1c26-4eff-bb22-9cf9c15e1f5f
# ╟─6529324b-49f7-4f01-b5cd-33668979b56b
# ╟─a9ba64a6-bd17-4b86-acb4-e7c7d8f87dc1
# ╠═6a6f5db5-7dd4-4a26-ac1a-ef3ca3a5d6be
# ╠═e11bc23d-03e2-4780-9488-93f0adfc2fcb
# ╠═7192f7d0-28b9-11eb-28dd-8db0019504ca
# ╠═e5890dca-368f-437e-bc45-6a395714dae2
# ╠═61a978a7-b2d2-43e7-b3f7-95bf4c4f54b9
# ╠═da7f85ed-73b0-44b9-a97f-4a8ed5e9a1d1
# ╠═fed7f592-be03-4baf-8ace-ff9ee78fee71
# ╟─9f31928e-c933-43bb-9a5c-bf55f29737cd
# ╟─1e6d4c72-f25d-4bf5-99d6-ae0446be25dd
# ╟─02035eed-6d49-4253-a755-68a999e7e90e
# ╠═c188b7cf-03e7-40b3-bf90-773978c7b410
# ╠═f60b51b6-79c3-4a3f-8bfb-8f198cdb92a8
# ╟─bfac7a55-0d82-4549-b999-7db3272eaa20
# ╠═80ed91ad-6f9a-4f12-a1ee-bc8dc6b58168
# ╠═2d8689c9-1b9e-4360-865c-23c2c8c5dcb9
# ╠═fa4ae713-a0b2-46ce-b5f7-558b26d82019
# ╠═760e8ca2-596f-435f-befb-140b3d204f7f
# ╠═75956488-76ab-4262-a32f-0c53e1b3dc17
# ╠═a55260dd-f05d-4a65-9fdf-cb3fd4638857
# ╠═ba164ff3-5ef4-4480-8af2-9dd9a99b1e8d
# ╠═46de28bc-73ed-495c-a2e0-f565c1ce5651
# ╟─e004e4bf-ca04-4b41-947b-9f1ecd058abb
# ╠═bb9d4d98-3aff-439e-b8a3-8347cf345839
# ╟─012fa54a-5171-4980-a0ba-474a22c25981
# ╠═f3eb8ad4-55a8-4e76-a581-07f9a7cc75f7
