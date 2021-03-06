### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# ╔═╡ 4e42f72c-9b68-4f32-99d3-04443e353d3c
begin
	import Pkg
	Pkg.activate(mktempdir())
    Pkg.add([
        Pkg.PackageSpec(name="PlutoUI"),
		Pkg.PackageSpec(name="Polynomials"),
		Pkg.PackageSpec(name="Plots")
    ])
end

# ╔═╡ 8d74e339-3df8-482a-a1c9-ce02369fa5e9
begin
	using PlutoUI, Random, Polynomials, Plots, LinearAlgebra
	plotly()
end

# ╔═╡ a6ce470f-9c2f-4dc8-9c92-878d1eaed316
TableOfContents(title="📚 Table of Contents", aside=true)

# ╔═╡ 68cc76f3-e124-4129-9339-98a21093bf1f
begin
	# Generate random points
	Random.seed!(125)
	n=6
	x=rand(n)
	y=rand(n)
	x₀=minimum(x)
	xₙ=maximum(x)
end

# ╔═╡ 8cb77244-7845-4e64-920e-11a024f880e6
begin
	# Functions to manipulate Vandermonde matrices
	import Base.getindex, Base.size
	struct Vandermonde{T} <: AbstractMatrix{T}
		c :: AbstractVector{T}
	end
	
	getindex(V::Vandermonde, i::Int, j::Int) = V.c[i]^(j-1)
	isassigned(V::Vandermonde, i::Int, j::Int) = isassigned(V.c, i)
	
	size(V::Vandermonde, r::Int) = (r==1 || r==2) ? length(V.c) :
	    throw(ArgumentError("Invalid dimension $r"))
	size(V::Vandermonde) = length(V.c), length(V.c)
	
	function Matrix(V::Vandermonde{T}) where T
		n=size(V, 1)
		M=Array{T}(undef,n, n)
		for i=1:n
			M[:,i] = V.c.^(i-1)
		end
		M
	end
	
end

# ╔═╡ 724b6f69-0f5c-4bcf-8363-42f308897070
md"""
# Interpolation Polynomials

We are given $n+1$ points in the plane

$$T_i=(x_i,y_i), \quad i=0,1,\ldots,n,\quad x_i\neq x_j.$$

## Standard basis

In the basis

$$1,x,x^2,x^3\ldots,x^n$$

through given points we can fit __interpolation polynomial__ $p_n(x)$,

$$p_n(x)={\displaystyle {\begin{aligned}a_{0}&+a_{1}x+a_{2}x^{2}+a_{3}x^{3}+\cdots +a_{n}x^{n}.\end{aligned}}}$$


Coefficients of this polynomial satisfy the system of (linear) equations 

$$p_n(x_i)=y_i, \quad i=0,\ldots,n,$$

or

$$\begin{pmatrix} 
1 & x_0 & x_0^2 & x_0^3 & \cdots & x_0^n \cr
1 & x_1 & x_1^2 & x_1^3 & \cdots & x_1^n \cr
\vdots & & & & \vdots \cr
1 & x_n & x_n^2 & x_n^3 & \cdots & x_n^n \cr
\end{pmatrix}
\begin{pmatrix}a_0\cr a_1 \cr \vdots \cr a_n\end{pmatrix}
=\begin{pmatrix} y_0 \cr y_1 \cr \vdots \cr y_n\end{pmatrix}.$$

The system matrix $A$ is a __Vandermonde matrix__. 
Its determinant is given by 

$$\mathop{\mathrm{det}}(A)= \prod_{0\leq j<i\leq n}(x_i-x_j).$$

Since all abscissas are different ($x_i\neq x_j$ for $i\neq j$), we have $\mathop{\mathrm{det}}(A)\neq 0$.
Therefore, the matrix $A$ is non-singular and the system has unique solution. We conclude that

> the interpolation polynomial is __unique__.
"""

# ╔═╡ 4f060d90-7b59-4447-ae59-856663acf3d7
A=Vandermonde(x)

# ╔═╡ 99d7f5e3-01c6-46e7-b16f-d42b46214738
a=A\y

# ╔═╡ 366d7a3e-ff6e-4bd4-a92d-935bd0d30c35
p=Polynomial(a)

# ╔═╡ f719cad8-cd2c-4013-9014-1caef03cc575
# Given points
scatter(x,y,label="Points")

# ╔═╡ f57091bd-a669-461f-ba5f-9ab7539a3c37
# Plot the polynomial 
plot!(p,label="Polynomial",xlims=(0,1),ylims=(-20,20))

# ╔═╡ 51668880-461a-45ce-82fe-204543677c75
begin
	# Plot the polynomial using our function
	xᵣ=range(x₀,stop=xₙ,length=100)
	pS=p.(xᵣ)
	plot(xᵣ,pS,label="Polynomial")
	scatter!(x,y,label="Points")
end

# ╔═╡ d29a09b2-7cc6-477b-82fd-0591f1f0ec8f
# Vandermonde matrix has large condition number
cond(A)

# ╔═╡ ea329cb4-5b4b-49e3-825a-eb19d54a4e91
md"""
Solution of the given system using standard Gaussian elimination requires $O(n^3)$ floating-point operations.
However, there are methods for solving Vandermonde systems using $O(n^2)$ operations.

The evaluation of the poynomial in each point requires 
$2n$ operations (Horner sheme).

Vandermonde matrices usually have big large condition number, so this way of computing ceofficients of the interpolation polynomial can be unstable.
For this reason we can also use other approached to computing and evaluation interpolation polynomials. 

## Lagrange interpolation polynomial

We define $n+1$ polynomials of degree $n$:

$$
L_j(x)=\prod_{\displaystyle {i=0}\atop {\displaystyle i\neq j}}^n \frac{x-x_i}{x_j-x_i}$$

and use it as basis. It holds

$$
L_j(x_i)=\begin{cases}0, \quad i\neq j \\ 1,\quad i=j \end{cases}$$

so 

$$
p_n(x)=y_0\, L_0(x)+y_1 \, L_1(x)+\cdots + y_n\,  L_n(x).$$

Computation of denominators initially requires $O(n^2)$ operations, but afterwards each $p_n(x)$ is computed using  $O(n)$ operations (__explain how!__). 

The following implementation of the algorithm is simple, but not optimal.
"""

# ╔═╡ 4195f039-de9e-4046-89c2-45e328b30478
L(t)=sum(y.*[prod(t .-x[[1:j-1;j+1:end]])/prod(x[j].-x[[1:j-1;j+1:end]]) 
        for j=1:n])

# ╔═╡ ef4cdc9e-7e9a-4b1a-b25c-c8e4266f8d6a
begin
	pL=Array{Float64}(undef,length(xᵣ))
	for i=1:length(xᵣ)
	    pL[i]=L(xᵣ[i])
	end
end

# ╔═╡ 70cadf28-cb7c-4f44-836d-df7d7af666c2
begin
	plot(xᵣ,pL,label="Polynomial")
	scatter!(x,y,label="Points")
end

# ╔═╡ 0d08175e-43ab-41d4-b278-697154aa1966
norm(pS-pL,Inf)

# ╔═╡ f4766fb0-d630-4ef0-80e2-dedb057596bb
norm(abs.((pS-pL)./pL),Inf)

# ╔═╡ 65bf04af-9c28-4d16-b71d-b6d15a9371e3
md"""
## Newton interpolation polynomial

Here we use the basis

$$
\begin{aligned}
& 1,\\
& x-x_0,\\
& (x-x_0)(x-x_1),\\
& (x-x_0)(x-x_1)(x-x_2),\\
& \qquad \vdots\\
& (x-x_0)(x-x_1)\cdots (x-x_{n-1})
\end{aligned}$$

so the polynomial is given by

$$
p_n(x)=c_0 + c_1(x-x_0)+c_2(x-x_0)(x-x_1)+\cdots +c_n(x-x_0)(x-x_1)\cdots (x-x_{n-1}).$$

The coefficients $c_0,c_1,\ldots, c_n$ are solution of the __triangular__ system of linear equations $Lc=y$,

$$
\begin{pmatrix} 
1 & 0 & 0 & 0 & \cdots & 0 \\
1 & x_1-x_0 & 0 & 0 & \cdots & 0 \\
1 & x_2-x_0 & (x_2-x_0)(x_2-x_1) & 0 & \cdots & 0 \\
\vdots & & & & \vdots \\
1 & x_n-x_0 & (x_n-x_0)(x_n-x_1) & (x_n-x_0)(x_n-x_1)(x_n-x_2) & \cdots & (x_n-x_0)\cdots (x_n-x_{n-1}) \\
\end{pmatrix}
\begin{pmatrix}c_0\\ c_1 \\ c_2 \\\vdots \\ a_n\end{pmatrix}
=\begin{pmatrix} y_0 \\ y_1 \\ y_2 \\ \vdots \\ y_n\end{pmatrix}.$$

Forming the lower triangular matrix $L$ requires $O(n^2)$ operations. Compuation of the coefficients  $c_0,\ldots,c_n$ requires $O(n^2)$ operations (solution of lower triangular system) and this solution is __stable__.

For evaluation of $p_n(x)$ we use an algorithm similar to Horner scheme. 
"""

# ╔═╡ d5755a6e-a155-4090-af0e-1076501397fa
# Computing coefficients c
function Newton(x,y)
    n=length(x)
    L=zeros(n,n)
    L[:,1]=ones(n)
    for i=2:n
        for j=2:i
            L[i,j]=prod([x[i]-x[k] for k=1:j-1])
        end
    end
    c=L\y
end  

# ╔═╡ 4976320f-d288-49a3-a1cf-44d2d223b488
c=Newton(x,y)

# ╔═╡ b85666e1-49c3-4f2d-b55b-9c48df27159a
# Evaluating Newton's polynomial through abscissas x and 
# coefficients c in point t 
function evalnewton(c,x,t::Number)
    p=c[end]
    for i=length(c)-1:-1:1
        p=p*(t-x[i])+c[i]
    end
    p
end

# ╔═╡ 3f6a6282-3363-4a51-a5a8-bacc3d7e8880
begin
	pN=Array{Float64}(undef,length(xᵣ))
	for i=1:length(xᵣ)
	    pN[i]=evalnewton(c,x,xᵣ[i])
	end
end

# ╔═╡ 2060e3ad-c95f-4cd6-99c6-673da7dd64ac
begin
	plot(xᵣ,pN,label="Polynomial")
	scatter!(x,y,label="Points")
end

# ╔═╡ 0e0637af-4d77-442b-af67-ae3b7fe9075a
norm(abs.((pS-pN)./pN),Inf)

# ╔═╡ da067ec6-5085-467c-ae0d-143e1ba1e922
norm(abs.((pL-pN)./pN),Inf)

# ╔═╡ 477e1427-cd04-40f0-9b18-00ce68b1b837
md"""
We see that values `pN` and `pL` are closer to each other than to `pS` so we conclude that they are more accurate.
"""

# ╔═╡ Cell order:
# ╠═4e42f72c-9b68-4f32-99d3-04443e353d3c
# ╠═8d74e339-3df8-482a-a1c9-ce02369fa5e9
# ╠═a6ce470f-9c2f-4dc8-9c92-878d1eaed316
# ╟─724b6f69-0f5c-4bcf-8363-42f308897070
# ╠═68cc76f3-e124-4129-9339-98a21093bf1f
# ╠═8cb77244-7845-4e64-920e-11a024f880e6
# ╠═4f060d90-7b59-4447-ae59-856663acf3d7
# ╠═99d7f5e3-01c6-46e7-b16f-d42b46214738
# ╠═366d7a3e-ff6e-4bd4-a92d-935bd0d30c35
# ╠═f719cad8-cd2c-4013-9014-1caef03cc575
# ╠═f57091bd-a669-461f-ba5f-9ab7539a3c37
# ╠═51668880-461a-45ce-82fe-204543677c75
# ╠═d29a09b2-7cc6-477b-82fd-0591f1f0ec8f
# ╟─ea329cb4-5b4b-49e3-825a-eb19d54a4e91
# ╠═4195f039-de9e-4046-89c2-45e328b30478
# ╠═ef4cdc9e-7e9a-4b1a-b25c-c8e4266f8d6a
# ╠═70cadf28-cb7c-4f44-836d-df7d7af666c2
# ╠═0d08175e-43ab-41d4-b278-697154aa1966
# ╠═f4766fb0-d630-4ef0-80e2-dedb057596bb
# ╟─65bf04af-9c28-4d16-b71d-b6d15a9371e3
# ╠═d5755a6e-a155-4090-af0e-1076501397fa
# ╠═4976320f-d288-49a3-a1cf-44d2d223b488
# ╠═b85666e1-49c3-4f2d-b55b-9c48df27159a
# ╠═3f6a6282-3363-4a51-a5a8-bacc3d7e8880
# ╠═2060e3ad-c95f-4cd6-99c6-673da7dd64ac
# ╠═0e0637af-4d77-442b-af67-ae3b7fe9075a
# ╠═da067ec6-5085-467c-ae0d-143e1ba1e922
# ╟─477e1427-cd04-40f0-9b18-00ce68b1b837
