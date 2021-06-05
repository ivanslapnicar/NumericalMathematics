### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# ╔═╡ b7cc077b-960b-41ad-bb20-51bd0cfbd4e6
begin
	import Pkg
	Pkg.activate(mktempdir())
    Pkg.add([
        Pkg.PackageSpec(name="PlutoUI"),
		Pkg.PackageSpec(name="Plots")
    ])
end

# ╔═╡ 881ed68b-ac65-434a-9d62-df963fb033b1
begin
	using PlutoUI, Random, LinearAlgebra, Plots
	Random.seed!(123)
	plotly()
end

# ╔═╡ 07eaaf8d-3326-49b0-9413-402497753ec2
TableOfContents(title="📚 Table of Contents", aside=true)

# ╔═╡ e0f03ec0-f136-411f-a270-bcc85fde0579
md"""
# Applications of QR Factorization


Orthogonal matrices have two important properties:

$Q^{-1}=Q^T$

$\|Qx\|_2=\|x\|_2,\quad \forall x$

The first property follows from the definition of orthogonal matrix, since $Q^TQ=I$, 
and the second property follows from 

$$\|Qx\|_2^2=(Qx)^T(Qx)=x^TQ^TQx=x^Tx=\|x\|_2^2.$$

## Systems of linear equations

QR factorization can be used to solve system of linear equations $Ax=b$: premultiplying $QRx=b$ by $Q^T$ yields $Rx=Q^Tb$, so it remains to solve triangular system.

With respect to the solution using Gaussian elimination, it holds:

* the number of floating-point operations doubles,
* the solution is somewhat more accurate, and
* there is no element growth (pivoting is not necessary).

### Example
"""

# ╔═╡ d0aa1de2-f2e5-4c59-b7ce-a0022bf6dbef
begin
	n=10
	A=rand(n,n)
	b=rand(n)
	Q,R=qr(A)
	c=transpose(Q)*b
	# Triangular system
	x=R\c
end

# ╔═╡ fbfa49ba-8c53-4fdf-ae26-d9b678fcce2e
# Residual
A*x-b

# ╔═╡ 45923270-4d63-4eaf-8d99-f9e9baab558e
md"""
## Least squares problem

Programs or functions which solve least squares problems usually use QR factorization. We have

$$
\|Ax-b\|^2_2=\|QRx-b\|_2^2=\|Q(Rx-Q^Tb)\|_2=\|Rx-Q^Tb\|_2^2.$$

Let

$$
R=\begin{bmatrix}R_0 \\ 0\end{bmatrix},\quad Q^Tb =\begin{bmatrix}c\\ d \end{bmatrix}.$$

Then,

$$
\|Rx-Q^T b\|_2^2 = \| R_0x-c\|_2^2+\|d\|_2^2,$$

so the solution of the triangular system

$$
R_0x=c$$
    
is the solution of the least squares problem. Let us solve the problem from the regression notebook:

### Small example
"""

# ╔═╡ 2e4f479e-9a3a-4495-abb9-52a3203a35c4
begin
	y₁=[1,3,2,4,3]
	A₁=transpose([1 2 3 6 7;1 1 1 1 1])
end

# ╔═╡ 743b75c7-d6a0-4fa0-bd99-6412378453e2
#?qr  # Observe the structure of the solution

# ╔═╡ 04bb12df-5d6a-460f-a381-23058dcfe0b7
begin
	F₁=qr(A₁)
	c₁=transpose(Matrix(F₁.Q))*y₁
	x₁=F₁.R\c₁
end

# ╔═╡ c4762bb2-56fd-4e2a-9255-4922f41addfe
# Built-in function
A₁\y₁

# ╔═╡ a565c881-4f9a-4ab4-92a1-04e6fe0a8315
md"
### Random example
"

# ╔═╡ 881daa94-e922-4b76-bca3-6b255049ba67
begin
	# Bigger example
	m₂=8
	n₂=5
	A₂=randn(m₂,n₂)
	b₂=randn(m₂)
	F₂=qr(A₂)
end

# ╔═╡ 16461985-6cb3-4ddb-9eab-bd13917b8b69
F₂.Q

# ╔═╡ d01748e8-169d-4d39-bf79-90c0c0d2df56
# Generators
F₂.factors

# ╔═╡ ab0d603e-4e1f-4ec9-a247-c9c78c9e6488
begin
	# Solution
	c₂=transpose(Matrix(F₂.Q))*b₂
	x₂=F₂.R\c₂
end

# ╔═╡ f97a14a2-1e13-4a2d-81a6-c848666dbbfe
# Built-in function
A₂\b₂

# ╔═╡ 0ad1e3dc-e5d4-416a-b3e0-21fd710ed241
md"""
## Numerical  "orthogonalization" of polynomials

Numerical orthogonalization of powers of vectors produces orthogonal polynomials.
"""

# ╔═╡ 24ac4e43-b99d-439e-aa6b-ba3bc494b8dd
begin
	# Standard basis
	xₒ=range(-1,stop=1,length=101)
	# Quasi Vandermonde matrix
	V=[xₒ.^0 xₒ.^1 xₒ.^2 xₒ.^3 xₒ.^4 xₒ.^5]
end

# ╔═╡ 9006213a-325f-4296-a0ec-d5c02b4ec5e8
plot(xₒ,V,title="Standard basis",legend=:bottomright,
	label=["1" "x" "x^2" "x^3" "x^4" "x^5"])

# ╔═╡ 2b294375-d2de-4aeb-bf01-7646bc9a935a
md"
### Legendre polynomials
"

# ╔═╡ f157bb76-8445-4994-b08e-260b2a51ad51
begin
	# Orthogonalization with the weight function ω(x)=1 produces Legendre polynomials.
	Fₒ=qr(V)
	Qₒ=Matrix(Fₒ.Q)*sign.(Diagonal(Fₒ.R))
end

# ╔═╡ 68d56539-22ab-4a89-ad51-fb3c2f37ee8d
plot(xₒ,Qₒ,title="Legendre polynomials",label=["L₀" "L₁" "L₂" "L₃" "L₄" "L₅"])

# ╔═╡ 4b343213-f8d2-45d3-bae8-999dee675089
md"""
### Chebyshev polynomials

Given normalized vectors are values of scaled Legendre polynomials from notebook [L12 Orthogonal Polynomials.jl](https://ivanslapnicar.github.io/NumericalMathematics/L12%20Orthogonal%20Polynomials.jl.html).

In order to obtain Chebyshev polynomials, we need to use weight function $\omega(x)=\displaystyle\frac{1}{\sqrt{1-x^2}}$ and modify function `GramSchmidtQR()` from notebook [L15 QR Factorization.jl](https://ivanslapnicar.github.io/NumericalMathematics/L15%20QR%20Factorization.jl.html) such that it computes __weighted scalar products__. The obtained normalized vectors are values of scaled Chebyshev polynomials.
"""

# ╔═╡ 830acec9-74f1-4b02-a5f7-b65b05f457ef
function WeightedGramSchmidtQR(A::Array,ω::Vector)
    m,n=size(A)
    R=zeros(n,n)
    Q=Array{Float64,2}(undef,m,n)
    R[1,1]=norm(A[:,1])
    Q[:,1]=A[:,1]/R[1,1]
    for k=2:n
        for i=1:k-1
            R[i,k]=Q[:,i]⋅(A[:,k].*ω)/(Q[:,i]⋅(Q[:,i].*ω))
        end
        t=A[:,k]-sum([R[i,k]*Q[:,i] for i=1:k-1])
        R[k,k]=norm(t)
        Q[:,k]=t/R[k,k]
    end
    return Q,R
end

# ╔═╡ e0bf68c5-1e30-4621-871b-47eb71bfd2bb
begin
	x₃=range(-0.99,stop=0.99,length=101)
	ω=1 ./(sqrt.(1.0.-x₃.^2))
	# Quasi Vandermonde matrix
	V₃=[x₃.^0 x₃.^1 x₃.^2 x₃.^3 x₃.^4 x₃.^5]
end

# ╔═╡ be2197e5-565e-4b13-acb9-a245f2f3780e
begin
	Q₃,R₃=WeightedGramSchmidtQR(V₃,ω)
	Q₃=Q₃*sign.(Diagonal(R₃))
end

# ╔═╡ 957a675e-a12a-47bf-b196-bd4f3431849a
plot(x₃,Q₃,title="Chebyshev polynomials",label=["T₀" "T₁" "T₂" "T₃" "T₄" "T₅"])

# ╔═╡ 656aa3d6-0ebc-4ca5-95d5-071af0f5c2c0
md"""
__Problem.__ Normalize columns of the matrix $Q$ such that the vectors attain all values from the interval $[-1,1]$.  
"""

# ╔═╡ Cell order:
# ╠═b7cc077b-960b-41ad-bb20-51bd0cfbd4e6
# ╠═881ed68b-ac65-434a-9d62-df963fb033b1
# ╠═07eaaf8d-3326-49b0-9413-402497753ec2
# ╟─e0f03ec0-f136-411f-a270-bcc85fde0579
# ╠═d0aa1de2-f2e5-4c59-b7ce-a0022bf6dbef
# ╠═fbfa49ba-8c53-4fdf-ae26-d9b678fcce2e
# ╟─45923270-4d63-4eaf-8d99-f9e9baab558e
# ╠═2e4f479e-9a3a-4495-abb9-52a3203a35c4
# ╠═743b75c7-d6a0-4fa0-bd99-6412378453e2
# ╠═04bb12df-5d6a-460f-a381-23058dcfe0b7
# ╠═c4762bb2-56fd-4e2a-9255-4922f41addfe
# ╠═a565c881-4f9a-4ab4-92a1-04e6fe0a8315
# ╠═881daa94-e922-4b76-bca3-6b255049ba67
# ╠═16461985-6cb3-4ddb-9eab-bd13917b8b69
# ╠═d01748e8-169d-4d39-bf79-90c0c0d2df56
# ╠═ab0d603e-4e1f-4ec9-a247-c9c78c9e6488
# ╠═f97a14a2-1e13-4a2d-81a6-c848666dbbfe
# ╟─0ad1e3dc-e5d4-416a-b3e0-21fd710ed241
# ╠═24ac4e43-b99d-439e-aa6b-ba3bc494b8dd
# ╠═9006213a-325f-4296-a0ec-d5c02b4ec5e8
# ╟─2b294375-d2de-4aeb-bf01-7646bc9a935a
# ╠═f157bb76-8445-4994-b08e-260b2a51ad51
# ╠═68d56539-22ab-4a89-ad51-fb3c2f37ee8d
# ╟─4b343213-f8d2-45d3-bae8-999dee675089
# ╠═830acec9-74f1-4b02-a5f7-b65b05f457ef
# ╠═e0bf68c5-1e30-4621-871b-47eb71bfd2bb
# ╠═be2197e5-565e-4b13-acb9-a245f2f3780e
# ╠═957a675e-a12a-47bf-b196-bd4f3431849a
# ╟─656aa3d6-0ebc-4ca5-95d5-071af0f5c2c0
