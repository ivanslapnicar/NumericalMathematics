### A Pluto.jl notebook ###
# v0.12.10

using Markdown
using InteractiveUtils

# ╔═╡ 6a0417bf-00aa-47db-bf22-5842f1cc736e
begin
	using Polynomials
	using Plots
end

# ╔═╡ df89a063-e647-4bfe-91eb-167be078ac0e
md"""
# Interpolating Functions


Consider a function $f(x)$ on an interval $[a,b]$.

Wer choose $n+1$ points $x_i,\ i=0,\ldots, n$, in the interval $[a,b]$ such that $x_i\neq x_j$
and construct interpolation polynomial through the points $T_i=(x_i,f(x_i))$.

__Theorem.__ For every $x\in[a,b]$ the __error bound__ holds (we assume that $f$ is $n+1$ times differentiable:) 

$$\begin{aligned}
f(x)-p_n(x)&=\frac{\omega(x)}{(n+1)!} \,f^{(n+1)}(\xi), \cr
\omega(x)&=\prod_{k=0}^n (x-x_k)=(x-x_0)(x-x_1)\cdots (x-x_n),\quad  \xi \in (a,b).
\end{aligned}$$


_Proof._ (See [Numerička matematika, p. 23](http://www.mathos.unios.hr/pim/Materijali/Num.pdf).) 

For $x=x_i$, the statement is obvious. Let $x\neq x_i$. We define auxilliary function

$$g(y)= f(y)-p_n(y)-k\omega(y),$$

where the constant $k$ is chosen such that $g(x)=0$.
In this manner $g$ has at least $n+2$ zeros, 
$g(x)=g(x_0)=g(x_1)=\cdots=g(x_n)=0$. 
By Rolle theorem, the derivative $g'$ nas at least $n+1$ zeros, $g''$ has at least $n$ zeros, etc. 
The function $g^{(n+1)}$ has at least one zero $\xi\in(a,b)$.

It holds $p_n^{(n+1)}(y)=0$ and $\omega^{(n+1)}(y)=(n+1)!$ (the leading coefficient of $\omega(x)$ is $1$). Inserting gives

$$0=g^{(n+1)}(\xi)=f^{(n+1)}(\xi)-k(n+1)!$$ 

so $k=\displaystyle\frac{f^{(n+1)}(\xi)}{(n+1)!}$, which proves the theorem.


## Example

Let

$$
f(x)=\sin(x), \quad x\in[0,\pi].$$
"""

# ╔═╡ c88bd3a2-a46d-4931-b4cf-0b940696aa89
begin
	# These function are used to manipulate Vandermonde matrices
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

# ╔═╡ 3a599b53-077b-47e7-a1cf-17e15da6aa1f
begin
	n=7
	a=0
	b=pi
	x=range(a,stop=b,length=n)
	f(x)=sin(x)
	y=f.(x)
end

# ╔═╡ de651f77-3953-4ada-ac74-48bd58d147c1
A=Vandermonde(x)

# ╔═╡ 7f63d94f-3792-4369-b8f4-edd3a5b0cd78
c=A\y

# ╔═╡ 7c1650eb-bf62-40af-b085-1b81da0f47bf
p=Polynomial(c)

# ╔═╡ af46e285-715a-4b5c-8745-d17d61806dbf
# Plot points on f(x)
scatter(x,y,label="Points")

# ╔═╡ 6a1c5393-f4c9-4a45-8a11-0b4c884406b7
begin
	# Plot f(x) and interpolating polinomial
	x₀=range(a,stop=b,length=100)
	p₀=p.(x₀)
	F₀=f.(x₀)
	plot!(x₀,[p₀ F₀],label=["Polynomial" "Function"])
end

# ╔═╡ 9e80a29c-ea79-489f-8383-60fab341f5be
begin
	# Maximal absolute and relative errors
	using LinearAlgebra
	norm(p₀[2:end-1]-F₀[2:end-1],Inf), 
	norm((p₀[2:end-1]-F₀[2:end-1])./F₀[2:end-1],Inf)
end

# ╔═╡ 58277216-ffc5-47e7-a8c6-b181e54ec870
md"""
## Chebyshev points

__Chebyshev polynomial__ $T_n(x)$ is the polynomial of order $n$ given by the formula

$$
T_n(x)=\cos(n\, \arccos x), \quad n=0,1,2,\ldots$$

Trigonometric addition formula 

$$
\cos(n+1)\varphi+\cos(n-1)\varphi=2\cos\varphi \cos n\varphi$$

with $\varphi=\arccos x$ implies recursive formulas:

$$
\begin{aligned}
T_0(x)&=1,\cr 
T_1(x)&=x,\cr 
T_{n+1}(x)&=2\,x\,T_n(x)-T_{n-1}(x), \quad n=1,2,\ldots
\end{aligned}$$

For example, 

$$
T_2(x)=2x^2-1,\quad T_3(x)=4x^3-3x, \ldots$$

The zeros of $T_n(x)$ are 

$$
x_k=\cos \bigg(\frac{2k-1}{2n}\pi\bigg), \quad k=1,2,\ldots,n,$$

and they all lie interval $[-1,1]$. In the points 

$$
\xi_k=\cos \bigg(\frac{k\pi}{n}\bigg), \quad k=0,1,2,\ldots,n,$$

$T_n(x)$ attains alternately local maxima and minima $1$ and $-1$,  respectively, on the interval $[-1,1]$.  

On the interval $[-1,1]$, the polynomial $T_n(x)$ attains values in the interval $[-1,1]$.



### Example
"""

# ╔═╡ 45533a9a-caa8-401b-9c22-40e9afb46bbe
T(n,x)=cos.(n*acos.(x))

# ╔═╡ 79a264ae-3635-466a-9a2c-4ffe612f417d
x₁=range(-1,stop=1,length=100)

# ╔═╡ 89136037-71ac-4dbe-9621-90f36c7b2702
y₁=T(10,x₁)

# ╔═╡ 62239acf-92c3-4b68-ae54-5e6a1dd4fb53
plot(x₁,y₁,label="Chebyshev polynomial")

# ╔═╡ 4f65fc49-4d4e-43d7-9775-d53e1e766805
xₙ=[cos((2*k-1)*pi/(2*10)) for k=10:-1:1]

# ╔═╡ 16f339d9-01ad-4409-8fec-4ff797e8fd02
yₙ=T(10,xₙ)

# ╔═╡ f2a33009-575a-4912-9c68-223e6862407e
scatter!(xₙ,yₙ,label="Zeros")

# ╔═╡ 2ee3a5fe-8da3-4a37-b61b-63fa0282d287
md"""
### Norms of functions

For functions 

$$f,g:[a,b]\to \mathbb{R}$$

define __scalar product__ (or __dot product__) as

$$
(f,g)=\int_a^b f(x)g(x)\, dx$$

and the __weighted scalar product__ with __weight__ $\omega(x)>0$ as

$$
(f,g)_\omega=\int_a^b f(x)g(x)\omega(x)\, dx.$$

Functions $f$ and $g$ are __orthogonal__ if $(f,g)=0$ or if $(f,g)_\omega=0$.

The following three __norms__ are natural generalizations of the corresponding vector norms:

$$
\begin{aligned}
\|f\|_2&=\sqrt{(f,f)}=\sqrt{\int_a^b f^2(x)\, dx} \cr
\|f\|_1 &= \int_a^b \big|f(x)\big|\, dx \cr
\|f\|_\infty&=\max_{x\in[a,b]} \big|f(x)\big|
\end{aligned}$$

"""

# ╔═╡ 5ab52c00-19e5-11eb-1f4a-7199bdd0be95
md"""
We have the following important theorem:

__Theorem__. Of all polynomials of degree at most $n$ with leading coefficient equal to $1$, the polynomial $\displaystyle\frac{1}{2^{n-1}}T_n(x)$ has the smallest norm
$\|\cdot\|_\infty$ on the interval $[-1,1]$ and this norm is exactly $\displaystyle\frac{1}{2^{n-1}}$.

_Proof_ : (See [Numerička matematika, p. 95](http://www.mathos.unios.hr/pim/Materijali/Num.pdf).)

The coefficient of the power $x^n$ of the Chebyshev polynomial $T_n(x)$ equals $2^{n-1}$.
Therefore, coefficient of the power $x^n$ of the polynomial $\displaystyle\frac{1}{2^{n-1}}T_n(x)$ is $1$. This and the properties of Chebishev polanomials imply

$$\left| \displaystyle\frac{1}{2^{n-1}}T_n(x) \right| \leq \displaystyle\frac{1}{2^{n-1}}$$

so

$$\left\| \displaystyle\frac{1}{2^{n-1}}T_n(x) \right\|_\infty = \displaystyle\frac{1}{2^{n-1}}.$$

Assume that the polynomial 

$$p_n(x)=x^n + \alpha_{n-1}x^{n-1} +\alpha_{n-2}x^{n-2}+\cdots \alpha_1 x+\alpha_0$$

satisfies $| p_n(x)|< \displaystyle\frac{1}{2^{n-1}}$ for each $x\in[-1,1]$. Let $\xi_0,\xi_1,\ldots,\xi_n$ be the points where $T_n(x)$ attains extremal values $-1$ or $1$. Then, we must have

$$
\begin{aligned}
p_n(\xi_0) & < \frac{1}{2^{n-1}}T_n(\xi_0) = \frac{1}{2^{n-1}} \cr
p_n(\xi_1) & > \frac{1}{2^{n-1}}T_n(\xi_1) = -\frac{1}{2^{n-1}} \cr
p_n(\xi_2) & < \frac{1}{2^{n-1}}T_n(\xi_2) = \frac{1}{2^{n-1}} \cr
& \vdots
\end{aligned}$$

This imples that the polynomial $p_n(x)-\frac{1}{2^{n-1}}T_n(x)$ of degree at most $n-1$ changes its sign $n$ times, which is impossibe. This completes the proof.


We conclude that the polynomial approximation (1) on the interval $[-1,1]$ will be the best in $\infty$-norm if we choose

$$\omega(x)=\frac{1}{2^{n}} T_{n+1}(x).$$ 

Thus, on the interval $[a,b]$ for interpolation points   
$x_0,x_1,\ldots,x_n$ we shall choose zeros of $T_{n+1}(x)$ mapped to the interval $[a,b]$.
"""

# ╔═╡ 7818f9c0-19e5-11eb-0d31-7d39ee46d05a
md"""
### Change of interval

The transformation 

$$
\gamma :[a,b]\to [-1,1],\quad \gamma(x)=\frac{2x}{b-a}-\frac{a+b}{b-a}$$

maps the system of orthogonal functions $\Phi_i$ defined on interval $[-1,1]$ to the system of orthogonal functions

$$
\Psi_i(x)=\Phi_i(\gamma(x))$$

on interval $[a,b]$.

We need the inverse transformation:

$$
x=\frac{a+b}{2}+\frac{b-a}{2}\gamma(x).$$
"""

# ╔═╡ 9bc340ed-75a2-4e50-b209-8075508d0611
# Interpolate sine using zeros of T(n,x)
xₜ=(a+b)/2 .+(b-a)/2*[cos((2*k-1)*pi/(2*n)) for k=n:-1:1]

# ╔═╡ 79d7ee5b-4dcc-47c4-b858-f6b6afb8b66a
yₜ=f.(xₜ)

# ╔═╡ e6ea8fd8-6fbf-4fc7-ac6a-7354b1e0be85
begin
	Aₜ=Vandermonde(xₜ)
	cₜ=Aₜ\yₜ
	pₜ=Polynomial(cₜ)
end

# ╔═╡ b2c12a8c-374a-4222-8e37-83998743988a
begin
	# x₂=range(a,stop=b,length=100)
	p₂=pₜ.(x₀)
	F₂=f.(x₀)
	scatter(xₜ,yₜ,label="Chebyshev points")
	plot!(x₀,[p₂ F₂],label=["Polynomial" "Function"])
end

# ╔═╡ e147f5fa-6826-4b21-80a4-1fbd14593c8c
# Maximal absolute and relative erros
norm(p₂[2:end-1]-F₂[2:end-1],Inf), 
norm((p₂[2:end-1]-F₂[2:end-1])./F₂[2:end-1],Inf)

# ╔═╡ 797d58e0-6aba-4171-82b2-574c08aae214
md"""
Let us plot true errors in both cases:
"""

# ╔═╡ af090d7e-155d-11eb-1b3b-77d83ad0b048
# Equidistant points
plot(x₀,p₀-F₀)

# ╔═╡ bba5f710-155d-11eb-1f00-892b62a32cb4
# Chebyshev points
plot(x₀,p₂-F₂)

# ╔═╡ 7aa65294-2803-4bca-aa69-ab11ed6e43a3
md"""
We se that the errors attained using Chebyshev points are smaller.

__Remark.__ For the sake of simplicity, we have used the least accurate method of computing interpolation polynomial.

### Example

Let us look at another interesting example (see [Numerička matematika, p. 24](http://www.mathos.unios.hr/pim/Materijali/Num.pdf)). Let us interpolate the function 

$$
f(x)=1-|x-1|,\quad x\in[0,2]$$.
 
using polynomials of degree 10. Large deviations of interpolation polynomial from the function on the edges are called __Runge's phenomenon__. The phenomenon dissapears when the function is interpolated using Chebyshev points.
"""

# ╔═╡ c6b0381f-dd60-42f3-8f6a-df8771bbf3c6
begin
	n₃=11
	a₃=0
	b₃=2
	f₃(x)=1 .-abs.(x .-1)
	# Equidistant points
	x₃=range(a₃,stop=b₃,length=n₃)
	y₃=f₃(x₃)
	A₃=Vandermonde(x₃)
	c₃=A₃\y₃
	p₃=Polynomial(c₃)
	xp₃=range(a₃,stop=b₃,length=100)
	pf₃=p₃.(xp₃)
	F₃=f₃(xp₃)
	scatter(x₃,y₃,label="Points")
	plot!(xp₃, [pf₃ F₃], label=["Polynomial" "Function"])
end

# ╔═╡ bcf8a882-eca1-4739-9789-43a9de96cb94
begin
	# Chebyshev points
	xt₃=(a₃+b₃)/2 .+(b₃-a₃)/2*[cos((2*k-1)*pi/(2*n₃)) for k=n₃+1:-1:1]
	yt₃=f₃(xt₃)
	At₃=Vandermonde(xt₃)
	ct₃=At₃\yt₃
	pc₃=Polynomial(ct₃)
	pCheb=pc₃.(xp₃)
	scatter(xt₃,yt₃,label="Points")
	plot!(xp₃,[pCheb F₃],label=["Polynomial" "Function"])
end

# ╔═╡ 3760da30-1560-11eb-083e-575c20a74101
xt₃

# ╔═╡ Cell order:
# ╟─df89a063-e647-4bfe-91eb-167be078ac0e
# ╠═6a0417bf-00aa-47db-bf22-5842f1cc736e
# ╠═c88bd3a2-a46d-4931-b4cf-0b940696aa89
# ╠═3a599b53-077b-47e7-a1cf-17e15da6aa1f
# ╠═de651f77-3953-4ada-ac74-48bd58d147c1
# ╠═7f63d94f-3792-4369-b8f4-edd3a5b0cd78
# ╠═7c1650eb-bf62-40af-b085-1b81da0f47bf
# ╠═af46e285-715a-4b5c-8745-d17d61806dbf
# ╠═6a1c5393-f4c9-4a45-8a11-0b4c884406b7
# ╠═9e80a29c-ea79-489f-8383-60fab341f5be
# ╟─58277216-ffc5-47e7-a8c6-b181e54ec870
# ╠═45533a9a-caa8-401b-9c22-40e9afb46bbe
# ╠═79a264ae-3635-466a-9a2c-4ffe612f417d
# ╠═89136037-71ac-4dbe-9621-90f36c7b2702
# ╠═62239acf-92c3-4b68-ae54-5e6a1dd4fb53
# ╠═4f65fc49-4d4e-43d7-9775-d53e1e766805
# ╠═16f339d9-01ad-4409-8fec-4ff797e8fd02
# ╠═f2a33009-575a-4912-9c68-223e6862407e
# ╟─2ee3a5fe-8da3-4a37-b61b-63fa0282d287
# ╟─5ab52c00-19e5-11eb-1f4a-7199bdd0be95
# ╟─7818f9c0-19e5-11eb-0d31-7d39ee46d05a
# ╠═9bc340ed-75a2-4e50-b209-8075508d0611
# ╠═79d7ee5b-4dcc-47c4-b858-f6b6afb8b66a
# ╠═e6ea8fd8-6fbf-4fc7-ac6a-7354b1e0be85
# ╠═b2c12a8c-374a-4222-8e37-83998743988a
# ╠═e147f5fa-6826-4b21-80a4-1fbd14593c8c
# ╟─797d58e0-6aba-4171-82b2-574c08aae214
# ╠═af090d7e-155d-11eb-1b3b-77d83ad0b048
# ╠═bba5f710-155d-11eb-1f00-892b62a32cb4
# ╟─7aa65294-2803-4bca-aa69-ab11ed6e43a3
# ╠═c6b0381f-dd60-42f3-8f6a-df8771bbf3c6
# ╠═bcf8a882-eca1-4739-9789-43a9de96cb94
# ╠═3760da30-1560-11eb-083e-575c20a74101
