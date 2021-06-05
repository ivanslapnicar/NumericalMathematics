### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# ╔═╡ d0d29fe8-6542-4968-a695-03371ba85543
begin
	import Pkg
	Pkg.activate(mktempdir())
    Pkg.add([
        Pkg.PackageSpec(name="PlutoUI"),
		Pkg.PackageSpec(name="Plots"),
		Pkg.PackageSpec(name="ForwardDiff")
    ])
end

# ╔═╡ 3a20c826-03cd-4c06-a0c6-7dc657067feb
begin
	using PlutoUI, Plots, ForwardDiff
	plotly()
end

# ╔═╡ 7bab2ccd-e6d7-479c-abef-d9ed8eefaf2a
TableOfContents(title="📚 Table of Contents", aside=true)

# ╔═╡ 1198a836-bb14-4b8e-9f30-160097bc4507
md"""
# Nonlinear Equations

__Problem.__ Find zeros of the function $f(x)$ in the closed interval $[a,b]$, that is, solve the equation

$$
f(x)=0, \quad x\in[a,b]. \qquad\qquad (1)$$

It holds:

If $f$ is __continuous__ on $[a,b]$ and if $f(a)\cdot f(b)<0$, then there exists at least one point $\xi\in(a,b)$ such

$$
f(\xi)=0.$$

If, additionally, $f'(x)\neq 0$ za $x\in(a,b)$, then $\xi$ is __unique__.

Thus, the equation (1) can be solved in two steps:

1. Find the interval $[a,b]$ in which function $f$ has unique zero $\xi$,
2. Approximate $\xi$ with accuracy given in advance.

We shall desribe four methods:

1. bisection,
2. simple iterations,
3. Newton's method (Tangent method), and 
4. Secant method.

Given starting approximation $x_0$, all methods   generate sequence of points $x_n$ which, under certain conditions, converges towards the solution
$\xi$. 

Method has __order of convergence__ equal to $r>0$ if there exists $A>0$ such that

$$
|\xi-x_{n+1}|\leq A|\xi-x_n|^r.$$

__Remark.__ Proofs of statements in this notebook and examples can be found in many textbooks.

"""

# ╔═╡ 7942e1ac-57f4-493c-9d09-de4154dcd7b9
md"""
# Bisection

Starting with the interval $[a,b]\equiv [a_0,b_0]$, construct a sequence of intervals  

$$
[a_0,b_0]\supset [a_1,b_1]\supset [a_2,b_2]\supset [a_3,b_3] \supset \cdots,$$

where $f(a_n)\cdot f(b_n)\leq 0$, and the sequence of points

$$
x_{n+1}=\frac{a_n+b_n}{2}.$$

__Speed of convergence__ is __linear__ since 

$$
|\xi-x_{n+1}|\leq \frac{1}{2}|\xi-x_n|,$$

and the __approximation error__ is bounded by

$$
|\xi-x_{n+1}|\leq \frac{1}{2}|a_n-b_n|.$$
"""

# ╔═╡ 9d707705-29a2-47fc-87b8-d39ff967efb5
function Bisection(f::Function,a::Number,b::Number,ϵ::Float64=1e-10)
    fa=f(a)
    fb=f(b)
    T=Float64
    x=T
    fx=T
    if fa*fb>zero(T)
        return "Wrong interval"
    end
    iter=0
    while b-a>ϵ && iter<1000
        x=(b+a)/2.0
        fx=f(x)
        if fa*fx<zero(T)
            b=x
            fb=fx
        else
            a=x
            fa=fx
        end
        iter+=1
        # @show x,fx
    end
    return x,fx,iter
end

# ╔═╡ 0233a80f-dfc6-420c-ae39-60828b19010a
md"""
## Examples

Given are functions and intervals:

$$
\begin{aligned}
f_1(x)&=e^x-x- \frac{5}{4},\quad &x\in [-2,2],\\
f_2(x)&=e^{-2x}\sin (6x)-\frac{2}{3}\,x-\frac{1}{2},\quad &x\in[-1,2],\\
f_3(x)&=x^3-6x+2,\quad &x\in[-4,4],\\
f_4(x)&=0.001\,x+0.5+\frac{\pi}{2}+\arctan(x),\quad &x\in[-1000,1000],\\
f_5(x)&=1000\,(x-4)-e^x,\quad &x\in[-10,10].
\end{aligned}$$
"""

# ╔═╡ 0be4ad52-497f-48a4-9f53-b6d1291beb58
begin
	f₁(x)=exp(x)-x-5.0/4
	(a₁,b₁)=(-1,1)
	f₂(x)=exp(-2x)*sin(6x)+2x/3-1.0/2
	(a₂,b₂)=(-1,2)
	f₃(x)=x^3-6*x+2
	(a₃,b₃)=(-4,4)
	f₄(x)=0.001x+0.5+π/2+atan(x)
	(a₄,b₄)=(-1000,1000)
	f₅(x)=1000(x-4)-exp(x)
	(a₅,b₅)=(-10,10)
end

# ╔═╡ 5a17ece9-a499-42c7-8dc8-cc9035f1b311
md"""
Using function graph we determine intervals which contain zeros, which we then compute and plot. 
"""

# ╔═╡ 910d7790-488d-48b8-a951-eee68774f79f
function Zeros(f,a,b,Intervals)
    plot(f,a,b,label="f(x)")
    for i=1:length(Intervals)
        iab=Intervals[i]
        x,y,iter=Bisection(f,iab[1],iab[2])
        scatter!([x],[y],label="Zero")
    end
    scatter!()
end

# ╔═╡ 5d233a49-774b-479c-8c80-5de46817614e
# Function f₁(x)
plot(f₁,a₁,b₁,label="f(x)")

# ╔═╡ aedcd7b9-8ef5-48a8-8ef8-d7dab4b733e0
Intervals₁=((-1,0),(0,1))

# ╔═╡ e7e92425-287a-482b-bb1e-8168b2b130f1
Zeros(f₁,a₁,b₁,Intervals₁)

# ╔═╡ 56a6bc1f-7ff8-464d-a864-ec7a19e3f9c3
# Function f₂(x)
plot(f₂,a₂,b₂,label="f(x)")

# ╔═╡ 7097ec9a-b556-4be7-a356-9c9d5003dd73
Intervals₂=((-1,-0.4),(-0.4,0.2),(0.2,0.6),(0.6,1))

# ╔═╡ bdb9d5d9-861c-4fd1-89d5-e0622120f944
Zeros(f₂,a₂,b₂,Intervals₂)

# ╔═╡ f7802590-744f-4175-ba60-5d801fa252ec
# Computation of an individual zero
Bisection(f₂,0.2,0.6)

# ╔═╡ d193fd54-c5a7-45f4-a61d-ae532bae4f44
# Function f₃(x)
plot(f₃,a₃,b₃,label="f(x)")

# ╔═╡ 52fbcdb5-ade7-48c0-bc09-0154f3d615b8
begin
	Intervals₃=((-4,-2),(0,1),(2,3))
	Zeros(f₃,a₃,b₃,Intervals₃)
end

# ╔═╡ e1f25e30-1ff7-4a91-b7cf-393c36c438b6
begin
	# Function f₄(x)
	Intervals₄=[(-600,-400)]
	Zeros(f₄,a₄,b₄,Intervals₄)
end

# ╔═╡ f33b37c9-2f88-422b-ac9f-40d7945ccb47
begin
	# Function f₅(x)
	Intervals₅=((0,5),(5,10))
	Zeros(f₅,a₅,b₅,Intervals₅)
end

# ╔═╡ 4c7dc70c-1832-461c-ad48-5e358b49661f
plot!(legend=:bottomright)

# ╔═╡ 69c1d80b-3d24-4ded-a226-391b337805e5
md"""
# Simple iterations

We are solving equation of the form 

$$
x=\varphi(x).  \qquad\qquad (2)$$

__Fixed Point Theorem. (Banach)__ Let

$$\varphi:[a,b]\to \mathbb{R}$$

be a __continuously differentiable function__ and let

$$\begin{aligned}
\varphi(x) &\in [a,b] \quad  \forall x\in [a,b], \\
|\varphi'(x)|&\leq q<1 \quad \forall x\in(a,b).
\end{aligned}\qquad\qquad (3)$$

Then there exists unique __fixed point__ $\xi \in [a,b]$ for which
$\xi=\varphi(\xi)$. 

Furthermore, for arbitrary starting point  $x_0\in[a,b]$, the sequence

$$
x_n=\varphi(x_{n-1}),\quad n=1,2,3,\ldots,$$

converges to $\xi$. The following __error bounds__ hold:

$$\begin{aligned}
|\xi-x_n|&\leq \displaystyle\frac{q^n}{1-q}|x_1-x_0|, \\
|\xi-x_n|&\leq \displaystyle\frac{q}{1-q}|x_n-x_{n-1}|, \\
|\xi-x_n|&\leq q|\xi-x_{n-1}|.
\end{aligned}$$

Thus, the convergence is __linear__.

__Problem.__ Find and study the proof of the theorem.
"""

# ╔═╡ d794399b-9c03-47a1-8a37-40f8545a0d46
function Iteration(φ::Function,x::Number,ϵ::Float64=1e-10)
    ξ=φ(x)
    iter=0
    while abs(x-ξ)>ϵ && iter<1000
        x=ξ
        ξ=φ(x)
        iter+=1
    end
    ξ,iter
end

# ╔═╡ 5cc494f2-a789-4455-9ac9-925428d274fe
md"""
In order to use simple iteration, we need to transform form (1) into form (2) such that the condition (3) holds.

To estimate derivative we can use the package `Calculus.jl` which approximates derivatives using finite differences, or the package
[`ForwardDiff.jl`](https://github.com/JuliaDiff/ForwardDiff.jl) which uses [automatic differentiation](https://en.wikipedia.org/wiki/Automatic_differentiation) and is more accurate. Symbolic computation with the package `SymPy.jl` can alse be used.
"""

# ╔═╡ 21495df3-1163-4d7f-afce-ad3600a5411a
varinfo(ForwardDiff.ForwardDiff)

# ╔═╡ 492a1939-ea1e-4912-b8d3-0330624f0ffd
md"""
## Example

Let us find zeros of the function $f_1(x)=e^x-x-\frac{5}{4}$. The form

$$
x=e^x-\frac{5}{4}\equiv \Phi(x)$$

can be used only to compute the negative zero, since in the vicinity of the positive zero $|\varphi'(x)|>1$.
For example, for  $x_0=1$, the sequence diverges very quickly, and for $x_0=0.6$, which is near positive zero, the sequence converges to the negative zero, without theoretical justification. 
"""

# ╔═╡ 7dab9934-414d-44dd-94b5-96a849184687
begin
	φ(x)=exp(x)-5.0/4
	plot([f₁,φ,x->ForwardDiff.derivative(φ,x)],-2.0,2.0,label=["f(x)" "φ(x)" "φ'(x)"])
end

# ╔═╡ 1e5be7af-8236-4a8c-ab10-4f988350bd6a
Iteration(φ,0.5)

# ╔═╡ ed701b40-1ca7-4eba-8ea4-7c2beef3e58d
Iteration(φ,1.0)

# ╔═╡ 8eb3ad86-6d79-4016-9f65-5797f51d779e
Iteration(φ,0.6)

# ╔═╡ 54214d6b-b630-43d7-afdc-96ebfd4800a1
md"""
The positive zero can be computed using the form

$$
x=\ln\big(x+\frac{5}{4}\big)\equiv \Psi(x).$$
"""

# ╔═╡ 70f9b3bc-ef11-4d7c-af8c-0abc262b92fa
begin
	Ψ(x)=log(x+5.0/4)
	plot([f₁,x->ForwardDiff.derivative(φ,x), x->ForwardDiff.derivative(Ψ,x)],
	    -1.0,1.0,label=["f(x)" "φ'(x)" "Ψ'(x)"])
end

# ╔═╡ 3ccbe741-478e-4b23-8a64-d6f3380d8b1b
scatter!([Iteration(φ,-0.5)[1],Iteration(Ψ,1.0)[1]],[0,0],label="Zeros")

# ╔═╡ 5dee83f3-2171-478a-b811-04c3ccef1a0e
md"""
## Square root of 2

Let us approximate $\sqrt{2}$, that is, find the positive solution of the equation 

$$
x^2-2=0.$$

The equation can be transformed into form (2) as

$$
x=\frac{2}{x},$$

but then $\varphi'(x)=-\displaystyle\frac{2}{x^2}$, so on the interval $[1,2]$ the condition (3) does not hold. Thus, we set

$$
\frac{x}{2}=\frac{1}{x},$$

or

$$
x=\frac{x}{2}+\frac{1}{x}=\frac{1}{2}\left(x+\frac{2}{x}\right)\equiv\varphi(x).$$

Very accurate approximation is attained after just 5 iterations!
"""

# ╔═╡ b27404ef-ac18-4dce-804c-07255ab001c2
begin
	φ₁(x)=(x+2.0/x)/2.0
	plot([φ₁,x->ForwardDiff.derivative(φ₁,x)],1.0,2.0,label=["φ₁(x)" "φ₁'(x)"])
end

# ╔═╡ 25640f91-4b47-4bcf-b360-499541077c80
Iteration(φ₁,1.0,1e-15), √2

# ╔═╡ 96619f40-3a10-11eb-1529-29d9aa539c03
# Manual computation using rational numbers
begin
	y=1//1
	y₁=(y+2//y)//2
	y₂=(y₁+2//y₁)//2
	y₃=(y₂+2//y₂)//2
	y₄=(y₃+2//y₃)//2
	y₅=(y₄+2//y₄)//2
end

# ╔═╡ bc84a732-3a10-11eb-11ed-1b1abc4f530b
y₅-√2

# ╔═╡ cf6ff8ea-2c1e-442a-bfc6-5eac4f7269dd
begin
	# Let us try √10
	φ₂(x)=(9x+10.0/x)/10.0
	plot([φ₂,x->ForwardDiff.derivative(φ₂,x)],3.0,4.0,label=["φ₂(x)" "φ₂'(x)"])
end

# ╔═╡ 4fa06267-97f2-409f-8025-c160362014a9
Iteration(φ₂,3.0,1e-10), √10 # 1e-15

# ╔═╡ 46bebd7f-69ba-402f-91aa-6ad864add60b
begin
	# Try sqrt(10) in a different way
	φ₃(x)=(4x+10.0/x)/5.0
	plot([φ₃,x->ForwardDiff.derivative(φ₃,x)],3.0,4.0)
end

# ╔═╡ 9b237357-f6dc-4054-bee8-1cac4a5aad07
Iteration(φ₃,3.0,1e-10), √10 # 1e-15

# ╔═╡ 4f3b3c3a-0cf4-4bcd-a8ad-9061c7653386
md"""
# Newton's method

__Newton's method__ or __Tangent method__ is based on the following idea: in the vicinity of the starting point $x_0$, the function $f(x)$ is approximated by the tangent line through the point $(x_0,f(x_0))$,

$$
f_1(x)=f(x_0)+f'(x_0)(x-x_0).$$

The next approximation is the intersection of the tangent line and the $x$-axis. In this way we obtain series of approximations:

$$
x_{n+1}=x_n-\frac{f(x_n)}{f'(x_n)},\quad n=0,1,2,\ldots \qquad\qquad (4)$$

We have the following:

__Theorem.__  Let the function $f:[a,b]\to \mathbb{R}$ satisfy:

*  $f''$ is continuous on  $(a,b)$,
*  $f(a)\cdot f(b)<0$,
*  $f'$ and $f''$ have consatnt sign on $(a,b)$, i 
*  for the starting approximation  $x_0\in [a,b]$ it holds $f(x_0)\cdot f''(x_0)>0$.

Then the sequence (4) converges to __unique__ solution $\xi$ of the equation $f(x)=0$. Also, the follwoing  __error bounds__ hold:

$$
\begin{aligned}
|\xi-x_n|&\leq \displaystyle\frac{M_2}{2m_1}(x_n-x_{n-1})^2, \\
|\xi-x_{n+1}|&\leq \displaystyle\frac{M_2}{2m_1}(\xi-x_{n})^2, \\
\end{aligned}$$

where

$$
M_2=\max_{x\in(a,b)}|f''(x)|,\quad
m_1=\min_{x\in(a,b)}|f'(x)|.$$

Therefore, the speed of convergence is __quadratic__.

## Example
"""

# ╔═╡ d8d1e5e8-13c6-43cd-b2ca-aa6e96c142e3
function Newton(f::Function,x::Number,ϵ::Float64=1e-10)
    ξ=x-f(x)/(x->ForwardDiff.derivative(f,x))(x)
    iter=0
    while abs(x-ξ)>ϵ && iter<100
        x=ξ
        ξ=x-f(x)/(x->ForwardDiff.derivative(f,x))(x)
        iter+=1
    end
    ξ,iter
end

# ╔═╡ 6258bcb3-2682-47e5-8db5-0aaf03f32807
begin
	f₆(x)=exp(-x)+x^2-2
	plot(f₆,-3,4,label="f(x)")
end

# ╔═╡ f9376e6e-0595-4045-be73-c59c07e7be34
md"""
Let us check the conditions of the theorem for the positive zero:
"""

# ╔═╡ dca6feb2-e2e3-4c1e-a568-794c610cc611
begin
	a₀=1
	b₀=2
	x₀=1.5
	plot([f₆,x->ForwardDiff.derivative(f₆,x),
	        x->ForwardDiff.derivative(x->ForwardDiff.derivative(f₆,x),x)],a₀,b₀, 
	    label=["f(x)" "f'(x)" "f''(x)"])
end

# ╔═╡ 0cf32332-f169-465a-9315-b9704a31badf
f₆(a₀)*f₆(b₀)<0, 
f₆(x₀)*(x->ForwardDiff.derivative(
        x->ForwardDiff.derivative(f₆,x),x))(x₀)>0

# ╔═╡ 6de6ad42-cad9-40b4-8503-78f1316b3b38
Newton(f₆,x₀) # 1e-15

# ╔═╡ 9ba461d0-3b8c-4ec7-9b14-8ea720a569c3
begin
	# Negative zero
	a=-1
	b=0
	x₁=-1.0
	f₆(a)*f₆(b)<0, 
	f₆(x₁)*(x->ForwardDiff.derivative(
	        x->ForwardDiff.derivative(f₆,x),x))(x₁)>0
end

# ╔═╡ 9609f70c-e5e9-4552-8ecb-8b95cf484d29
Newton(f₆,x₁)

# ╔═╡ 518e1f14-3254-4ca6-9d87-c695008ae79e
begin
	plot(f₆,-1.0,2)
	scatter!([Newton(f₆,x₀)[1],Newton(f₆,x₁)[1]],[0,0])
end

# ╔═╡ b7109008-05de-45db-a383-78b8a6973733
md"""
__Remark.__ If we choose starting values $x_0=1$ and $x_0=0$, respectively, the iterations will also converge towards desired zeros, albeit without theoretical justification: 
"""

# ╔═╡ 0fbd88ae-6c19-41ac-983b-83ba8bd97965
Newton(f₆,1)

# ╔═╡ 917b4fc3-e485-4801-8d16-0e32a22649b5
Newton(f₆,0)

# ╔═╡ 1db0ba3d-e82e-4e4d-8b5d-c1100d04fd2d
md"""
# Secant method

If, in the formula (4), the derivative $f'(x_n)$ is approximated by the finite difference (secant line) though the  __two__ previously obtained points,

$$
f'(x_n)\approx \frac{f(x_n)-f(x_{n-1})}{x_n-x_{n-1}}, $$

we obtain the sequence

$$
x_{n+1}=\frac{x_{n-1}f(x_n)-x_nf(x_{n-1})}{f(x_n)-f(x_{n-1})},\qquad f(x_n)\neq f(x_{n-1}), \quad n=1,2,3,\ldots.$$

Obviously, we need  __two__ initial approximations, $x_0,x_1\in[a,b]$. The convergence properties are similar to those of the Newton's method.
"""

# ╔═╡ dce97964-5846-4642-859e-f041b8a72d76
function Secant(f::Function,x::Number,ζ::Number,ϵ::Float64=1e-10)
    ξ=(x*f(ζ)-ζ*f(x))/(f(ζ)-f(x))
    iter=0
    while abs(ζ-ξ)>ϵ && iter<100
        x=ζ
        ζ=ξ
        ξ=(x*f(ζ)-ζ*f(x))/(f(ζ)-f(x))
        iter+=1
    end
    ξ,iter
end

# ╔═╡ 33402e76-139d-4062-b72d-57b45e1c3eec
Secant(f₆,-1,0), Secant(f₆,1,2)

# ╔═╡ Cell order:
# ╠═d0d29fe8-6542-4968-a695-03371ba85543
# ╠═3a20c826-03cd-4c06-a0c6-7dc657067feb
# ╠═7bab2ccd-e6d7-479c-abef-d9ed8eefaf2a
# ╟─1198a836-bb14-4b8e-9f30-160097bc4507
# ╟─7942e1ac-57f4-493c-9d09-de4154dcd7b9
# ╠═9d707705-29a2-47fc-87b8-d39ff967efb5
# ╟─0233a80f-dfc6-420c-ae39-60828b19010a
# ╠═0be4ad52-497f-48a4-9f53-b6d1291beb58
# ╟─5a17ece9-a499-42c7-8dc8-cc9035f1b311
# ╠═910d7790-488d-48b8-a951-eee68774f79f
# ╠═5d233a49-774b-479c-8c80-5de46817614e
# ╠═aedcd7b9-8ef5-48a8-8ef8-d7dab4b733e0
# ╠═e7e92425-287a-482b-bb1e-8168b2b130f1
# ╠═56a6bc1f-7ff8-464d-a864-ec7a19e3f9c3
# ╠═7097ec9a-b556-4be7-a356-9c9d5003dd73
# ╠═bdb9d5d9-861c-4fd1-89d5-e0622120f944
# ╠═f7802590-744f-4175-ba60-5d801fa252ec
# ╠═d193fd54-c5a7-45f4-a61d-ae532bae4f44
# ╠═52fbcdb5-ade7-48c0-bc09-0154f3d615b8
# ╠═e1f25e30-1ff7-4a91-b7cf-393c36c438b6
# ╠═f33b37c9-2f88-422b-ac9f-40d7945ccb47
# ╠═4c7dc70c-1832-461c-ad48-5e358b49661f
# ╟─69c1d80b-3d24-4ded-a226-391b337805e5
# ╠═d794399b-9c03-47a1-8a37-40f8545a0d46
# ╟─5cc494f2-a789-4455-9ac9-925428d274fe
# ╠═21495df3-1163-4d7f-afce-ad3600a5411a
# ╟─492a1939-ea1e-4912-b8d3-0330624f0ffd
# ╠═7dab9934-414d-44dd-94b5-96a849184687
# ╠═1e5be7af-8236-4a8c-ab10-4f988350bd6a
# ╠═ed701b40-1ca7-4eba-8ea4-7c2beef3e58d
# ╠═8eb3ad86-6d79-4016-9f65-5797f51d779e
# ╟─54214d6b-b630-43d7-afdc-96ebfd4800a1
# ╠═70f9b3bc-ef11-4d7c-af8c-0abc262b92fa
# ╠═3ccbe741-478e-4b23-8a64-d6f3380d8b1b
# ╟─5dee83f3-2171-478a-b811-04c3ccef1a0e
# ╠═b27404ef-ac18-4dce-804c-07255ab001c2
# ╠═25640f91-4b47-4bcf-b360-499541077c80
# ╠═96619f40-3a10-11eb-1529-29d9aa539c03
# ╠═bc84a732-3a10-11eb-11ed-1b1abc4f530b
# ╠═cf6ff8ea-2c1e-442a-bfc6-5eac4f7269dd
# ╠═4fa06267-97f2-409f-8025-c160362014a9
# ╠═46bebd7f-69ba-402f-91aa-6ad864add60b
# ╠═9b237357-f6dc-4054-bee8-1cac4a5aad07
# ╟─4f3b3c3a-0cf4-4bcd-a8ad-9061c7653386
# ╠═d8d1e5e8-13c6-43cd-b2ca-aa6e96c142e3
# ╠═6258bcb3-2682-47e5-8db5-0aaf03f32807
# ╟─f9376e6e-0595-4045-be73-c59c07e7be34
# ╠═dca6feb2-e2e3-4c1e-a568-794c610cc611
# ╠═0cf32332-f169-465a-9315-b9704a31badf
# ╠═6de6ad42-cad9-40b4-8503-78f1316b3b38
# ╠═9ba461d0-3b8c-4ec7-9b14-8ea720a569c3
# ╠═9609f70c-e5e9-4552-8ecb-8b95cf484d29
# ╠═518e1f14-3254-4ca6-9d87-c695008ae79e
# ╟─b7109008-05de-45db-a383-78b8a6973733
# ╠═0fbd88ae-6c19-41ac-983b-83ba8bd97965
# ╠═917b4fc3-e485-4801-8d16-0e32a22649b5
# ╟─1db0ba3d-e82e-4e4d-8b5d-c1100d04fd2d
# ╠═dce97964-5846-4642-859e-f041b8a72d76
# ╠═33402e76-139d-4062-b72d-57b45e1c3eec
