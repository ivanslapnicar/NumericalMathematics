### A Pluto.jl notebook ###
# v0.16.0

using Markdown
using InteractiveUtils

# ╔═╡ 0d4fbf25-a629-416a-bb67-b087144a5862
begin
	using PlutoUI, Plots, ODE
	plotly()
end

# ╔═╡ ac1342d7-a80c-4942-b424-6cb44904a07e
TableOfContents(title="📚 Table of Contents", aside=true)

# ╔═╡ b5d87000-5af5-11eb-0594-57573e58e809
md"""
# Differential Equations


Consider first order differential equation

$$
\frac{d}{dx} y(x)=f(x,y(x)),$$

with the given __initial condition__

$$
y(x_0)=y_0.$$

__Theorem.__ If the functions $f$ and $\displaystyle \frac{\partial f}{\partial y}$ are continuous in some neighbourhood of the point $(x_0,y_0)$, then, for some $\varepsilon > 0$, the given initial value problem has unique solution $y(x)$ on the interval $\displaystyle [x_{0}-\varepsilon ,x_{0}+\varepsilon ]$.

__Remark.__ The independent variable is often denoted by $t$ (time).
"""

# ╔═╡ 68a29859-4ba6-46ed-97ef-fdffec416003
md"""
# Euler's method

Taylor's formula in the neighbourhood of the point $x$ can be written as
(with the assumption that $y'''(x)$ is bounded)

$$
y(x+h)=y(x)+h\cdot y'(x) +\frac{1}{2}h^2 \cdot y''(x)+ O(h^3). \tag{1}$$

Starting from the point $x_0$ and the initial condition, for equally spaced points

$$
x_{k+1}=x_{k}+h,$$

__Euler's method__ approximates the value of the function $y$ using the first two terms of the Taylor series arround the point $x_k$:

$$
y_{k+1}=y_k+h \cdot f(x_k,y_k). \tag{2}$$

If $y(x_k)=y_k$, then the __local error__ of the method is the error in one step of the Taylor formula,

$$
y(x_{k+1})-y_{k+1}=\displaystyle\frac{1}{2}h^2 \cdot y''(x_k)+ O(h^3).$$

__Condition.__ During computations we generally have $y(x_k)\neq y_k$. Lagrange Mean Value Theorem implies

$$
f(x_k,y_k)=f(x_k,y(x_k))+f_y(x_k,\eta)(y_k-y(x_k)) \tag{3}$$

for some $\eta$ located between $y_k$ and $y(x_k)$. Using Taylor's formula (1) for $x=x_k$ and Euler's formula (2), we have

$$
y(x_{k+1})-y_{k+1}=y(x_k)+h\cdot f(x_k,y(x_k)) +\frac{1}{2}h^2 \cdot y''(x_k)+ O(h^3)
-y_k-h \cdot f(x_k,y_k).$$

Inserting Lagrange formula (3) yields

$$
y(x_{k+1})-y_{k+1}=y(x_k)-y_k + h\cdot [f(x_k,y_k)-f_y(x_k,\eta)(y_k-y(x_k))] - h \cdot f(x_k,y_k)+O(h^2),$$

or

$$
y(x_{k+1})-y_{k+1}=(y(x_k)-y_k) [1+h\cdot f_y(x_k,\eta)]+O(h^2). \tag{4}$$

We conclude that during the propagation error decreases if $f_y<0$, and increases if $f_y>0$, in which case the problem is badly conditioned.

__Global error__ is the error in the point $x$, after $n$ steps which were needed for the method to reach $x$ from the initial point. The global error is comprised of error in $y_k$ which is due to accumulation of previous local errors, and error in the Taylor formula in the $n$-th step. If the increment is $h$, then, obviously, $n=\displaystyle\frac{x-x_0}{h}$ so we expect that the global error is (approximately) proportional to $h$. This can also be concluded from the bound (4).
We say that Euler's method is a __first-order method__.
"""

# ╔═╡ fce92627-56f1-482e-983a-083b7281fbcb
function Euler(f::Function,y₀::T,x::T1) where {T,T1}
    h=x[2]-x[1]
    y=Array{T}(undef,length(x))
    y[1]=y₀
    for k=2:length(x)
        y[k]=y[k-1]+h*f(x[k-1],y[k-1])
    end
    y
end

# ╔═╡ 934a054a-f31d-497e-bcf8-295dac6db83b
md"""
## Example 1

The solution of the initial value problem

$$
y'=x+y,\quad y(0)=1,$$

is

$$
y=2e^x-x-1.$$

Compute the solution $y$ for $x\in[0,1]$.
"""

# ╔═╡ 4f30e023-2c8d-42b7-a124-c240e5464769
begin
	# 10 sub-intervals on the interval [0,1]
	x₁=range(0,stop=1,length=11)
	f₁(x,y)=x+y
	y₁=Euler(f₁,1.0,x₁)
end

# ╔═╡ 6f9f8bbd-1318-4d66-bdd6-6d02f778f8b8
begin
	# Plot exact solution and the computed points
	solution₁(x)=2*exp(x)-x-1
	plot(solution₁,0,1,xlabel="x",ylabel="y",
	    label="Exact solution",legend=:topleft)
	plot!(x₁,y₁,label="Euler()")
	scatter!(x₁,y₁,label="Computed points")
end

# ╔═╡ 5e6c51c1-c7a9-4542-89fd-9c7b0846437c
md"""
## Error bound

Notice that the given problem is badly conditioned.
Accurate bound for the global error is given by the following theorem:

__Theorem.__ Let the function $f(x,y)$ have continuous first partial derivatives od some rectangle $D\subseteq \mathbb{R}^2$ and let

$$
K=\max_{(x,y)\in D}|f_y(x, y)|<\infty,\quad   M= \max_{(x,y)\in D}|(f_x+f\cdot f_y)(x,y)|<\infty.$$

If

$$(x_k,y_k)\in D\quad (x_k,y(x_k))\in D, \quad k=0,1,2,\ldots n,$$

and if $Kh<1$, then the error is bounded by

$$
|y(x)-y_n|\leq \frac{M}{2K}\left(e^{Kx}−1\right)h. \tag{5}$$

_Proof._ See [Glenn Ledder, Error in Euler’s Method](http://www.math.unl.edu/~gledder1/Math447/EulerError).
"""

# ╔═╡ 54c26825-46b9-44f6-ae65-31f4be69bf3b
md"""
__Problem.__ Determine the number of steps $n$ such that the value $y(1)$ in Example 1 is computed with an error less than $\epsilon=0.01$.

It holds $f_y(x,y)=1$ so $K=1$. For positive $x$ and $y$ we have $y'>0$ so $y$ is an increasing function. We can estimate the value $y(1)$ with $4$ since in the previous computation with $10$ sub-intervals $y_{10}\approx 3.19$. Thus, $M=5$. Inserting $K$, $M$ and $h=\displaystyle\frac{1-0}{n}$ into the formula (5), we have

$$
\frac{5}{2}(2.7183-1)\frac{1}{n} \approx 4.3 \frac{1}{n} < 0.01,$$

so we can take $n=500$.
"""

# ╔═╡ 22f0b250-0b37-4593-bdea-4bd13c19f2bf
begin
	# 500 sub-intervals on the interval [0,1]
	xx₁=range(0,stop=1,length=501)
	yy₁=Euler(f₁,1.0,xx₁)
	plot(solution₁,0,1,xlabel="x",ylabel="y",
	    label="Exact solution",legend=:topleft)
	plot!(xx₁,yy₁,label="Euler()")
end

# ╔═╡ bfc718ab-e392-4fb8-82c6-0c889b1977e9
# Check error
solution₁(1)-yy₁[end]

# ╔═╡ 12ceb8e5-f976-4639-9175-047e84321c85
md"""
## Example 2

The solution of the problem

$$
y'=30(\sin x-y), \quad y(0)=0,$$

is

$$
y(x)=\frac{30}{901}(30\sin x-\cos x+e^{-30x}).$$

Compute the solution $y$ for $x\in[0,1]$.
"""

# ╔═╡ 1e36310b-6f42-44a1-8b04-e82391cb33aa
begin
	# 100 sub-intervals on the interval [0,1]
	f₂(x,y)=30(sin(x)-y)
	x₂=range(0,stop=1,length=101)
	y₂=Euler(f₂,0.0,x₂)
	solution₂(x)=30(30*sin(x)-cos(x)+exp(-30x))/901
	plot(solution₂,0,1,xlabel="x",ylabel="y",
	    label="Exact solution",legend=:topleft)
	plot!(x₂,y₂,label="Euler()")
end

# ╔═╡ 8899d0fd-e0d9-424f-98e4-5547d029a36a
md"""
__Problem.__ Estimate the accuracy of the computed value $y(1)$ using the bound (5).
"""

# ╔═╡ 2d44d600-5af7-11eb-0958-b3eab7ae938d
md"""
# Runge-Kutta methods

Euler's method is a first-order method, and it not accurate enough. Thus, in practice we use methods of higher order which approximate $y(x)$ in the point $x_{k+1}$ using values of the function $f(x,y)$ in several points from the interval

$$[x_k,x_{k+1}]\equiv[x_k,x_k+h].$$

## Heun's method

$$\begin{aligned}
k_1&=hf(x_k,y_k),\\
k_2&=hf(x_k+h,y_k+k_1),\\
y_{k+1}&=y_k+\frac{1}{2}(k_1+k_2).
\end{aligned}$$

Heun's method is derived as ffollows: integrating given differential equation gives

$$
\int_{x_k}^{x_{k+1}} \frac{dy}{dx}\, dx = y(x_{k+1})-y(x_k)= \int_{x_k}^{x_{k+1}}
f(x,y(x))\, dx.$$

Trapezoid rule gives

$$
y(x_{k+1})= y(x_k) +\frac{h}{2}[f(x_k,y(x_k))+f(x_{k+1},y(x_{k+1}))].$$

Setting $y(x_k)=y_k$, $y(x_{k+1})=y_{k+1}$, and applying Euler's formula (2) to the right hand side gives Heun's formulas.

The Trapezoid formula has an error of the order of magnitude $O(h^2)$, so the global error of Heun's method is also of the order of magnitude $O(h^2)$. Notice that the function $f(x,y)$ is evaluated twice in each step.
"""

# ╔═╡ 8b933ce6-7fea-4190-940a-bd7720fc93be
md"""

## Standard Runge-Kutta method

$$\begin{aligned}
k_1&=hf(x_k,y_k),\\
k_2&=hf\big(x_k+\frac{h}{2},y_k+\frac{k_1}{2}\big),\\
k_3&=hf\big(x_k+\frac{h}{2},y_k+\frac{k_2}{2}\big),\\
k_4&=hf(x_k+h,y_k+k_3),\\
y_{k+1}&=y_k+\frac{1}{6}(k_1+2k_2+2k_3+k_4).
\end{aligned}$$

Standard Runge-Kutta method is the method of order $4$ - order of magnitude of the global error is $O(h^4)$, and the function $f(x,y)$ is evaluated four times in each step.

Local error of the standard Runge-Kutta method is $O(h^5)$.
"""

# ╔═╡ cb8721aa-7f4b-4672-a802-da84d9c67b4a
function RungeKutta4(f::Function,y₀::T,x::T1) where {T,T1}
    h=x[2]-x[1]
    y=Array{T}(undef,length(x))
    y[1]=y₀
    for k=2:length(x)
        ξ=x[k-1]
        η=y[k-1]
        k₁=h*f(ξ,η)
        k₂=h*f(ξ+h/2,η+k₁/2)
        k₃=h*f(ξ+h/2,η+k₂/2)
        k₄=h*f(ξ+h,η+k₃)
        y[k]=η+(k₁+2*k₂+2*k₃+k₄)/6.0
    end
    y
end

# ╔═╡ ad54ce7f-8caa-46f8-9a57-f3710ca89f40
md"""
## Example

Let us solve problems from Examples 1 and 2. For Example 1 the numerical solution graphically overlaps the exact solution. For Example 2 the solution using the function `RungeKutta4()` is an order of magnitude more accurate than the solution obtained using the function  red `Euler()`.
"""

# ╔═╡ 4adf6672-5ae1-42e1-85e9-e522d6069fab
begin
	y₃=RungeKutta4(f₁,1.0,x₁)
	plot(solution₁,0,1,xlabel="x",ylabel="y",
	    label="Exact solution",legend=:topleft)
	plot!(x₁,y₃,label="RunkeKutta4()")
end

# ╔═╡ a262689c-74f3-46bb-bd55-1aa7eec5379f
begin
	x₄=range(0,stop=1,length=21)
	yEuler=Euler(f₂,0.0,x₄)
	yRK4=RungeKutta4(f₂,0.0,x₄)
	plot(solution₂,0,1,xlabel="x",ylabel="y",
	    label="Exact solution",legend=:topleft)
	plot!(x₄,[yEuler,yRK4],label=["Euler()" "RungeKutta4()"])
end

# ╔═╡ 481e7da8-95f2-4972-863c-fcd7a79cf417
solution₂(1), yEuler[end],yRK4[end]

# ╔═╡ d52999d7-2644-4793-a36a-413ad8011e09
md"""
# Existing routines

Most programming languages have built-in routines for numerical solution of ordinary differential equations. For example,

* Matlab has commands `ode*` (see [Matlab, Ordinary Differential Equations](https://www.mathworks.com/help/matlab/ordinary-differential-equations.html?searchHighlight=ordinary%20differential&s_tid=srchtitle)), a
* Julia has the package [ODE.jl](https://github.com/JuliaODE/ODE.jl).

Standard Runge-Kutta method is implemented in the function `ode4()`, and Heun's method is implemented in the function `ODE.ode2_heun()`.

__Remark.__ Function `ODE.ode2_heun()` is not visible with the function `varinfo()` since it is not being exported, but it can be found in the file `runge_kutta.jl`.
"""

# ╔═╡ 1a0815c9-900d-4ba5-855a-7fa582e8d498
# Check!
# varinfo(ODE)

# ╔═╡ d8e15baa-34a9-447d-90d2-8dfe8065ee1b
methods(ode4)

# ╔═╡ f35f1bb2-bc2c-4bb7-9dc9-367bc5a6c15d
methods(ODE.ode2_heun)

# ╔═╡ e50cf9c2-972d-4916-a3ca-6c20a0c73519
# Let us solve the problem from Example 2.
# Computed values yₖ are the first elements of the output vectors.
yode4=ode4(f₂,0.0,range(0,stop=1,length=21))[2]

# ╔═╡ f4d37a14-e6ef-4400-bb30-a3a4405bb675
yode2=ODE.ode2_heun(f₂,0.0,range(0,stop=1,length=21))[2];

# ╔═╡ a7841a3d-c7e4-4c09-9012-a3890c9e1195
# Compare solutions
[yRK4 yode4 yRK4-yode4 yode2 yode4-yode2]

# ╔═╡ 471ddc14-d7ba-4d94-a902-9ca459be5ce3
md"""
# Systems of differential equations

Let us solve the system of $n$ equations

$$\begin{aligned}
y_1'(x)&=f_1(x,y_1,y_2,\ldots,y_n),\\
y_2'(x)&=f_2(x,y_1,y_2,\ldots,y_n),\\
&\vdots \\
y_n'(x)&=f_2(x,y_1,y_2,\ldots,y_n)
\end{aligned}$$

and $n$ unknown functions $y_1,y_2,\ldots,y_n$ satisfying initial conditions

$$
y_i(x_0)=\zeta_i.$$

Using notation

$$
f=\begin{bmatrix}f_1\\ f_2\\ \vdots\\f_n\end{bmatrix},\quad
y=\begin{bmatrix}y_1\\ y_2\\ \vdots \\y_n\end{bmatrix},\quad
\zeta=\begin{bmatrix}\zeta_1\\ \zeta_2\\ \vdots\\ \zeta_n\end{bmatrix},$$

we can write the given system in the vector form as

$$
y'(x)=f(x,y),\quad y(x_0)=\zeta.$$

The system is successfully solved using Euler's method and Runge-Kutta methods in vector form.
"""

# ╔═╡ 705ebcb9-bfa6-4357-8610-27c43c0c3f96
md"""
## Lotka-Volterra equations

Modeling of the __predator-prey__ system gives __Lotka-Volterra__ equations (see  [Lotka-Volterra equations](https://en.wikipedia.org/wiki/Lotka%E2%80%93Volterra_equations)):

$$\begin{aligned}
\frac{dZ}{dt}&=z\,Z-a\, Z\, V = Z\,(z-a\, V), \qquad(6)\\
\frac{dV}{dt}&=-v\,V+b\, Z\, V = V\,(-v+b\, Z), \quad v,z,a,b>0,
\end{aligned}$$

with initial conditions

$$
V(t_0)=V_0,\qquad Z(t_0)=Z_0.$$

__Stable states__ are the states in which there is no change, that is, the states where both derivatives are zero. These are the __trivial__ stable state, $V=Z=0$, and

$$
V=\frac{z}{a},\qquad Z=\frac{b}{v}.$$

In __phase space__, elimation of the independent variable $t$ yields one differential equation:

$$
\frac{dV}{dZ}=\frac{\displaystyle\frac{dV}{dt}}{\displaystyle\frac{dZ}{dt}}=
\frac{V\,(-v+b\, Z)}{Z\,(z-a\, V)}.$$

This is a separable differential equations which has implicit solution:

$$
V^z \, Z^v = C\,  e^{aV}\, e^{bZ},\qquad
C=\frac{V_0^{z}\, Z_0^{v}}{\displaystyle e^{a V_0}e^{b Z_0}}.  \tag{7} $$

Let us solve the system for populations of predators (say, wolves) $V$ and prey (say, hares) $Z$ with

$$
v=0.02, \quad z=0.06,\quad a=0.001,\quad b=0.00002,\quad V(0)=30, \quad Z(0)=800.$$

__Remark.__ Functions `Euler()`, `RungeKutta4()` and the functions from the package `ODE.jl` are already capable of solving systems in the vector form.
"""

# ╔═╡ a6cd2db9-f4b6-46cb-aaa0-59d53430056c
begin
	# Given problem: y=[V,Z], t₀=0, y₀=y(0)=[30,800]
	# 356 days with step of 1/10 of day
	t=range(0,stop=365,length=3651)
	# System parameters
	v=0.02
	z=0.06
	a=0.001
	b=0.00002
	# Initial population of wolves
	V₀=30.0
	# Initial population of hares
	Z₀=800.0
	# Starting point
	y₀=[V₀,Z₀]
	# Vector function
	fVZ(t,y)=[y[1]*(-v+b*y[2]),y[2]*(z-a*y[1])]
	# Solutions
	yₗ=Euler(fVZ,y₀,t)
end

# ╔═╡ 78eae1ea-5ce9-4bda-81d3-4b056ee7df1e
begin
	# Scale Z to Z/10 to make the graph more readable
	V=[yₗ[i][1] for i=1:length(yₗ)]
	Z=[yₗ[i][2] for i=1:length(yₗ)]
	plot(t,[V,Z/10],xlabel="t ( days )", ylabel="Population",
	    label=["Wolves" "Hares / 10"])
end

# ╔═╡ 00796b0e-f2d5-4646-af33-06d81d13fc0c
begin
	# Compare solutions obtained using RungeKutta4() and ode4()
	yRK4ₗ=RungeKutta4(fVZ,y₀,t)
	yode4ₗ=ode4(fVZ,y₀,t)[2]
	[yₗ yRK4ₗ yode4ₗ]
end

# ╔═╡ 0b5ad338-3628-477e-bd0b-8cdda13d345d
# Plot solution in the phase space
plot(Z,V,xlabel="Z",ylabel="V", label=:none)

# ╔═╡ 5308aa28-6da6-4943-b393-c08592b98c27
md"""
## Scaled Lotka-Volterra equations

Plotting the exact solution (7) in the phase space is impractical since plotting of implicitly defined functions on a large area is time consuming. However, using transformations (see [Modeling Complex Systems, section 2.1](http://www.springer.com/us/book/9781441965615))

$$
X=\frac{b}{v}Z,\quad Y=\frac{a}{z}V,\quad \tau=\sqrt{z\cdot v}\, t,\quad \rho=\sqrt{\displaystyle\frac{z}{v}},$$

the system (6) can be converted to __dimensionless variables__
in the __scaled time__ $\tau$:

$$\begin{aligned}
\frac{dX}{d\tau}&=\rho\, X\,(1-Y),\\
\frac{dY}{d\tau}&=-\frac{1}{\rho}\, Y\,(1-X).
\end{aligned}$$

This system depends on only __one__ parameter $\rho$. The system has non-trivial stable state $X=Y=1$, and the solution (7) in the phase space is

$$
Y \, X^{1/\rho^2} = C\,  e^{Y}\, e^{X/\rho^2},\qquad
C=\frac{Y_0\, X_0^{1/\rho^2}}{\displaystyle e^{Y_0}e^{X_0/\rho^2}}.$$


Let us solve the system in the dimensionless form and compare solutions graphically:

"""

# ╔═╡ 00f084ce-a5ff-4e18-913f-9d14ae346bd9
begin
	ρ=sqrt(z/v)
	τ=range(0,stop=365*sqrt(z*v),length=3651)
	Y₀=[Z₀*b/v,V₀*a/z]
	fXY(τ,y)=[ρ*y[1]*(1-y[2]),-y[2]*(1-y[1])/ρ]
	Y=Euler(fXY,Y₀,τ);
end

# ╔═╡ 0ab04df9-d9a0-406a-b45f-299a0212fd9b
begin
	X₁=[Y[i][1] for i=1:length(Y)]
	Y₁=[Y[i][2] for i=1:length(Y)]
	# RThe solutions overlap
	plot(t,[V,Z/10],xlabel="t ( days )", ylabel="Population",
	    label=["Wolves" "Hares / 10"])
	plot!(τ/sqrt(z*v),[Y₁*z/a,X₁*v/(10b)],
	    label=["Wolves (dimensionless)" "Hares / 10 (dimensionless)"])
end

# ╔═╡ e580ec4b-73f9-4522-9cb3-4d24ebe1e027
# Plot dimensionless solution in the phase space
plot(X₁,Y₁,xlabel="Z",ylabel="V", label=:none)

# ╔═╡ 3ada3bbf-cc3b-46b4-adbd-c3685568781a
md"""
# Differential equations of higher order

Using substitutions, differential equation of higher order can be reduced to a system of differential equations of the first order.

## Example

The solution of the initial value problem

$$
y'''+y''=x,\qquad y(0)=0,\quad y'(0)=0,\quad y''(0)=0$$

is

$$
y(x)=-1+x+e^{-x}+\frac{x^3}{6}-\frac{x^2}{2}.$$

Substitutions

$$
y'=u,\quad y''=v,$$

yield the system

$$\begin{aligned}
y'&=u \\
u'&=v \\
v'&=-v+x\\
\end{aligned}$$

with initial conditions

$$
y(0)=0,\quad u(0)=0,\quad v(0)=0.$$

"""

# ╔═╡ 39d0b178-70bd-4815-bee8-3738bdaea250
begin
	x₅=range(0,stop=2,length=201)
	# Vector function
	f₅(x,y)=[y[2],y[3],-y[3]+x]
	# Initial conditions
	y₅=[0.0,0,0]
	# Compute solution using three methods
	yEuler₅=Euler(f₅,y₅,x₅)
	yRK4₅=RungeKutta4(f₅,y₅,x₅)
	yode4₅=ode4(f₅,y₅,x₅)
	# Computed values yₖ are the first elements of the output vectors
	YEuler₅=[yEuler₅[i][1] for i=1:length(yEuler₅)]
	YRK4₅=[yRK4₅[i][1] for i=1:length(yEuler₅)]
	# Computed values yₖ are the first elements of the vectors
	# of the second output array
	Yode4₅=[yode4₅[2][i][1] for i=1:length(yEuler₅)]
	# Exact solution
	solution₅(x)=-1+x+exp(-x)+x^3/6-x^2/2
	# Plotting
	plot(solution₅,0,2,xlabel="x",ylabel="y",label="Exact solution",legend=:topleft)
	plot!(x₅,[YEuler₅,YRK4₅,Yode4₅],xlabel="x",ylabel="y",label=["Euler()" "RungeKutta4()" "ode4()"])
end

# ╔═╡ 12c6cb02-6bc0-4558-bc10-9ca1d8cd029a
# Norm of errors in the used points
# import LinearAlgebra; LinearAlgebra.norm(solution₅.(x)-Y)
√(sum((solution₅.(x₅)-YRK4₅).^2))

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
ODE = "c030b06c-0b6d-57c2-b091-7029874bd033"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
ODE = "~2.13.0"
Plots = "~1.22.1"
PlutoUI = "~0.7.9"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "84918055d15b3114ede17ac6a7182f68870c16f7"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.1"

[[ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[ArrayInterface]]
deps = ["Compat", "IfElse", "LinearAlgebra", "Requires", "SparseArrays", "Static"]
git-tree-sha1 = "b8d49c34c3da35f220e7295659cd0bab8e739fed"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "3.1.33"

[[Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[BitTwiddlingConvenienceFunctions]]
deps = ["Static"]
git-tree-sha1 = "652aab0fc0d6d4db4cc726425cadf700e9f473f1"
uuid = "62783981-4cbd-42fc-bca8-16325de8dc4b"
version = "0.1.0"

[[Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[CPUSummary]]
deps = ["Hwloc", "IfElse", "Static"]
git-tree-sha1 = "ed720e2622820bf584d4ad90e6fcb93d95170b44"
uuid = "2a0fbf3d-bb9c-48f3-b0a9-814d99fd7ab9"
version = "0.1.3"

[[Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "f2202b55d816427cd385a9a4f3ffb226bee80f99"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+0"

[[ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "4ce9393e871aca86cc457d9f66976c3da6902ea7"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.4.0"

[[CloseOpenIntervals]]
deps = ["ArrayInterface", "Static"]
git-tree-sha1 = "ce9c0d07ed6e1a4fecd2df6ace144cbd29ba6f37"
uuid = "fb6a15b2-703c-40df-9091-08a04967cfa9"
version = "0.1.2"

[[ColorSchemes]]
deps = ["ColorTypes", "Colors", "FixedPointNumbers", "Random"]
git-tree-sha1 = "9995eb3977fbf67b86d0a0a0508e83017ded03f2"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.14.0"

[[ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "024fe24d83e4a5bf5fc80501a314ce0d1aa35597"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.0"

[[Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[CommonSolve]]
git-tree-sha1 = "68a0743f578349ada8bc911a5cbd5a2ef6ed6d1f"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.0"

[[CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "4866e381721b30fac8dda4c8cb1d9db45c8d2994"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.37.0"

[[CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f74e9d5388b8620b4cee35d4c5a618dd4dc547f4"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.3.0"

[[Contour]]
deps = ["StaticArrays"]
git-tree-sha1 = "9f02045d934dc030edad45944ea80dbd1f0ebea7"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.5.7"

[[DEDataArrays]]
deps = ["ArrayInterface", "DocStringExtensions", "LinearAlgebra", "RecursiveArrayTools", "SciMLBase", "StaticArrays"]
git-tree-sha1 = "31186e61936fbbccb41d809ad4338c9f7addf7ae"
uuid = "754358af-613d-5f8d-9788-280bf1605d4c"
version = "0.2.0"

[[DataAPI]]
git-tree-sha1 = "bec2532f8adb82005476c141ec23e921fc20971b"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.8.0"

[[DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "7d9d316f04214f7efdbb6398d545446e246eff02"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.10"

[[DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[DiffEqBase]]
deps = ["ArrayInterface", "ChainRulesCore", "DEDataArrays", "DataStructures", "Distributions", "DocStringExtensions", "FastBroadcast", "ForwardDiff", "FunctionWrappers", "IterativeSolvers", "LabelledArrays", "LinearAlgebra", "Logging", "MuladdMacro", "NonlinearSolve", "Parameters", "PreallocationTools", "Printf", "RecursiveArrayTools", "RecursiveFactorization", "Reexport", "Requires", "SciMLBase", "Setfield", "SparseArrays", "StaticArrays", "Statistics", "SuiteSparse", "ZygoteRules"]
git-tree-sha1 = "420ad175d5e420e2c55a0ed8a9c18556e6735f80"
uuid = "2b5f629d-d688-5b77-993f-72d75c75574e"
version = "6.73.2"

[[DiffResults]]
deps = ["StaticArrays"]
git-tree-sha1 = "c18e98cba888c6c25d1c3b048e4b3380ca956805"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.0.3"

[[DiffRules]]
deps = ["NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "7220bc21c33e990c14f4a9a319b1d242ebc5b269"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.3.1"

[[Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[Distributions]]
deps = ["ChainRulesCore", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SparseArrays", "SpecialFunctions", "Statistics", "StatsBase", "StatsFuns"]
git-tree-sha1 = "f4efaa4b5157e0cdb8283ae0b5428bc9208436ed"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.16"

[[DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "a32185f5428d3986f47c2ab78b1f216d5e6cc96f"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.5"

[[Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3f3a2501fa7236e9b911e0f7a588c657e822bb6d"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.3+0"

[[Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b3bfd02e98aedfa5cf885665493c5598c350cd2f"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.2.10+0"

[[ExprTools]]
git-tree-sha1 = "b7e3d17636b348f005f11040025ae8c6f645fe92"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.6"

[[FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "Pkg", "Zlib_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "d8a578692e3077ac998b50c0217dfd67f21d1e5f"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.0+0"

[[FastBroadcast]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "26be48918640ce002f5833e8fc537b2ba7ed0234"
uuid = "7034ab61-46d4-4ed7-9d0f-46aef9175898"
version = "0.1.8"

[[FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "Statistics"]
git-tree-sha1 = "caf289224e622f518c9dbfe832cdafa17d7c80a6"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "0.12.4"

[[FiniteDiff]]
deps = ["ArrayInterface", "LinearAlgebra", "Requires", "SparseArrays", "StaticArrays"]
git-tree-sha1 = "8b3c09b56acaf3c0e581c66638b85c8650ee9dca"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.8.1"

[[FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "NaNMath", "Printf", "Random", "SpecialFunctions", "StaticArrays"]
git-tree-sha1 = "b5e930ac60b613ef3406da6d4f42c35d8dc51419"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.19"

[[FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "87eb71354d8ec1a96d4a7636bd57a7347dde3ef9"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.4+0"

[[FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[FunctionWrappers]]
git-tree-sha1 = "241552bc2209f0fa068b6415b1942cc0aa486bcc"
uuid = "069b7b12-0de2-55c6-9aab-29f3d0a68a2e"
version = "1.1.2"

[[Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "dba1e8614e98949abfa60480b13653813d8f0157"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.5+0"

[[GR]]
deps = ["Base64", "DelimitedFiles", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Printf", "Random", "Serialization", "Sockets", "Test", "UUIDs"]
git-tree-sha1 = "c2178cfbc0a5a552e16d097fae508f2024de61a3"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.59.0"

[[GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Pkg", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "ef49a187604f865f4708c90e3f431890724e9012"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.59.0+0"

[[GeometryBasics]]
deps = ["EarCut_jll", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "58bcdf5ebc057b085e58d95c138725628dd7453c"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.1"

[[Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "7bf67e9a481712b3dbe9cb3dac852dc4b1162e02"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.68.3+0"

[[Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[HTTP]]
deps = ["Base64", "Dates", "IniFile", "Logging", "MbedTLS", "NetworkOptions", "Sockets", "URIs"]
git-tree-sha1 = "60ed5f1643927479f845b0135bb369b031b541fa"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "0.9.14"

[[HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "8a954fed8ac097d5be04921d595f741115c1b2ad"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+0"

[[HostCPUFeatures]]
deps = ["BitTwiddlingConvenienceFunctions", "IfElse", "Libdl", "Static"]
git-tree-sha1 = "3169c8b31863f9a409be1d17693751314241e3eb"
uuid = "3e5b6fbb-0976-4d2c-9146-d79de83f2fb0"
version = "0.1.4"

[[Hwloc]]
deps = ["Hwloc_jll"]
git-tree-sha1 = "92d99146066c5c6888d5a3abc871e6a214388b91"
uuid = "0e44f5e4-bd66-52a0-8798-143a42290a1d"
version = "2.0.0"

[[Hwloc_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3395d4d4aeb3c9d31f5929d32760d8baeee88aaf"
uuid = "e33a78d0-f292-5ffc-b300-72abe9b543c8"
version = "2.5.0+0"

[[IfElse]]
git-tree-sha1 = "28e837ff3e7a6c3cdb252ce49fb412c8eb3caeef"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.0"

[[IniFile]]
deps = ["Test"]
git-tree-sha1 = "098e4d2c533924c921f9f9847274f2ad89e018b8"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.0"

[[InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[Intervals]]
deps = ["Dates", "Printf", "RecipesBase", "Serialization", "TimeZones"]
git-tree-sha1 = "323a38ed1952d30586d0fe03412cde9399d3618b"
uuid = "d8418881-c3e1-53bb-8760-2df7ec849ed5"
version = "1.5.0"

[[IrrationalConstants]]
git-tree-sha1 = "f76424439413893a832026ca355fe273e93bce94"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.0"

[[IterTools]]
git-tree-sha1 = "05110a2ab1fc5f932622ffea2a003221f4782c18"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.3.0"

[[IterativeSolvers]]
deps = ["LinearAlgebra", "Printf", "Random", "RecipesBase", "SparseArrays"]
git-tree-sha1 = "1a8c6237e78b714e901e406c096fc8a65528af7d"
uuid = "42fd0dbc-a981-5370-80f2-aaf504508153"
version = "0.9.1"

[[IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "642a199af8b68253517b80bd3bfd17eb4e84df6e"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.3.0"

[[JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "8076680b162ada2a031f707ac7b4953e30667a37"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.2"

[[JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d735490ac75c5cb9f1b00d8b5509c11984dc6943"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.0+0"

[[LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[LaTeXStrings]]
git-tree-sha1 = "c7f1c695e06c01b95a67f0cd1d34994f3e7db104"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.2.1"

[[LabelledArrays]]
deps = ["ArrayInterface", "LinearAlgebra", "MacroTools", "StaticArrays"]
git-tree-sha1 = "bdde43e002847c34c206735b1cf860bc3abd35e7"
uuid = "2ee39098-c373-598a-b85f-a56591580800"
version = "1.6.4"

[[Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "Printf", "Requires"]
git-tree-sha1 = "a4b12a1bd2ebade87891ab7e36fdbce582301a92"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.6"

[[LayoutPointers]]
deps = ["ArrayInterface", "LinearAlgebra", "ManualMemory", "SIMDTypes", "Static"]
git-tree-sha1 = "d2bda6aa0b03ce6f141a2dc73d0bcb7070131adc"
uuid = "10f19ff3-798f-405d-979b-55457f8fc047"
version = "0.1.3"

[[LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "761a393aeccd6aa92ec3515e428c26bf99575b3b"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+0"

[[Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "7739f837d6447403596a75d19ed01fd08d6f56bf"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.3.0+3"

[[Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "42b62845d70a619f063a7da093d995ec8e15e778"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+1"

[[Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "340e257aada13f95f98ee352d316c3bed37c8ab9"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.3.0+0"

[[Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[LinearAlgebra]]
deps = ["Libdl"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[LogExpFunctions]]
deps = ["ChainRulesCore", "DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "34dc30f868e368f8a17b728a1238f3fcda43931a"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.3"

[[Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[LoopVectorization]]
deps = ["ArrayInterface", "CPUSummary", "CloseOpenIntervals", "DocStringExtensions", "HostCPUFeatures", "IfElse", "LayoutPointers", "LinearAlgebra", "OffsetArrays", "PolyesterWeave", "Requires", "SLEEFPirates", "Static", "ThreadingUtilities", "UnPack", "VectorizationBase"]
git-tree-sha1 = "d4046cd65d9c9cb328741d66bccff8eec019f81b"
uuid = "bdcacae8-1622-11e9-2a5c-532679323890"
version = "0.12.74"

[[MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "5a5bc6bf062f0f95e62d0fe0a2d99699fed82dd9"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.8"

[[ManualMemory]]
git-tree-sha1 = "9cb207b18148b2199db259adfa923b45593fe08e"
uuid = "d125e4d3-2237-4719-b19c-fa641b8a4667"
version = "0.1.6"

[[Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "Random", "Sockets"]
git-tree-sha1 = "1c38e51c3d08ef2278062ebceade0e46cefc96fe"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.0.3"

[[MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[Measures]]
git-tree-sha1 = "e498ddeee6f9fdb4551ce855a46f54dbd900245f"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.1"

[[Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[Mocking]]
deps = ["ExprTools"]
git-tree-sha1 = "748f6e1e4de814b101911e64cc12d83a6af66782"
uuid = "78c3b35d-d492-501b-9361-3d52fe80e533"
version = "0.7.2"

[[MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[MuladdMacro]]
git-tree-sha1 = "c6190f9a7fc5d9d5915ab29f2134421b12d24a68"
uuid = "46d2c3a1-f734-5fdb-9937-b9b9aeba4221"
version = "0.2.2"

[[MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "3927848ccebcc165952dc0d9ac9aa274a87bfe01"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "0.2.20"

[[NaNMath]]
git-tree-sha1 = "bfe47e760d60b82b66b61d2d44128b62e3a369fb"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "0.3.5"

[[NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[NonlinearSolve]]
deps = ["ArrayInterface", "FiniteDiff", "ForwardDiff", "IterativeSolvers", "LinearAlgebra", "RecursiveArrayTools", "RecursiveFactorization", "Reexport", "SciMLBase", "Setfield", "StaticArrays", "UnPack"]
git-tree-sha1 = "e9ffc92217b8709e0cf7b8808f6223a4a0936c95"
uuid = "8913a72c-1f9b-4ce2-8d82-65094dcecaec"
version = "0.3.11"

[[ODE]]
deps = ["Compat", "DiffEqBase", "LinearAlgebra", "Polynomials", "Reexport"]
git-tree-sha1 = "3a1b05ab7115ea428a71ef61c2f2b9597b6fb4b9"
uuid = "c030b06c-0b6d-57c2-b091-7029874bd033"
version = "2.13.0"

[[OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "c0e9e582987d36d5a61e650e6e543b9e44d9914b"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.10.7"

[[Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7937eda4681660b4d6aeeecc2f7e1c81c8ee4e2f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+0"

[[OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "15003dcb7d8db3c6c857fda14891a539a8f2705a"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.10+0"

[[OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[PCRE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b2a7af664e098055a7529ad1a900ded962bca488"
uuid = "2f80f16e-611a-54ab-bc61-aa92de5b98fc"
version = "8.44.0+0"

[[PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "4dd403333bcf0909341cfe57ec115152f937d7d8"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.1"

[[Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "2276ac65f1e236e0a6ea70baff3f62ad4c625345"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.2"

[[Parsers]]
deps = ["Dates"]
git-tree-sha1 = "438d35d2d95ae2c5e8780b330592b6de8494e779"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.0.3"

[[Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b4f5d02549a10e20780a24fce72bea96b6329e29"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.1+0"

[[Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[PlotThemes]]
deps = ["PlotUtils", "Requires", "Statistics"]
git-tree-sha1 = "a3a964ce9dc7898193536002a6dd892b1b5a6f1d"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "2.0.1"

[[PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "2537ed3c0ed5e03896927187f5f2ee6a4ab342db"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.0.14"

[[Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "GeometryBasics", "JSON", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "PlotThemes", "PlotUtils", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs"]
git-tree-sha1 = "4c2637482176b1c2fb99af4d83cb2ff0328fc33c"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.22.1"

[[PlutoUI]]
deps = ["Base64", "Dates", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "Suppressor"]
git-tree-sha1 = "44e225d5837e2a2345e69a1d1e01ac2443ff9fcb"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.9"

[[Polyester]]
deps = ["ArrayInterface", "BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "ManualMemory", "PolyesterWeave", "Requires", "Static", "StrideArraysCore", "ThreadingUtilities"]
git-tree-sha1 = "74d358e649e0450cb5d3ff54ca7c8d806ed62765"
uuid = "f517fe37-dbe3-4b94-8317-1923a5111588"
version = "0.5.1"

[[PolyesterWeave]]
deps = ["BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "Static", "ThreadingUtilities"]
git-tree-sha1 = "371a19bb801c1b420b29141750f3a34d6c6634b9"
uuid = "1d0040c9-8b98-4ee7-8388-3f51789ca0ad"
version = "0.1.0"

[[Polynomials]]
deps = ["Intervals", "LinearAlgebra", "MutableArithmetics", "RecipesBase"]
git-tree-sha1 = "0bbfdcd8cda81b8144de4be8a67f5717e959a005"
uuid = "f27b6e38-b328-58d1-80ce-0feddd5e7a45"
version = "2.0.14"

[[PreallocationTools]]
deps = ["ArrayInterface", "ForwardDiff", "LabelledArrays"]
git-tree-sha1 = "9e917b108c4aaf47e8606542325bd2ccbcac7ca4"
uuid = "d236fae5-4411-538c-8e31-a6e3d9e00b46"
version = "0.1.0"

[[Preferences]]
deps = ["TOML"]
git-tree-sha1 = "00cfd92944ca9c760982747e9a1d0d5d86ab1e5a"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.2.2"

[[Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[Qt5Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "ad368663a5e20dbb8d6dc2fddeefe4dae0781ae8"
uuid = "ea2cea3b-5b76-57ae-a6ef-0a8af62496e1"
version = "5.15.3+0"

[[QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "78aadffb3efd2155af139781b8a8df1ef279ea39"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.4.2"

[[REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[Random]]
deps = ["Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[RecipesBase]]
git-tree-sha1 = "44a75aa7a527910ee3d1751d1f0e4148698add9e"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.1.2"

[[RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "RecipesBase"]
git-tree-sha1 = "7ad0dfa8d03b7bcf8c597f59f5292801730c55b8"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.4.1"

[[RecursiveArrayTools]]
deps = ["ArrayInterface", "ChainRulesCore", "DocStringExtensions", "LinearAlgebra", "RecipesBase", "Requires", "StaticArrays", "Statistics", "ZygoteRules"]
git-tree-sha1 = "00bede2eb099dcc1ddc3f9ec02180c326b420ee2"
uuid = "731186ca-8d62-57ce-b412-fbd966d074cd"
version = "2.17.2"

[[RecursiveFactorization]]
deps = ["LinearAlgebra", "LoopVectorization", "Polyester", "StrideArraysCore", "TriangularSolve"]
git-tree-sha1 = "575c18c6b00ce409f75d96fefe33ebe01575457a"
uuid = "f2c3362d-daeb-58d1-803e-2bc74f2840b4"
version = "0.2.4"

[[Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "4036a3bd08ac7e968e27c203d45f5fff15020621"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.1.3"

[[Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "bf3188feca147ce108c76ad82c2792c57abe7b1f"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.0"

[[Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "68db32dff12bb6127bac73c209881191bf0efbb7"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.3.0+0"

[[SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[SIMDTypes]]
git-tree-sha1 = "330289636fb8107c5f32088d2741e9fd7a061a5c"
uuid = "94e857df-77ce-4151-89e5-788b33177be4"
version = "0.1.0"

[[SLEEFPirates]]
deps = ["IfElse", "Static", "VectorizationBase"]
git-tree-sha1 = "2e8150c7d2a14ac68537c7aac25faa6577aff046"
uuid = "476501e8-09a2-5ece-8869-fb82de89a1fa"
version = "0.6.27"

[[SciMLBase]]
deps = ["ArrayInterface", "CommonSolve", "ConstructionBase", "Distributed", "DocStringExtensions", "IteratorInterfaceExtensions", "LinearAlgebra", "Logging", "RecipesBase", "RecursiveArrayTools", "StaticArrays", "Statistics", "Tables", "TreeViews"]
git-tree-sha1 = "ff686e0c79dbe91767f4c1e44257621a5455b1c6"
uuid = "0bca4576-84f4-4d90-8ffe-ffa030f20462"
version = "1.18.7"

[[Scratch]]
deps = ["Dates"]
git-tree-sha1 = "0b4b7f1393cff97c33891da2a0bf69c6ed241fda"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.1.0"

[[Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "Requires"]
git-tree-sha1 = "fca29e68c5062722b5b4435594c3d1ba557072a3"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "0.7.1"

[[SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "b3363d7460f7d098ca0912c69b082f75625d7508"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.0.1"

[[SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[SpecialFunctions]]
deps = ["ChainRulesCore", "LogExpFunctions", "OpenSpecFun_jll"]
git-tree-sha1 = "a322a9493e49c5f3a10b50df3aedaf1cdb3244b7"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "1.6.1"

[[Static]]
deps = ["IfElse"]
git-tree-sha1 = "a8f30abc7c64a39d389680b74e749cf33f872a70"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "0.3.3"

[[StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "3240808c6d463ac46f1c1cd7638375cd22abbccb"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.2.12"

[[Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[StatsAPI]]
git-tree-sha1 = "1958272568dc176a1d881acb797beb909c785510"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.0.0"

[[StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "8cbbc098554648c84f79a463c9ff0fd277144b6c"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.10"

[[StatsFuns]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "46d7ccc7104860c38b11966dd1f72ff042f382e4"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "0.9.10"

[[StrideArraysCore]]
deps = ["ArrayInterface", "CloseOpenIntervals", "IfElse", "LayoutPointers", "ManualMemory", "Requires", "SIMDTypes", "Static", "ThreadingUtilities"]
git-tree-sha1 = "1258e25e171aec339866f283a11e7d75867e77d7"
uuid = "7792a7ef-975c-4747-a70f-980b88e8d1da"
version = "0.2.4"

[[StructArrays]]
deps = ["Adapt", "DataAPI", "StaticArrays", "Tables"]
git-tree-sha1 = "2ce41e0d042c60ecd131e9fb7154a3bfadbf50d3"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.3"

[[SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[Suppressor]]
git-tree-sha1 = "a819d77f31f83e5792a76081eee1ea6342ab8787"
uuid = "fd094767-a336-5f1f-9728-57cf17d0bbfb"
version = "0.2.0"

[[TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "TableTraits", "Test"]
git-tree-sha1 = "1162ce4a6c4b7e31e0e6b14486a6986951c73be9"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.5.2"

[[Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[ThreadingUtilities]]
deps = ["ManualMemory"]
git-tree-sha1 = "03013c6ae7f1824131b2ae2fc1d49793b51e8394"
uuid = "8290d209-cae3-49c0-8002-c8c24d57dab5"
version = "0.4.6"

[[TimeZones]]
deps = ["Dates", "Future", "LazyArtifacts", "Mocking", "Pkg", "Printf", "RecipesBase", "Serialization", "Unicode"]
git-tree-sha1 = "6c9040665b2da00d30143261aea22c7427aada1c"
uuid = "f269a46b-ccf7-5d73-abea-4c690281aa53"
version = "1.5.7"

[[TreeViews]]
deps = ["Test"]
git-tree-sha1 = "8d0d7a3fe2f30d6a7f833a5f19f7c7a5b396eae6"
uuid = "a2a6695c-b41b-5b7d-aed9-dbfdeacea5d7"
version = "0.3.0"

[[TriangularSolve]]
deps = ["CloseOpenIntervals", "IfElse", "LayoutPointers", "LinearAlgebra", "LoopVectorization", "Polyester", "Static", "VectorizationBase"]
git-tree-sha1 = "ed55426a514db35f58d36c3812aae89cfc057401"
uuid = "d5829a12-d9aa-46ab-831f-fb7c9ab06edf"
version = "0.1.6"

[[URIs]]
git-tree-sha1 = "97bbe755a53fe859669cd907f2d96aee8d2c1355"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.3.0"

[[UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[VectorizationBase]]
deps = ["ArrayInterface", "CPUSummary", "HostCPUFeatures", "Hwloc", "IfElse", "LayoutPointers", "Libdl", "LinearAlgebra", "SIMDTypes", "Static"]
git-tree-sha1 = "a5324cccb9ebab2e8bfc9bb8eb684394de2517e1"
uuid = "3d5dd08c-fd9d-11e8-17fa-ed2836048c2f"
version = "0.21.9"

[[Wayland_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "3e61f0b86f90dacb0bc0e73a0c5a83f6a8636e23"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.19.0+0"

[[Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll"]
git-tree-sha1 = "2839f1c1296940218e35df0bbb220f2a79686670"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.18.0+4"

[[XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "1acf5bdf07aa0907e0a37d3718bb88d4b687b74a"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.9.12+0"

[[XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "5be649d550f3f4b95308bf0183b82e2582876527"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.6.9+4"

[[Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4e490d5c960c314f33885790ed410ff3a94ce67e"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.9+4"

[[Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fe47bd2247248125c428978740e18a681372dd4"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.3+4"

[[Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6783737e45d3c59a4a4c4091f5f88cdcf0908cbb"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.0+3"

[[Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "daf17f441228e7a3833846cd048892861cff16d6"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.13.0+3"

[[Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "926af861744212db0eb001d9e40b5d16292080b2"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.0+4"

[[Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "4bcbf660f6c2e714f87e960a171b119d06ee163b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.2+4"

[[Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "5c8424f8a67c3f2209646d4425f3d415fee5931d"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.27.0+4"

[[Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "79c31e7844f6ecf779705fbc12146eb190b7d845"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.4.0+3"

[[Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "cc4bf3fdde8b7e3e9fa0351bdeedba1cf3b7f6e6"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.0+0"

[[ZygoteRules]]
deps = ["MacroTools"]
git-tree-sha1 = "9e7a1e8ca60b742e508a315c17eef5211e7fbfd7"
uuid = "700de1a5-db45-46bc-99cf-38207098b444"
version = "0.2.1"

[[libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "c45f4e40e7aafe9d086379e5578947ec8b95a8fb"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+0"

[[nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"

[[x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "ece2350174195bb31de1a63bea3a41ae1aa593b6"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "0.9.1+5"
"""

# ╔═╡ Cell order:
# ╠═0d4fbf25-a629-416a-bb67-b087144a5862
# ╠═ac1342d7-a80c-4942-b424-6cb44904a07e
# ╟─b5d87000-5af5-11eb-0594-57573e58e809
# ╟─68a29859-4ba6-46ed-97ef-fdffec416003
# ╠═fce92627-56f1-482e-983a-083b7281fbcb
# ╟─934a054a-f31d-497e-bcf8-295dac6db83b
# ╠═4f30e023-2c8d-42b7-a124-c240e5464769
# ╠═6f9f8bbd-1318-4d66-bdd6-6d02f778f8b8
# ╟─5e6c51c1-c7a9-4542-89fd-9c7b0846437c
# ╟─54c26825-46b9-44f6-ae65-31f4be69bf3b
# ╠═22f0b250-0b37-4593-bdea-4bd13c19f2bf
# ╠═bfc718ab-e392-4fb8-82c6-0c889b1977e9
# ╟─12ceb8e5-f976-4639-9175-047e84321c85
# ╠═1e36310b-6f42-44a1-8b04-e82391cb33aa
# ╟─8899d0fd-e0d9-424f-98e4-5547d029a36a
# ╟─2d44d600-5af7-11eb-0958-b3eab7ae938d
# ╟─8b933ce6-7fea-4190-940a-bd7720fc93be
# ╠═cb8721aa-7f4b-4672-a802-da84d9c67b4a
# ╟─ad54ce7f-8caa-46f8-9a57-f3710ca89f40
# ╠═4adf6672-5ae1-42e1-85e9-e522d6069fab
# ╠═a262689c-74f3-46bb-bd55-1aa7eec5379f
# ╠═481e7da8-95f2-4972-863c-fcd7a79cf417
# ╟─d52999d7-2644-4793-a36a-413ad8011e09
# ╠═1a0815c9-900d-4ba5-855a-7fa582e8d498
# ╠═d8e15baa-34a9-447d-90d2-8dfe8065ee1b
# ╠═f35f1bb2-bc2c-4bb7-9dc9-367bc5a6c15d
# ╠═e50cf9c2-972d-4916-a3ca-6c20a0c73519
# ╠═f4d37a14-e6ef-4400-bb30-a3a4405bb675
# ╠═a7841a3d-c7e4-4c09-9012-a3890c9e1195
# ╟─471ddc14-d7ba-4d94-a902-9ca459be5ce3
# ╟─705ebcb9-bfa6-4357-8610-27c43c0c3f96
# ╠═a6cd2db9-f4b6-46cb-aaa0-59d53430056c
# ╠═78eae1ea-5ce9-4bda-81d3-4b056ee7df1e
# ╠═00796b0e-f2d5-4646-af33-06d81d13fc0c
# ╠═0b5ad338-3628-477e-bd0b-8cdda13d345d
# ╟─5308aa28-6da6-4943-b393-c08592b98c27
# ╠═00f084ce-a5ff-4e18-913f-9d14ae346bd9
# ╠═0ab04df9-d9a0-406a-b45f-299a0212fd9b
# ╠═e580ec4b-73f9-4522-9cb3-4d24ebe1e027
# ╟─3ada3bbf-cc3b-46b4-adbd-c3685568781a
# ╠═39d0b178-70bd-4815-bee8-3738bdaea250
# ╠═12c6cb02-6bc0-4558-bc10-9ca1d8cd029a
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
