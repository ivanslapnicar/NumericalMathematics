### A Pluto.jl notebook ###
# v0.12.18

using Markdown
using InteractiveUtils

# ╔═╡ 0d4fbf25-a629-416a-bb67-b087144a5862
using Plots

# ╔═╡ 69cef7ae-8c14-46fc-8583-942c1bdec167
using ODE

# ╔═╡ b5d87000-5af5-11eb-0594-57573e58e809
md"""
# Differential Equations


Consider first order differential equation

$$
\frac{d}{dx} y(x)=f(x,y(x)),$$

with the given __initial condition__

$$
y(x_0)=y_0.$$

__Theorem.__ If the functions $f$ and $\displaystyle \frac{\partial f}{\partial y}$ are continuous in some neighbourhood of the point $(x_0,y_0)$, then, for some $\varepsilon > 0$, the given initial value pbolem has unique solution $y(x)$ on the interval $\displaystyle [x_{0}-\varepsilon ,x_{0}+\varepsilon ]$.

__Remark.__ The independent variable is often denoted by $t$ (time). 
"""

# ╔═╡ 68a29859-4ba6-46ed-97ef-fdffec416003
md"""
## Euler's method

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
### Example 1

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
### Example 2

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
## Runge-Kutta methods

Euler's method is a first-order method, and it not accurate enough. Thus, in practice we use methods of higher order which approximate $y(x)$ in the point $x_{k+1}$ using values of the function $f(x,y)$ in several points from the interval 

$$[x_k,x_{k+1}]\equiv[x_k,x_k+h].$$

### Heun's method

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

### Standard Runge-Kutta method

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
### Example 3

Let us solve problems from Examples 1 and 2. For Example 1 the numerical solution graphically overlaps the exact solution. 
For Example 2 the solution using the function `RungeKutta4()` is an order of magnitude more accurate than the solution obtained using the function  red `Euler()`.  
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
### Existing routines

Most programming languages have built-in routines for numerical solution of ordinary differential equations.
For example, 

* Matlab has commands `ode*` (see [Matlab, Ordinary Differential Equations](https://www.mathworks.com/help/matlab/ordinary-differential-equations.html?searchHighlight=ordinary%20differential&s_tid=srchtitle)), a 
* Julia has the package [ODE.jl](https://github.com/JuliaODE/ODE.jl).

Standard Runge-Kutta method is implemented in the function `ode4()`, and Heun's method is implemented in the function 
`ODE.ode2_heun()`. 

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
# Let us solve the problem of Example 2.
# Computed values yₖ are the first elements of the output vectors.
yode4=ode4(f₂,0.0,range(0,stop=1,length=21))[2]

# ╔═╡ f4d37a14-e6ef-4400-bb30-a3a4405bb675
yode2=ODE.ode2_heun(f₂,0.0,range(0,stop=1,length=21))[2];

# ╔═╡ a7841a3d-c7e4-4c09-9012-a3890c9e1195
# Compare solutions
[yRK4 yode4 yRK4-yode4 yode2 yode4-yode2]

# ╔═╡ 471ddc14-d7ba-4d94-a902-9ca459be5ce3
md"""
## Systems of differential equations

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
### Lotka-Volterra equations

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
### Scaled Lotka-Volterra equations

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
## Differential equations of higher order

Using substitutions, differential equation of higher order can be reduced to a system of differential equations of the first order.

### Example 4

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
sqrt(sum((solution₅.(x₅)-YRK4₅).^2))

# ╔═╡ ae54d9f1-bed6-4f56-83a2-832e819b5068


# ╔═╡ Cell order:
# ╟─b5d87000-5af5-11eb-0594-57573e58e809
# ╟─68a29859-4ba6-46ed-97ef-fdffec416003
# ╠═0d4fbf25-a629-416a-bb67-b087144a5862
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
# ╠═69cef7ae-8c14-46fc-8583-942c1bdec167
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
# ╠═ae54d9f1-bed6-4f56-83a2-832e819b5068
