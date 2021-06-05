### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# ╔═╡ 271cba4c-ed08-497e-8142-c99aea337cd8
# On your computer, comment this cell
begin
	import Pkg
	Pkg.activate(mktempdir())
    Pkg.add([
        Pkg.PackageSpec(name="PlutoUI"),
		Pkg.PackageSpec(name="ForwardDiff"),
		Pkg.PackageSpec(name="Plots"),
		Pkg.PackageSpec(name="NLsolve"),
		Pkg.PackageSpec(name="Optim")
    ])
end

# ╔═╡ 57657193-e72f-432f-965b-13f36c104874
begin
	using PlutoUI, ForwardDiff, Plots, LinearAlgebra, NLsolve, Optim
	plotly()
end

# ╔═╡ 63690bb6-5749-429d-8c39-dd9c727bc938
TableOfContents(title="📚 Table of Contents", aside=true)

# ╔═╡ 94da0d5c-6dfe-487a-9711-9d101f20e97a
md"""
# Systems of Nonlinear Equations


__Problem.__ Let us find a solution $\xi=(\xi_1,\xi_2,\ldots,\xi_n)$ of the system of $n$ equations 

$$\begin{aligned}
f_1(x)&=0,\cr
f_2(x)&=0,\cr
&\vdots \cr
f_n(x)&=0,
\end{aligned}$$

and $n$ unknowns $x=(x_1,x_2,\ldots,x_n)$. Denoting $f=(f_1,f_2,\ldots,f_n)^T$, we can write this system as

$$
f(x)=0.$$

We shall describe the __Newton's method__ and three __quasi-Newton__ methods:

2. __Broyden's__ method,
3. __Davidon-Fletcher-Powell__ method, and 
3. __Broyden-Fletcher-Goldfarb-Schano__ method.

Given starting approximation $x^{(0)}$, all methods generate sequence of points $x^{(n)}$ which, under certain condditions, converges towards the solution $\xi$. 

__Remark.__ For detailed description of the methods and examples see [Numerical Methods for Unconstrained Optimization and Nonlinear Equations](https://epubs.siam.org/doi/book/10.1137/1.9781611971200?mobileUi=0) or [Numerička matematika, section 4.4](http://www.mathos.unios.hr/pim/Materijali/Num.pdf).
"""

# ╔═╡ 125c9230-cc68-493f-92c5-8a83a78b5863
md"""
# Newton's method

__Jacobian__ or __Jacobi matrix__ of functions $f$ in the point $x$ is the matrix of first partial derivatives

$$
J(f,x)=\begin{pmatrix} \displaystyle\frac{\partial f_1(x)}{\partial x_1} & \displaystyle\frac{\partial f_1(x)}{\partial x_2} & \cdots &
\displaystyle\frac{\partial f_1(x)}{\partial x_n} \\
\displaystyle\frac{\partial f_2(x)}{\partial x_1} & \displaystyle\frac{\partial f_2(x)}{\partial x_2} & \cdots &
\displaystyle\frac{\partial f_2(x)}{\partial x_n} \\
\vdots & \vdots & \ddots & \vdots \\
\displaystyle\frac{\partial f_n(x)}{\partial x_1} & \displaystyle\frac{\partial f_n(x)}{\partial x_2} & \cdots &
\displaystyle\frac{\partial f_n(x)}{\partial x_n} 
\end{pmatrix}.$$

Given starting approximation $x^{(0)}$, we compute the sequence of points

$$
x^{(k+1)}=x^{(k)}-s^{(k)}, \quad k=0,1,2,\ldots,$$

where $s^{(k)}$ is the solution of the system of linear equations

$$
J\big(f,x^{(k)}\big)\cdot s=f\big(x^{(k)}\big).$$

Jacobians are computed using the package [`ForwardDiff.jl`](http://www.juliadiff.org/ForwardDiff.jl/perf_diff.html#derivatives). Plotting is done using the package `Plots.jl`.
"""

# ╔═╡ 77095a05-eb64-47e2-a2e0-5053a4fbc837
function Newton(f::Function,J::Function,x::Vector{T},ϵ::Float64=1e-10) where T
    iter=0
    s=ones(T,size(x))
    ξ=x
    while norm(s)>ϵ && iter<100
        s=J(x)\f(x)
        ξ=x-s
        iter+=1
        x=ξ
    end
    ξ,iter
end

# ╔═╡ ed8971fd-1420-41af-b736-e052439a2c65
md"""
## Examples

### $2\times 2$ nonlinear system

The solutions of the system

$$\begin{aligned}
2(x+y)^2+(x-y)^2-8&=0\\
5x^2+(y-3)^2-9&=0
\end{aligned}$$

are the points $T_1=(1,1)$ and $T_2\approx(-1.18,1.59)$.
"""

# ╔═╡ 96eb13a8-ca15-447b-a252-76d481b30912
# Vector function
f₁(x)=[2(x[1]+x[2])^2+(x[1]-x[2])^2-8,5*x[1]^2+(x[2]-3)^2-9]

# ╔═╡ 7366924f-29e6-46a6-a5f1-e1940daa15dd
f₁((1.0,2))

# ╔═╡ d7930d33-99f3-4551-8258-500bbde17aff
md"""
Let us plot the functions and the contours in order to approximately locate the zeros:
"""

# ╔═╡ 46060676-7192-4e7d-9038-9492524c1485
begin
	# Number of points
	m=101
	X=range(-2,stop=3,length=m)
	Y=range(-2,stop=3,length=m)
	# First applicate
	surface(X,Y,(x,y)->f₁([x,y])[1],xlabel="x",ylabel="y",colorbar=false)
	# Second applicate
	surface!(X,Y,(x,y)->f₁([x,y])[2],seriescolor=:blues)
end

# ╔═╡ 0234d953-b83d-48a1-9c32-a0e2b8bec20a
begin
	# Locate the solutions using contour plot
	contour(X,Y,(x,y)->f₁([x,y])[1],contour_labels=true)
	contour!(X,Y,(x,y)->f₁([x,y])[2],contour_labels=true)
end

# ╔═╡ 7d4c9549-1233-4e03-ace4-1a1c87147036
# Clearer picture
contour!(clims=(0,0.01),xlabel="x",ylabel="y",colorbar=:none)

# ╔═╡ e2b8b699-39e0-408a-8165-2f606ea25740
md"""
We see that the approximate zeros are $T_1=(-1.2,1.5)$ i $T_2=(1,1)$. Moreover, $T_2$ is exactly equal $(1,1)$ (one iteration in the third computation). Furthermore, the method does not always converge (fourth computation).   
"""

# ╔═╡ 08505d4c-b799-4833-af0d-4efea3f925c1
J₁(x)=ForwardDiff.jacobian(f₁,x)

# ╔═╡ a5631cd4-c6b7-4a7b-9411-adf83f7c15b3
# For example
J₁([1.0,2])

# ╔═╡ 68d9ec70-cf1d-4260-8a66-5cce42c42eea
Newton(f₁,J₁,[-1.0,0.0]), Newton(f₁,J₁,[0.5,1.1]), 
Newton(f₁,J₁,[1.0,1.0]), Newton(f₁,J₁,[0.0,0.0])

# ╔═╡ 3cfa6d61-75f0-4897-97d1-b94c8de2458e
md"""
### $2\times 2$ nonlinear system

The solutions of the system

$$\begin{aligned}
x_1^2-x_2^2-2&=0\\
e^{x_1-1}+x_2^3-2&=0
\end{aligned}$$

are the points $T_1=(1,1)$ and $T_2\approx (-0.71,1.22)$.
"""

# ╔═╡ 2a1b5d9b-672f-4e00-8620-daec0196d6b7
begin
	f₂(x)=[x[1]^2+x[2]^2-2,exp(x[1]-1)+x[2]^3-2]
	contour(X,Y,(x,y)->f₁([x,y])[1],contour_labels=true)
	contour!(X,Y,(x,y)->f₁([x,y])[2],contour_labels=true)
	contour!(clims=(0,0.01),xlabel="x",ylabel="y",colorbar=:none)
end

# ╔═╡ 6c8b54c2-d4ae-4ada-b0d9-cd83eafa9cc1
begin
	J₂(x)=ForwardDiff.jacobian(f₂,x)
	Newton(f₂,J₂,[-1.0,1]), Newton(f₂,J₂,[0.8,1.2])
end

# ╔═╡ 2f2a6ae5-755d-4f65-9c16-4db2e1aa48f4
md"""
### $3\times 3$ nonlinear system

Solve $f(x)=0$, where

$$
f(x)=\begin{bmatrix}x_1 \\ x_2^2-x_2 \\ e^{x_3}-1 \end{bmatrix}.$$

The exact solutions are $T_1=(0,0,0)$ and $T_2=(0,-1,0)$. We shall compute the solutions using several starting points.
"""

# ╔═╡ c6f98f59-312d-498d-9854-4c78c969011f
begin
	f₃(x)=[x[1],x[2]^2+x[2],exp(x[3])-1]
	J₃(x)=ForwardDiff.jacobian(f₃,x)
end

# ╔═╡ ee2d1d89-d399-4d1f-8d08-c99e372d5537
Newton(f₃,J₃,[-1.0,1.0,0.0]),Newton(f₃,J₃,[1.0,1,1]),
Newton(f₃,J₃,[-1.0,1,-10]),Newton(f₃,J₃,[0.5,-1.5,0])

# ╔═╡ e195e6ad-f3ea-4a74-a91a-d16d0384a054
md"""
### Rosenbrock parabolic valley

Given is the function 

$$
f(x)=100\,(x_2-x_1)^2+(1-x_1)^2.$$

We want to find possible extremal points, thet is, we want to solve the equation

$$
\mathop{\mathrm{grad}} f(x)=0.$$
"""

# ╔═╡ 6e1b5ed8-3d23-4aad-b140-1567c8dc7336
f₄(x)=100(x[2]-x[1]^2)^2+(1-x[1])^2

# ╔═╡ d9704e21-9e89-4b5d-bb4e-7dd24d6b506a
# Plot the function using X and Y from Example 1
surface(X,Y,(x,y)->f₄([x,y]), seriescolor=:blues, xlabel="x", ylabel="y")

# ╔═╡ 008c4d22-4c2f-4087-8be7-47599600de47
begin
	# This function is very demanding w.r.t. finding extremal points
	g₄(x)=ForwardDiff.gradient(f₄,collect(x))
	contour(X,Y,(x,y)->g₄([x,y])[1],contour_labels=true)
	contour!(X,Y,(x,y)->g₄([x,y])[2],contour_labels=true)
	contour!(clims=(-0.5,0.5),xlabel="x",ylabel="y",colorbar=:none)
end

# ╔═╡ ac8f76d2-9801-4f64-8add-278c8ada3135
md"""
By observing contours, we conclude that this example is numerically demanding, while it is easy to see analytically that the only zero is $T=(1,1)$.

Here the vector function is given as the gradient of the scalar function $f$, so the Jacobi matrix is computed using function `FowardDiff.hessian()` which approximates the matrix of second partial derivatives of the scalar function $f$.  
"""

# ╔═╡ 34576988-89ce-431b-9e1c-78b1a98e002e
Newton(g₄,x->ForwardDiff.hessian(f₄,x),[-1.0,2.0])

# ╔═╡ 9625fb8d-5e0e-4e07-bc49-b4cb9ab259ef
md"""
### Least-squares fitting of sum of exponentials

Let

$$
f(x)=\sum_{i=1}^{11} \bigg(x_3 \cdot \exp\bigg(-\frac{(t_i-x_1)^2}{x_2}\bigg)-y_i\bigg)^2,$$

where the pairs  $(t_i,y_i)$ are given by the table:

$$
\begin{array}{c|c|c|c|c|c|c|c|c|c|c|c}
 i & 1 & 2 & 3 & 4 & 5 & 6 & 7 & 8 & 9 & 10 & 11 \\ \hline
t_i & 0 & 1 & 2 & 3 & 4 & 5 & 6 & 7 & 8 & 9 &10 \\
 y_i & 0.001 & .01 & .04 & .12 & .21 & .25 & .21 & .12 & .04 & .01 & .001 
\end{array}$$

We want to solve the equation

$$
\mathop{\mathrm{grad}} f(x)=0.$$

Unlike previous examples, where the condition number is $\kappa(J)=O(10)$ in the first three examples and $\kappa(J)=O(1000)$ in the previous example, here $\kappa(J)>O(10^6)$ so the method is inaccurate and does not converge towards the exact solution $x=(4.93,2.62,0.28)$.
"""

# ╔═╡ 33ca0d51-92fe-4313-ba67-7792af17b11e
begin
	t=collect(0:10)
	y=[0.001,0.01,0.04,0.12,0.21,0.25,0.21,0.12,0.04,0.01,0.001]
	f₅(x)=sum([( x[3]*exp(-((t[i]-x[1])^2/x[2]))-y[1])^2 for i=1:11])
end

# ╔═╡ e5d2cee7-c2d4-45a5-901c-206be7d2ff58
begin
	# Starting point is very near the solution
	x₀=[4.9,2.63,0.28]
	f₅(x₀)
	g₅(x)=ForwardDiff.gradient(f₅,x)
	J₅(x)=ForwardDiff.hessian(f₅,x)
	f₅(x₀), g₅(x₀), cond(J₅(x₀))
end

# ╔═╡ ab87421c-e498-44ff-af23-13ea14a0c5d3
x₅,iter₅=Newton(g₅,J₅,x₀,1e-8)

# ╔═╡ d9df727c-a96f-4af4-ab5a-60715662adce
g₅(x₅)

# ╔═╡ 82fffd1c-9ec0-459c-b033-d926641268dc
Newton(g₅,J₅,[4.9,2.62,0.28],1e-8)

# ╔═╡ e823348a-1ea3-45b7-b1ad-d51096f224dc
md"""
# Broyden's method

Given starting approximation $x_0$ and matrix $B_0$, for each $k=0,1,2,\ldots$, we compute:

$$\begin{aligned}
B_k \cdot s_k & = -f(x_k) \quad \textrm{(sustav)}\\
x_{k+1}&=x_{k}+s_k\\
y_k&=f(x_{k+1})-f(x_{k})\\
B_{k+1}&=B_k+\frac{(y_k-B_ks_k)s_k^T}{s_k\cdot s_k}
\end{aligned}$$

In this manner, we avoid the computation of the Jacobi matrix in each step. We can take $B_0=J(x_0)$, but also some other matrix.
"""

# ╔═╡ 73ad34ac-e4fc-4e35-9d96-25969b37a4d6
function Broyden(f::Function,B::Matrix,x::Vector{T},ϵ::Float64=1e-10) where T
    iter=0
    s=ones(T,length(x))
    ξ=x
    while norm(s)>ϵ && iter<100
        s=-(B\f(x))
        ξ=x+s
        y=f(ξ)-f(x)
        B=B+(y-B*s)*(s/(s⋅s))'
        x=ξ
        iter+=1
    end
    ξ,iter
end

# ╔═╡ 19dc36a5-6d13-4d61-a9b4-863e7fbe60fe
md"
## Examples
"

# ╔═╡ 6059b96d-8ba0-475c-9929-8cc613355562
# Example 1
Broyden(f₁,J₁([-1.0,0.0]),[-1.0,0.0]), 
Broyden(f₁,J₁([1.0,1.5]),[1.0,1.5])

# ╔═╡ 64d397ac-ed6f-426c-9139-af21222723da
begin
	# Explain behaviour the method when we set B₀=I
	eye(n)=Matrix{Float64}(I,n,n)
	Broyden(f₁,eye(2),[-1.0,0.0]), Broyden(f₁,eye(2),[1.0,1.5]),
	Broyden(f₁,eye(2),[-1,1.5])
end

# ╔═╡ ba09dd4f-b114-45bb-9b6a-41c2ed9dbec5
begin
	# Example 2
	x0=[-1.0,1]
	x1=[0.8,1.2]
	Broyden(f₂,J₂([-1.0,1]),[-1.0,1]), Broyden(f₂,J₂([0.8,1.2]),[0.8,1.2])
end

# ╔═╡ 52506f6c-069a-4883-8160-1cf415e80c36
# Example 3
Broyden(f₃,J₃([-1.0,1,0]),[-1.0,1,0]), Broyden(f₃,J₃([0.5,-1.5,0]),[0.5,-1.5,0])

# ╔═╡ c0fb2bfc-7221-428d-b13d-ebb483bd0c1c
# Example 4
Broyden(g₄,(x->ForwardDiff.hessian(f₄,x))([-1.0,2]),[-1.0,2]), # ali
Broyden(g₄,(x->ForwardDiff.hessian(f₄,x))([1,2.0]),[-1.0,2]),
Broyden(g₄,(x->ForwardDiff.hessian(f₄,x))([0.8,0.5]),[0.8,0.5])

# ╔═╡ 1fe960aa-f12a-4d3e-994a-1a5d5c10f29e
begin
	# Example 5
	x₆,iter₆=Broyden(g₅,(x->ForwardDiff.hessian(f₅,x))([4.9,2.6,0.2]), [4.9,2.6,0.2])
end

# ╔═╡ b4e246f8-72b9-4042-8160-a1eda2084e45
g₅(x₆)

# ╔═╡ 11259db5-f754-47c7-9e92-600adf6896f8
md"""
# Davidon-Fletcher-Powell (DFP) method

DFP is an optimization method which finds extremal points of the function $F:\mathbb{R}^n \to \mathbb{R}$, in which case $f(x)=\mathop{\mathrm{grad}}F(x)$.

Given initial approximation $x_0$ and matrix $H_0$, for $k=0,1,2,\ldots$, we compute:

$$\begin{aligned}
s_k&=-H_k f(x_k)\\
\beta_k&=\mathop{\mathrm{arg\ min}}_\beta F(x_{k}+\beta s_k) \\
s_k&=\beta_k s_k\\
x_{k+1}&=x_{k}+s_k \\
y_k&=f(x_{k+1})-f(x_{k})\\
H_{k+1}&=H_k+ \frac{s_k s_k^T}{y_k\cdot s_k}-\frac{H_k y_k y_k^T H_k}{y_k\cdot (H_k y_k)}.
\end{aligned}$$

We can take $H_0=I$. Notice that the iteration step does not require solving a system linear equations. Instead, all updates are performed using $O(n^2)$ operations.

The one-dimensional minimzation along the line $x_{k}+\beta s_k$ is preformed by finding zeros of the directional derivative using bisection.
"""

# ╔═╡ ff33e971-4871-4533-95ef-4fd9b14f32d1
function Bisection(f::Function,a::T,b::T,ϵ::Float64=1e-10) where T
    fa=f(a)
    fb=f(b)
    x=T
    fx=T
    if fa*fb>zero(T)
        # return "Incorrect interval"
        if abs(fa)>abs(fb)
            return b,fb,0
        else
            return a,fa,0
        end
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
    x,fx,iter
end

# ╔═╡ f29daf79-022a-4f7d-8616-6cac0e7ea9e8
function DFP(f::Function,H::Matrix,x::Vector{T},ϵ::Float64=1e-10) where T
    iter=0
    s=ones(T,length(x))
    ξ=x
    while norm(s)>ϵ && iter<50
        s=-H*f(x)
        s0=s/norm(s)
        F(ζ)=f(x+ζ*s)⋅s0
        β,fx,iterb=Bisection(F,0.0,1.0,10*eps())
        s*=β
        ξ=x+s
        y=f(ξ)-f(x)
        z=H*y
        H=H+(s/(y⋅s))*s'-(z/(y⋅z))*z'
        x=ξ
        iter+=1
    end
    ξ,iter
end

# ╔═╡ 7f7672da-787f-4628-83aa-e9297f7da322
md"""
### Minimum of a function of $2$ variables

Let us find extremal points of the function 

$$
f(x,y)=(x+2y-7)^2+(2x+y-5)^2.$$

The function has minimum at $(1,3)$.
"""

# ╔═╡ 8f4fc906-99a7-4e06-8620-83cf427da8b9
f₆(x) = (x[1] + 2*x[2]-7)^2 + (2*x[1] + x[2]-5)^2

# ╔═╡ b687e934-69d9-4608-b73f-0f3676c210ec
f₆([1,2])

# ╔═╡ 7ebf3fb9-85c7-4ea0-a18a-a4617375082d
g₆(x)=ForwardDiff.gradient(f₆,x)

# ╔═╡ 1f6187a2-41a2-43e0-9cd3-322a069222aa
DFP(g₆,eye(2),[0.8,2.7],eps())

# ╔═╡ d1e5bfe4-1ac9-4dbe-aafb-f81822a72d3a
# Example 4
DFP(g₄,eye(2),[0.9,1.1])

# ╔═╡ 02ef6f4f-3259-447e-878a-b116a3d0eb4a
# Example 5
DFP(g₅,eye(3),[4.9,2.6,0.2])

# ╔═╡ f7483843-9cc9-4799-ad56-1e74bd0a08e7
md"""
# Broyden-Fletcher-Goldfarb-Schano (BFGS) method

BFGS is an optimization method for finding extremal points of the function $F:\mathbb{R}^n \to \mathbb{R}$, in which case $f(x)=\mathop{\mathrm{grad}}F(x)$.

The method is similar to the DFP method, with somewhat better convergence properties.

Let $F:\mathbb{R}^n \to \mathbb{R}$, whose minmum we seek, and let $f(x)=\mathop{\mathrm{grad}} F(x)$.

Given initial approximation $x_0$ and matrix $H_0$, for $k=0,1,2,\ldots$, we compute:

$$\begin{aligned}
s_k&=-H_k f(x_k)\\
\beta_k&=\mathop{\mathrm{arg\ min}}F(x_{k}+\beta_k s_k) \\
s_k&=\beta_k s_k\\
x_{k+1}&=x_{k}+s_k \\
y_k&=f(x_{k+1})-f(x_{k})\\
H_{k+1}&=\bigg(I-\frac{s_k y_k^T}{y_k\cdot s_k}\bigg)H_k
\bigg( I-\frac{y_k s_k^T}{y_k\cdot s_k}\bigg)+\frac{s_k s_k^T}{y_k\cdot s_k}.
\end{aligned}$$

We can take $H_0=I$. Notice that the iteration step does not require solving a system linear equations. Instead, all updates are performed using $O(n^2)$ operations.

The one-dimensional minimzation along the line $x_{k}+\beta s_k$ is preformed by finding zeros of the directional derivative using bisection.
"""

# ╔═╡ b985704f-8f07-43b6-baa9-6b3d4541de27
function BFGS(f::Function,H::Matrix,x::Vector{T},ϵ::Float64=1e-10) where T
    iter=0
    s=ones(T,length(x))
    ξ=x
    while norm(s)>ϵ && iter<50
        s=-H*f(x)
        s0=s/norm(s)
        F(ζ)=f(x+ζ*s)⋅s0
        β,fx,iterb=Bisection(F,0.0,1.0,10*eps())
        s*=β
        ξ=x+s
        y=f(ξ)-f(x)
        z=H*y
        α=y⋅s
        s1=s/α
        H=H-s1*z'-z*s1'+s1*(y⋅z)*s1'+s1*s'
        x=ξ
        iter+=1
    end
    ξ,iter
end

# ╔═╡ 31a1b2aa-7b5b-46b0-82ec-f75111b72a0d
BFGS(g₆,eye(2),[0.8,2.7],eps())

# ╔═╡ a0e4f6af-d4ab-4b04-9863-49de95eec947
# Example 4
BFGS(g₄,eye(2),[0.9,1.1])

# ╔═╡ c20593d7-885f-4f36-984c-ecded0884ecc
# Example 5
BFGS(g₅,eye(3),[4.9,2.6,0.2])

# ╔═╡ cef37fa4-df23-4267-bde4-c1fe6ddade23
md"""
# Julia packages

Previous programs are simple illustrative implementations of the mentioned algorithms.
Julia package [NLsolve.jl](https://github.com/JuliaNLSolvers/NLsolve.jl) is used for solving systems of linear equations and the package [Optim.jl](https://github.com/JuliaNLSolvers/Optim.jl) is used for nonlinear optimization.

## Examples
"""

# ╔═╡ 92d27828-c97b-41d0-85e9-0d3b1738dfed
# Example 1
function f₁!(fvec,x)
    fvec[1] = 2(x[1]+x[2])^2+(x[1]-x[2])^2-8
    fvec[2] = 5*x[1]^2+(x[2]-3)^2-9
end

# ╔═╡ 6055ba96-dd58-420b-854e-277b60a3c1e7
nlsolve(f₁!,[-1.0,0])

# ╔═╡ ca176d00-70f9-42b4-a3f4-fd5d37c9f31a
nlsolve(f₁!,[0.5,1.1])

# ╔═╡ a8080f5e-6797-482e-84ca-7ddfcd40264d
begin
	# Example 2
	function f₂!(fvec,x)
	    fvec[1] = x[1]^2+x[2]^2-2
	    fvec[2] = exp(x[1]-1)+x[2]^3-2
	end
	nlsolve(f₂!,[-1.0,1]), nlsolve(f₂!,[0.8,1.2])
end

# ╔═╡ 2a644f10-2683-4ea8-bd25-95343abb7e48
begin
	# Example 3
	function f₃!(fvec,x)
	    fvec[1] = x[1]
	    fvec[2] = x[2]^2+x[2]
	    fvec[3] = exp(x[3])-1
	end
	nlsolve(f₃!,[-1.0,1.0,0.0]), nlsolve(f₃!,[1.0,1,1]),
	nlsolve(f₃!,[-1.0,1,-10]), nlsolve(f₃!,[0.5,-1.5,0])
end

# ╔═╡ bc477197-2ac5-4798-a62b-62f89791fb94
# Example 4
optimize(f₄,[-1.0,2],Optim.BFGS())

# ╔═╡ 90043461-3f6a-45e6-a661-57c8df3a2a93
# Example 5 - again there is no convergence
optimize(f₅,[4.9,2.6,0.2],Optim.BFGS())

# ╔═╡ Cell order:
# ╠═271cba4c-ed08-497e-8142-c99aea337cd8
# ╠═57657193-e72f-432f-965b-13f36c104874
# ╠═63690bb6-5749-429d-8c39-dd9c727bc938
# ╟─94da0d5c-6dfe-487a-9711-9d101f20e97a
# ╟─125c9230-cc68-493f-92c5-8a83a78b5863
# ╠═77095a05-eb64-47e2-a2e0-5053a4fbc837
# ╟─ed8971fd-1420-41af-b736-e052439a2c65
# ╠═96eb13a8-ca15-447b-a252-76d481b30912
# ╠═7366924f-29e6-46a6-a5f1-e1940daa15dd
# ╟─d7930d33-99f3-4551-8258-500bbde17aff
# ╠═46060676-7192-4e7d-9038-9492524c1485
# ╠═0234d953-b83d-48a1-9c32-a0e2b8bec20a
# ╠═7d4c9549-1233-4e03-ace4-1a1c87147036
# ╟─e2b8b699-39e0-408a-8165-2f606ea25740
# ╠═08505d4c-b799-4833-af0d-4efea3f925c1
# ╠═a5631cd4-c6b7-4a7b-9411-adf83f7c15b3
# ╠═68d9ec70-cf1d-4260-8a66-5cce42c42eea
# ╟─3cfa6d61-75f0-4897-97d1-b94c8de2458e
# ╠═2a1b5d9b-672f-4e00-8620-daec0196d6b7
# ╠═6c8b54c2-d4ae-4ada-b0d9-cd83eafa9cc1
# ╟─2f2a6ae5-755d-4f65-9c16-4db2e1aa48f4
# ╠═c6f98f59-312d-498d-9854-4c78c969011f
# ╠═ee2d1d89-d399-4d1f-8d08-c99e372d5537
# ╟─e195e6ad-f3ea-4a74-a91a-d16d0384a054
# ╠═6e1b5ed8-3d23-4aad-b140-1567c8dc7336
# ╠═d9704e21-9e89-4b5d-bb4e-7dd24d6b506a
# ╠═008c4d22-4c2f-4087-8be7-47599600de47
# ╟─ac8f76d2-9801-4f64-8add-278c8ada3135
# ╠═34576988-89ce-431b-9e1c-78b1a98e002e
# ╟─9625fb8d-5e0e-4e07-bc49-b4cb9ab259ef
# ╠═33ca0d51-92fe-4313-ba67-7792af17b11e
# ╠═e5d2cee7-c2d4-45a5-901c-206be7d2ff58
# ╠═ab87421c-e498-44ff-af23-13ea14a0c5d3
# ╠═d9df727c-a96f-4af4-ab5a-60715662adce
# ╠═82fffd1c-9ec0-459c-b033-d926641268dc
# ╟─e823348a-1ea3-45b7-b1ad-d51096f224dc
# ╠═73ad34ac-e4fc-4e35-9d96-25969b37a4d6
# ╟─19dc36a5-6d13-4d61-a9b4-863e7fbe60fe
# ╠═6059b96d-8ba0-475c-9929-8cc613355562
# ╠═64d397ac-ed6f-426c-9139-af21222723da
# ╠═ba09dd4f-b114-45bb-9b6a-41c2ed9dbec5
# ╠═52506f6c-069a-4883-8160-1cf415e80c36
# ╠═c0fb2bfc-7221-428d-b13d-ebb483bd0c1c
# ╠═1fe960aa-f12a-4d3e-994a-1a5d5c10f29e
# ╠═b4e246f8-72b9-4042-8160-a1eda2084e45
# ╟─11259db5-f754-47c7-9e92-600adf6896f8
# ╠═ff33e971-4871-4533-95ef-4fd9b14f32d1
# ╠═f29daf79-022a-4f7d-8616-6cac0e7ea9e8
# ╟─7f7672da-787f-4628-83aa-e9297f7da322
# ╠═8f4fc906-99a7-4e06-8620-83cf427da8b9
# ╠═b687e934-69d9-4608-b73f-0f3676c210ec
# ╠═7ebf3fb9-85c7-4ea0-a18a-a4617375082d
# ╠═1f6187a2-41a2-43e0-9cd3-322a069222aa
# ╠═d1e5bfe4-1ac9-4dbe-aafb-f81822a72d3a
# ╠═02ef6f4f-3259-447e-878a-b116a3d0eb4a
# ╟─f7483843-9cc9-4799-ad56-1e74bd0a08e7
# ╠═b985704f-8f07-43b6-baa9-6b3d4541de27
# ╠═31a1b2aa-7b5b-46b0-82ec-f75111b72a0d
# ╠═a0e4f6af-d4ab-4b04-9863-49de95eec947
# ╠═c20593d7-885f-4f36-984c-ecded0884ecc
# ╟─cef37fa4-df23-4267-bde4-c1fe6ddade23
# ╠═92d27828-c97b-41d0-85e9-0d3b1738dfed
# ╠═6055ba96-dd58-420b-854e-277b60a3c1e7
# ╠═ca176d00-70f9-42b4-a3f4-fd5d37c9f31a
# ╠═a8080f5e-6797-482e-84ca-7ddfcd40264d
# ╠═2a644f10-2683-4ea8-bd25-95343abb7e48
# ╠═bc477197-2ac5-4798-a62b-62f89791fb94
# ╠═90043461-3f6a-45e6-a661-57c8df3a2a93
