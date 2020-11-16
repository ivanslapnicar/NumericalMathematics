### A Pluto.jl notebook ###
# v0.12.10

using Markdown
using InteractiveUtils

# ╔═╡ 2d895d76-dcb9-4d43-b985-b7e9c97d0cd0
begin
	using Polynomials
	using Plots
	using SpecialMatrices
end

# ╔═╡ cf690c86-f1eb-418f-992f-e569937ff499
begin
	# Generate points
	using Random, LinearAlgebra
	Random.seed!(123)
	n=5
	x=sort(rand(n+1))
	y=rand(n+1)
end

# ╔═╡ 5ec7a38e-cecd-45c0-94a1-39d58e437add
md"""
# Natural Cubic Spline


Given is a function $f(x)$ on an interval $[a,b]$.

Choose $n+1$ point 

$$
a\equiv x_0<x_1<x_2<\cdots <x_n\equiv b$$ 

and compute the values 

$$
y_i=f(x_i), \quad i=0,1,\ldots,n.$$
 
On the interval $[x_{i-1},x_i]$ the function $f$ is approximated by cubic polynomial $C_i$.
Therefore, on the entire interval $[a,b]$ function $f$ is approximated by the function 

$$
C(x)=C_i(x), \quad x\in[x_{i-1},x_i],\quad i=1,\ldots,n.$$

We require that $C(x)$ 

* is __continous__,
* its __first derivative is continious__, and
* its __second derivative is continious__.

Therefore,

$$\begin{aligned}
C_i(x_{i-1})&=y_{i-1}, \quad &i=1,\ldots,n, \\
C_i(x_{i})&=y_{i} \quad &i=1,\ldots, n,\\
C'_i(x_i)&=C'_{i+1}(x_i), \quad &i=1,\ldots,n-1, \\
C''_i(x_i)&=C''_{i+1}(x_i), \quad &i=1,\ldots,n-1,
\end{aligned}$$

which is a system of $4n-2$ equations and $4n$ unknowns (each of $n$ cubic polynomials has 4 coefficients).

We have:

$$
C_i(x)=y_{i-1}-s_{i-1}\frac{h_i^2}{6}+b_i(x-x_{i-1})+\frac{s_{i-1}}{6h_i}(x_i-x)^3
+\frac{s_i}{6h_i}(x-x_{i-1})^3,$$

where

$$\begin{aligned}
b_i&=d_i-(s_i-s_{i-1})\frac{h_i}{6},\\
d_i&=\frac{y_i-y_{i-1}}{h_i},\\
h_i&=x_i-x_{i-1},
\end{aligned}$$

and numbers $s_i$, $i=0,1,\ldots,n$, satisfy the system of linear equations  

$$
s_{i-1}h_i+2s_i(h_i+h_{i+1})+s_{i+1}h_{i+1}=6(d_{i+1}-d_i),\quad i=1,\ldots,n-1.$$

If we fix $s_0$ and $s_n$, the system will have unique solution. 

We most often use __natural conditions__:

$$
s_0=0, \quad s_n=0,$$

hence the name __natural spline__.

In this case, $s_1,\ldots,s_{n-1}$ are solutions of the system

$$
{\small
\begin{pmatrix} 2(h_1+h_2) & h_2 & 0 & \cdots & 0 & 0 \\
h_2 & 2(h_2+h_3) & h_3 & \cdots & 0 & 0 \\
0 & h_3 & 2(h_3+h_4) & \cdots & 0 & 0 \\
\vdots & \vdots & \vdots & \ddots & \vdots & \vdots \\
0 & 0 & 0 & \cdots & 2(h_{n-2}+h_{n-1}) & h_{n-1} \\
0 & 0 & 0 & \cdots & h_{n-1}  & 2(h_{n-1}+h_{n})\\
\end{pmatrix}
\begin{pmatrix}
s_1\\ s_2 \\ s_3 \\ \vdots \\ s_{n-2} \\ s_{n-1}
\end{pmatrix}
= 
\begin{pmatrix}
6(d_2-d_1)\\
6(d_3-d_2)\\
6(d_4-d_3) \\
\vdots \\
6(d_{n-1}-d_{n-2}\\
6(d_n-d_{n-1})
\end{pmatrix}.}$$


The proof can be found in [Numerical Analysis, Section 6.4, p. 349](https://books.google.hr/books/about/Numerical_Analysis.html?id=x69Q226WR8kC&redir_esc=y).

The system matrix is  __tridiagonal__ and __positive definite__, so the system can be solved using Cholesky factorization (without pivoting) using $O(n)$ operations. 


The following  __error bounds__ hold:

$$\begin{aligned}
\max |f(x)-C(x)| &\leq \frac{5}{384} \max h_i^4 \\
\max |f'(x)-C'(x)| &\leq \frac{1}{24} \max h_i^3 \\
\max |f''(x)-C''(x)| &\leq \frac{3}{8} \max h_i^2. 
\end{aligned}$$

These bound also can be used separately on each of the subintervals.
"""

# ╔═╡ c12ee2d1-c5fa-44de-92e5-3c4e096e0187
md"""
## Interpolation of random points
"""

# ╔═╡ a990adad-3abe-4dc3-be83-02e3ce5c46d4
function Spline(x,y)
    h=x[2:end]-x[1:end-1]
    d=(y[2:end]-y[1:end-1])./h
    H=SymTridiagonal(2*(h[1:end-1]+h[2:end]),h[2:end-1])
    b₀=6*(d[2:end]-d[1:end-1])
    s=H\b₀
    s=[0;s;0]
    # Define polynomials
    b=d-(s[2:end]-s[1:end-1]).*h/6
	n=length(x)-1
    C=Array{Any}(undef,n)
    C=[xx -> 
        y[i]-s[i]*h[i]^2/6+b[i]*(xx-x[i])+s[i]*(x[i+1]-xx)^3/(6*h[i])+s[i+1]*(xx-x[i])^3/(6*h[i]) 
        for i=1:n]
    return C
end 

# ╔═╡ e0f49f72-b3c9-40e8-94ba-6646edb0966d
function plotspline(C,x,xx)
	# Points on spline
	ySpline=Array{Float64}(undef,length(xx))
	for i=1:length(xx)
		for k=1:length(C)
			if xx[i]<=x[k+1]
				ySpline[i]=C[k](xx[i])
	            break
	        end
	    end
	end
	return ySpline
end

# ╔═╡ 1e6ed7c3-a1bb-4d85-8c48-18410365aa47
C=Spline(x,y)

# ╔═╡ bca31410-9afb-4b45-bc46-d4a98d926301
begin
	# Plot
	lsize=200
	xx=range(x[1],x[end],length=lsize)
	scatter(x,y,label="Points")
	ySpline=plotspline(C,x,xx)
	plot!(xx,ySpline,label="Spline")
end

# ╔═╡ 4afd29e7-8e7c-44b6-8897-62dd7317b983
md"""
Compare spline and interpolation polynomial:
"""

# ╔═╡ e667e8e4-79ca-43cb-b0c1-2e673e13cfe3
begin
	A=Vandermonde(x)
	p=Polynomial(A\y)
	yPoly=p.(xx)
	scatter(x,y,label="Points")
	plot!(xx,[ySpline yPoly],label=["Spline" "Polynomial"])
end

# ╔═╡ 4e9a4998-2a8b-4476-8b6f-d086bb5ca6c3
md"""
## Interpolating functions

Let us compare spline and interpolation polynomial for functions

$$f(x)=\sin(x), \quad x\in[0,\pi]$$

and

$$f(x)=1-|x-1|,\quad  x\in[0,2].$$

"""

# ╔═╡ 20f11b0f-110e-465a-84d5-0f78e99b14ff
begin
	n₁=6; a=0; b=pi; f(x)=sin.(x)
	# n₁=10; a=0; b=2; f(x)=1 .-abs.(x .-1)
	
	x₁=collect(range(a,stop=b,length=n₁+1))
	y₁=f(x₁)
	xx₁=collect(range(a,stop=b,length=lsize))
	
	# Polynomial
	A₁=Vandermonde(x₁)
	p₁=Polynomial(A₁\y₁)
	yPoly₁=p₁.(xx₁)
	
	# Spline
	C₁=Spline(collect(x₁),y₁)
	ySpline₁=plotspline(C₁,x₁,xx₁)
	
	# Funkction 
	yFun=f(xx₁)
	
	# Plot
	scatter(x₁,y₁,label="Points")
	plot!(xx₁,[ySpline₁ yPoly₁ yFun],label=["Spline" "Polynomial" "Function"])
end

# ╔═╡ Cell order:
# ╟─5ec7a38e-cecd-45c0-94a1-39d58e437add
# ╟─c12ee2d1-c5fa-44de-92e5-3c4e096e0187
# ╠═2d895d76-dcb9-4d43-b985-b7e9c97d0cd0
# ╠═cf690c86-f1eb-418f-992f-e569937ff499
# ╠═a990adad-3abe-4dc3-be83-02e3ce5c46d4
# ╠═e0f49f72-b3c9-40e8-94ba-6646edb0966d
# ╠═1e6ed7c3-a1bb-4d85-8c48-18410365aa47
# ╠═bca31410-9afb-4b45-bc46-d4a98d926301
# ╟─4afd29e7-8e7c-44b6-8897-62dd7317b983
# ╠═e667e8e4-79ca-43cb-b0c1-2e673e13cfe3
# ╟─4e9a4998-2a8b-4476-8b6f-d086bb5ca6c3
# ╠═20f11b0f-110e-465a-84d5-0f78e99b14ff
