### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# ╔═╡ 7be18986-4e71-48cd-a897-586acfda6c5d
# When running on your computer, you can comment this cell
begin
	import Pkg
	Pkg.activate(mktempdir())
    Pkg.add([
        Pkg.PackageSpec(name="PlutoUI"),
		Pkg.PackageSpec(name="Plots")
    ])
end

# ╔═╡ 25c68fb0-2422-11eb-3214-6fbc147f322d
begin
	using PlutoUI, LinearAlgebra, Plots
	plotly()
end

# ╔═╡ bbc5eb02-6f79-4758-9bae-2591458dda5d
TableOfContents(title="📚 Table of Contents", aside=true)

# ╔═╡ ec0b9159-1c36-4112-8ce7-ec69c6709861
md"""
# Regression

__Regression__ is fitting a function $f$, which depends on $n$ parameters, through points 

$$(x_i,y_i),\quad i=1,2,\ldots, m,$$

where $m>n$, so that the __norm of deviations is minimal__:

$$
\| f(\mathbf{x})-\mathbf{y}\|_{1,2,\infty}\to \min,\quad f(\mathbf{x})=\begin{pmatrix}f(x_1) \\ f(x_2) \\ \vdots \\ f(x_m)\end{pmatrix},\quad
\mathbf{y}=\begin{pmatrix}y_1 \\ y_2 \\ \vdots \\ y_m\end{pmatrix}.$$

Regression in the __least squares sense__ is

$$
\| f(\mathbf{x})-\mathbf{y}\|_{2}\to \min.$$

When $f$ is a line,

$$
f(x)=kx+l,$$

we talk about __linear regression__. In this case we obtain system of linear equations

$$
k x_i + l=y_i, \quad i=1,2,\ldots,m.$$

If all points __are not on the same line__, the system is not solvable, so we are looking for the least squares solution.
"""

# ╔═╡ 108def07-d919-498c-82c2-7c43c7fd27f3
md"""
## Linear regression

Let us pass a line through the points 

$$(1,1), \ (2,3),\ (4,2), \ (6,4), \ (7,3),$$

and compute the relative residual. 
"""

# ╔═╡ 64ebb9c6-9c26-42bb-9141-45f5efab3d36
begin
	x=[1,2,4,6,7]
	y=[1,3,2,4,3]
	A=[x ones(Int,length(x))]
end

# ╔═╡ 5a57b1b6-eabb-46be-9986-1c48dcaf1785
# Coefficients of the regression line
(k,l)=A\y

# ╔═╡ 748be150-4eee-4c29-93fe-bdc1058cb3e8
begin
	# Plot the points and the regression line
	f(x)=k*x+l
	scatter(x,y,label="Points",legend=:bottomright)
	plot!(x->x,f,x[1],x[end],label="Regression line")
end

# ╔═╡ 1b894260-f180-4407-8c68-ad93ad40c55b
# Relative residual
sqrt(norm(A*[k;l]-y)/norm(y))

# ╔═╡ 986e6cb9-4a13-49e2-9891-29020262d002
md"""
## Quadratic regression

Through the given points, we can fit a quadratic polynomial $y=ax^2+bx+c$. If all points __do not lie on the same curve__, the system of linear equations

$$
ax_i^2+bx_i+c=y_i, \quad i=1,\ldots,m,$$

is not solvable, so we are looking for the least squares solution.

Let us fit quadratic polynomial through the points

$$
(1,0),\ (2,1), \ (4,4),\ (5,8), \ (6,14).$$
"""

# ╔═╡ a786bb3e-e442-4deb-a7ce-ed25fbf6bd75
begin
	x₁=[1,2,4,5,6]
	y₁=[0,1,4,8,14]
	A₁=[x₁.^2 x₁ ones(Int,length(x₁))]
end

# ╔═╡ a98c4043-4f99-4aa4-977a-2f5f777e762b
# Coefficients of the regression polynomial
(a,b,c)=A₁\y₁

# ╔═╡ 544441c1-f921-4767-859b-804f5fbd9508
# Plot the points and the parabola
g(x)=a*x^2+b*x+c

# ╔═╡ 93fdb6ef-09d2-47b7-b5af-6561529f520f
begin
	scatter(x₁,y₁,label="Points",legend=:topleft)
	plot!(x->x,g,x₁[1],x₁[end],label="Regression parabola")
end

# ╔═╡ 65f7a8eb-df16-4a01-91ff-0214ac092e61
# Relative residual
sqrt(norm(A₁*[a;b;c]-y₁)/norm(y₁))

# ╔═╡ b6166388-1d3b-4dce-969a-b7fce38426c1
md"""
## Growth of World population

The growth of World population is given in the following table, see [http://en.wikipedia.org/wiki/World_population](http://en.wikipedia.org/wiki/World_population): 


$$
\begin{array}{c|c|c|c|c|c|c|c|c|c}
\textrm{Year} & 1750 & 1800 & 1850 & 1900 & 1950 & 1999 & 2008 & 2010 & 2012 \\ \hline
\textrm{Population (milllions)} & 791 & 978 & 1262 & 1650 & 2521 & 5978 & 6707 & 6896 & 7052 
\end{array}$$




Let us approximate the population growth with an exponantial function, 
 

$$
P(t)=Ce^{kt},$$

and predict the population in the year 2050.

By taking logarithms, the system of equations 

$$
Ce^{kt_i}=P_i, \quad i=1,2,\ldots, 9,$$

is transformed to the system of linear equations

$$
k \,t_i + \ln C =\ln P_i.$$

__Not all points lie on the same curve__, so system is not solvable, and we are looking for the least squares solution.
"""

# ╔═╡ d86cd66b-5e87-4f06-b991-a9c734e5b9fb
begin
	nₚ=9
	t=[1750,1800,1850,1900,1950,1999,2008,2010,2012]
	P=[791,978,1262,1650,2521,5978,6707,6896,7052]
	Aₚ=[t ones(Int,length(t))]
	(kₚ,C)=Aₚ\log.(P)
end

# ╔═╡ f1b3366d-c0c2-499c-86b7-c34f86d7998a
# Values on the curve
P₁(t)=exp(C)*exp(kₚ*t)

# ╔═╡ 998e1154-e985-49a9-8cce-4b436a887e2f
begin
	# Plot the points and the regression curve
	scatter(t,P,label="Population",legend=:topleft)
	plot!(t->t,P₁,t[1],t[end],label="Regression curve")
end

# ╔═╡ f48de770-1eb0-11eb-0d18-8f9690cd89c6
md"
__Question.__ What caused the knick in the population?
"

# ╔═╡ 1ac0491b-6cc9-4a14-af99-c6c6defc3074
# Prediction for the year 2050.
P₁(2050)

# ╔═╡ 23e34ea5-8dc6-4444-bcd6-94229abad290
md"
The computed prediciton is smaller than the one on the web. If we restrict ourselves the the period since the year 1950,  we have
"

# ╔═╡ 52e78105-19f7-4830-bf59-4e336daa9272
begin
	A₂= [t[5:end] ones(5)]
	(k₂,C₂)=A₂\log.(P[5:end])
	P₂(t)=exp(C₂)*exp(k₂*t)
	scatter(t[5:end],P[5:end],label="Population",legend=:topleft)
	plot!(t->t,P₂,t[5],t[end],label="Regression curve")
end

# ╔═╡ e8911abc-74bd-4daa-9336-5b10270616b1
# Prediction for the year 2050.
P₂(2050)

# ╔═╡ Cell order:
# ╠═7be18986-4e71-48cd-a897-586acfda6c5d
# ╠═25c68fb0-2422-11eb-3214-6fbc147f322d
# ╠═bbc5eb02-6f79-4758-9bae-2591458dda5d
# ╟─ec0b9159-1c36-4112-8ce7-ec69c6709861
# ╟─108def07-d919-498c-82c2-7c43c7fd27f3
# ╠═64ebb9c6-9c26-42bb-9141-45f5efab3d36
# ╠═5a57b1b6-eabb-46be-9986-1c48dcaf1785
# ╠═748be150-4eee-4c29-93fe-bdc1058cb3e8
# ╠═1b894260-f180-4407-8c68-ad93ad40c55b
# ╟─986e6cb9-4a13-49e2-9891-29020262d002
# ╠═a786bb3e-e442-4deb-a7ce-ed25fbf6bd75
# ╠═a98c4043-4f99-4aa4-977a-2f5f777e762b
# ╠═544441c1-f921-4767-859b-804f5fbd9508
# ╠═93fdb6ef-09d2-47b7-b5af-6561529f520f
# ╠═65f7a8eb-df16-4a01-91ff-0214ac092e61
# ╟─b6166388-1d3b-4dce-969a-b7fce38426c1
# ╠═d86cd66b-5e87-4f06-b991-a9c734e5b9fb
# ╠═f1b3366d-c0c2-499c-86b7-c34f86d7998a
# ╠═998e1154-e985-49a9-8cce-4b436a887e2f
# ╟─f48de770-1eb0-11eb-0d18-8f9690cd89c6
# ╠═1ac0491b-6cc9-4a14-af99-c6c6defc3074
# ╟─23e34ea5-8dc6-4444-bcd6-94229abad290
# ╠═52e78105-19f7-4830-bf59-4e336daa9272
# ╠═e8911abc-74bd-4daa-9336-5b10270616b1
