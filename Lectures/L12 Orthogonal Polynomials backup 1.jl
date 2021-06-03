### A Pluto.jl notebook ###
# v0.10.0

using Markdown

# ╔═╡ 2bc2db29-cb52-4c0e-bf2c-161e55983850
md"""
# Orthogonal Polynomials

---


Let 

$$
L(x_0,x_1,\ldots,x_n)
$$

be a (sub)space spanned by linearly independent vectors (or functions) $x_0,x_1,\ldots,x_n$.

In other words, $L$ is the set of all linear combinations of the given vectors.

Using __Gram-Schmidt orthogonalization procedure__ we can compute __orthogonal basis__ of this (sub)space,

$$
y_0,y_1,\ldots,y_n, 
$$

for which it holds

$$
(y_i,y_j)=0,\quad i\neq j. \tag{1}
$$

Let 

\begin{align}
y_0&=x_0\\
y_1&=x_1-\frac{(x_1,y_0)}{(y_0,y_0)}y_0\\
y_2&=x_2-\frac{(x_2,y_0)}{(y_0,y_0)}y_0-\frac{(x_2,y_1)}{(y_1,y_1)}y_1\\
& \vdots \\
y_n&=x_n-\sum_{j=0}^{n-1} \frac{(x_n,y_j)}{(y_j,y_j)}y_j.
\end{align}

Each $y_j$ is a linear combination of $x_0,x_1,\ldots,x_j$ so $y_j$ are linearly independent and

$$
L(x_0,x_1,\ldots,x_n)=L(y_0,y_1,\ldots,y_n).
$$

The equations (1)are verified directly.

__Weighted scalar product__  of functions $f$ and $g$ on interval $[a,b]$ with __weight__ $\omega(x)>0$ is

$$
(f,g)_\omega=\int_a^b f(x)g(x)\omega(x)\, dx.
$$

Non-zero functions $f$ and $g$ are __orthogonal__ if $(f,g)_\omega=0$.

__Orthogonal polynomials__ are obtained by orthogonalizing polynomials

$$
1,x,x^2,x^3,\ldots,x^n. \tag{2}
$$

Different choices of the wieght functions yield different systems of orthogonal polynomials.  
"""

# ╔═╡ 5d5d666c-eea0-4ded-bf69-c6e62e26ca0a
md"""
## Legendre polynomials 

These polynomials are obtained by orthogonalizing  system (2) with

$$
[a,b]=[-1,1], \quad \omega(x)=1.
$$

We use the package  `SymPy.jl` for symbolic computation. 
"""

# ╔═╡ f30cef79-1d34-41c6-aea8-59ab921eca16
using SymPy

# ╔═╡ cc114cfe-d006-45b1-a740-c0100091fd96
a=-1
b=1
n=8
P=Array{Any,1}(undef,n)
x=Sym("x")
P[1]=x^0
ω(x)=1
for k=2:n
    P[k]=x^(k-1)
    for j=1:k-1
        P[k]=P[k]-SymPy.integrate(x->x^(k-1)*P[j]*ω(x),a,b)*P[j]/
        SymPy.integrate(x->P[j]*P[j]*ω(x),a,b)
    end
end

# ╔═╡ 66182305-6b82-4744-a706-6776b36374b2
md"""
Julia starts indexing with 1, so all indexes are shifted, that is

$$
P_0(x)=P[1], \ P_1(x)=P[2], \ldots
$$
"""

# ╔═╡ 36eb7445-30ce-4277-ac13-ceeb968acb6c
P[1]

# ╔═╡ 962e15fd-49d2-41ed-8f05-56ab5c4ed014
P[4]

# ╔═╡ 71b7a33b-dcac-4c90-a85f-4d290ccf7bd2
P[6]

# ╔═╡ 030c6c1c-1b9c-4637-94be-6d8881d394ef
P[7]

# ╔═╡ 42eb8ada-1a04-4e20-a7b3-6f42bcf74f59
P[8]

# ╔═╡ 160a71b1-8a49-46eb-aaa9-eca617faee4d
md"""
Polynomials $P_n$ are, up to a multiplication by a constant, equal to __Legendre__ polyanomials,

$$
L_n(x)=\frac{1}{2^n n!}\frac{d^n}{dx^n}(x^2-1)^n, \quad n=0,1,2,3,\ldots
$$
"""

# ╔═╡ a02227da-499f-47b0-856c-ee87aaaef311
n=8
L=Array{Any,1}(missing,n)
L[1]=x^0
for k=1:n-1
    L[k+1]=expand(diff((x^2-1)^k/(2^k*factorial(k)),x,k))
end

# ╔═╡ e93db0d5-2715-4570-980d-64d535f5f998
L[1], P[1]

# ╔═╡ 43baaf85-622a-4ee0-8351-c41c5c06a656
L[2],P[2]

# ╔═╡ 69501389-653c-43a1-8ae6-ff0c267de178
L[4],P[4]

# ╔═╡ ed231f8e-7c46-4888-baac-ed027b78f1ea
L[7]

# ╔═╡ 36d88df7-5d6b-4548-a0ff-1b093141c99b
P[7]

# ╔═╡ 285540ab-1f8f-4975-a2d5-c48b803cad51
L[7]*16/231

# ╔═╡ 3d272428-89c3-4705-ac0e-f3e75436eca4
md"""
Besides orthogonality, the following properties hold:

* $L_n(x)$  has $n$ distinct zeros on the interval $[-1,1]$, 
* __three term recurrence formula__ holds: 

$$ 
L_{n+1}(x)=\frac{2n+1}{n+1}\,x\, L_n(x)-\frac{n}{n+1} L_{n-1}(x).
$$

Let us compute the polynomials numerically and plot them:
"""

# ╔═╡ 1eed107a-c382-4c4b-9bf6-d1298040245f
using Polynomials
using Plots

# ╔═╡ 34869b1a-3c62-4579-a60b-19ef65ab9269
n=40
L=Array{Any,1}(missing,n)
L[1]=Polynomial([1])
L[2]=Polynomial([0,1])
for i=3:n
    L[i]=(2*i-3)*L[2]*L[i-1]/(i-1)-(i-2)*L[i-2]/(i-1)
    # @show i, length(L[i])
end

# ╔═╡ 39f42034-812b-4a37-8261-4ba4c3a34463
L[7]

# ╔═╡ 3fdd941c-3c82-4c03-bfff-1d3a0cef6234
xx=range(-1,stop=1,length=300)
# Try different values of k from 1 to 40
k=20
yy=L[k].(xx)
plot(xx,yy)

# ╔═╡ 3104b5e4-6e78-4c9d-9ed4-b62607319460
md"""
## Chebyshev  polynomials

__Chebyshev polynomials__ $T_n(x)$ are obtained by orthogonalizing system (2) with

$$
[a,b]=[-1,1], \quad \omega(x)=\frac{1}{\sqrt{1-x^2}}.
$$

Chebyshev polynomials have the following properties:

* it holds

$$
T_n(x)=\cos(n\arccos x),\quad n=0,1,2,3,\ldots,
$$

* $T_n(x)$  has $n$ distinct zeros on the interval $[-1,1]$, 

$$ 
x_k=\cos \bigg(\frac{2k-1}{n}\frac{\pi}{2} \bigg), \quad k=1,\ldots,n,
$$

* __three term recurrence formula__ holds: 

\begin{align*}
T_0(x)&=1,\\
T_1(x)&=x, \\ 
T_{n+1}(x)&=2\,x\,T_n(x)-T_{n-1}(x),\quad n=1,2,3,\ldots.
\end{align*}
 
__Remark.__

The recurrence formula follows fron the trigonometric __addition formula__

$$
\cos(n+1)\varphi+\cos(n-1)\varphi=2\cos\varphi\cos n\varphi,
$$

by substituting

$$
\arccos x=\varphi.
$$

This is also used to prove orthogonality. 
"""

# ╔═╡ 29c029db-aa12-42d9-ae25-8db084c55374
# Symbolic
n=8
T=Array{Any,1}(missing,n)
T[1]=x^0
T[2]=x
for k=2:n-1
    T[k+1]=expand(2*x*T[k]-T[k-1])
end

# ╔═╡ e78bdea9-0617-4fbf-8641-9edfa70065ce
T[3]

# ╔═╡ 7ba0aea6-9ba8-450b-8978-fc1230cab4e7
T[7]

# ╔═╡ b088c3c5-8125-4ec0-91ff-aa7259a647e8
T[8]

# ╔═╡ 0309f52d-ce1d-4e58-a455-0f4b202f64a1
# Numeric
n=50
T=Array{Any,1}(missing,n)
T[1]=Polynomial([1])
T[2]=Polynomial([0,1])
for i=3:n
    T[i]=2*T[2]*T[i-1]-T[i-2]
    # @show i, length(T[i])
end

# ╔═╡ 76b829ef-7b11-462d-b7f1-a7c5882569c1
xx=range(-1,stop=1,length=300)
# Try different values of k from 1 to 50
k=27
yy=T[k].(xx)
plot(xx,yy)

# ╔═╡ e118607b-1f8a-493f-bd55-acd19f13d3ff
md"""
## Change of interval

Orthogonal system of functions $\Phi_i$ on the interval $[-1,1]$ using transformation

$$
\gamma :[a,b]\to [-1,1],\quad \gamma(x)=\frac{2x}{b-a}-\frac{a+b}{b-a}
$$

becomes orthogonal system on the interval $[a,b]$,

$$
\Psi_i(x)=\Phi_i(\gamma(x)).
$$
"""

# ╔═╡ 12d87ac2-7ae6-4dd5-8e5b-5772769cb1e9
a=1
b=4
xx=collect(range(a,stop=b,length=300))
γ=2*xx/(b-a).-(b+a)/(b-a)
# Probajte razne vrijednosti k od 1 do 50
k=17
yy=T[k].(γ)
plot(xx,yy)

# ╔═╡ be6dc196-84d8-4416-a06d-df2a68c538b1


# ╔═╡ Cell order:
# ╟─2bc2db29-cb52-4c0e-bf2c-161e55983850
# ╟─5d5d666c-eea0-4ded-bf69-c6e62e26ca0a
# ╠═f30cef79-1d34-41c6-aea8-59ab921eca16
# ╠═cc114cfe-d006-45b1-a740-c0100091fd96
# ╟─66182305-6b82-4744-a706-6776b36374b2
# ╠═36eb7445-30ce-4277-ac13-ceeb968acb6c
# ╠═962e15fd-49d2-41ed-8f05-56ab5c4ed014
# ╠═71b7a33b-dcac-4c90-a85f-4d290ccf7bd2
# ╠═030c6c1c-1b9c-4637-94be-6d8881d394ef
# ╠═42eb8ada-1a04-4e20-a7b3-6f42bcf74f59
# ╟─160a71b1-8a49-46eb-aaa9-eca617faee4d
# ╠═a02227da-499f-47b0-856c-ee87aaaef311
# ╠═e93db0d5-2715-4570-980d-64d535f5f998
# ╠═43baaf85-622a-4ee0-8351-c41c5c06a656
# ╠═69501389-653c-43a1-8ae6-ff0c267de178
# ╠═ed231f8e-7c46-4888-baac-ed027b78f1ea
# ╠═36d88df7-5d6b-4548-a0ff-1b093141c99b
# ╠═285540ab-1f8f-4975-a2d5-c48b803cad51
# ╟─3d272428-89c3-4705-ac0e-f3e75436eca4
# ╠═1eed107a-c382-4c4b-9bf6-d1298040245f
# ╠═34869b1a-3c62-4579-a60b-19ef65ab9269
# ╠═39f42034-812b-4a37-8261-4ba4c3a34463
# ╠═3fdd941c-3c82-4c03-bfff-1d3a0cef6234
# ╟─3104b5e4-6e78-4c9d-9ed4-b62607319460
# ╠═29c029db-aa12-42d9-ae25-8db084c55374
# ╠═e78bdea9-0617-4fbf-8641-9edfa70065ce
# ╠═7ba0aea6-9ba8-450b-8978-fc1230cab4e7
# ╠═b088c3c5-8125-4ec0-91ff-aa7259a647e8
# ╠═0309f52d-ce1d-4e58-a455-0f4b202f64a1
# ╠═76b829ef-7b11-462d-b7f1-a7c5882569c1
# ╟─e118607b-1f8a-493f-bd55-acd19f13d3ff
# ╠═12d87ac2-7ae6-4dd5-8e5b-5772769cb1e9
# ╠═be6dc196-84d8-4416-a06d-df2a68c538b1
