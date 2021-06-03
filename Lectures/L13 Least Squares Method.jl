### A Pluto.jl notebook ###
# v0.12.8

using Markdown
using InteractiveUtils

# ╔═╡ bb659710-2421-11eb-25d4-af6542e470eb
begin
	import Random
	using LinearAlgebra
end

# ╔═╡ fcc44b72-e162-4351-8601-f7402e2ed694
md"""
# Least Squares Method

Consider a system of linear equations with more equations than unknowns:

$$Ax=b, \quad m>n.$$

If the system has a solution $x$, then $Ax-b=0$, so $\| Ax-b\|=0$ for every vector norm.

If the system does not have solution, then it is natural to seek $x$ such that

$$
\|Ax-b \|_{1,2,\infty}\to \min$$

for the chosen vecor norm.
"""

# ╔═╡ a062b872-1eaa-11eb-005f-9d66fad5ee28
md"
If $\mathop{\mathrm{rank}} A=n$, then the __unique__ $x$ for which 

$$
\|Ax-b \|_{2}\to \min$$

is the solution of the system of  __normal equations__:

$$
A^T A x=A^T b. \tag{*}$$
"

# ╔═╡ a6fc5380-1eaa-11eb-11a9-3544f456255c
md"
_Proof._  Define

$$
Q(x)=\|Ax-b\|_2^2=(x^TA^T-b^T)(Ax-b)=x^TA^T A x -2x^T A^T b+b^Tb.$$

It holds

$$\begin{aligned}
Q(x+h)&=(x^T+h^T)A^TA(x+h)-2(x^T+h^T)A^Tb+b^Tb \\
&=Q(x) +2h^T(A^TAx-A^Tb)+h^TA^TAh\\ &= Q(x)+\|Ah\|_2^2 \\
&\geq Q(x),
\end{aligned}$$

so the minimum is indeed attained at $x$.

This solution is unique since $Q(x)=Q(y)$ implies $\|Ah\|_2=0$, so either $h=0$ or $\mathop{\mathrm{rang}} A<n$, which is a contradiction, and the theorem is proved.
"

# ╔═╡ b4ddcb00-1eaa-11eb-23c6-d15643cd207a
md"
__Geometrical interpretation.__ Vektori $Ax$ i $Ax-b$ are mutually orthogonal, 

$$
(Ax)^T\cdot (Ax - b)=x^T (A^TAx - A^Tb)=0.$$ 

Therefore, $Ax$ is the orthogonal projection of the vector $b$ onto the set $\{Ay:\ y \textrm{ arbitrary}\}$.

The solution $x$ is called  __least squares solution__ 
of the system $A x=b$. __Relative residual__

$$
q=\sqrt{\frac{Q(x)}{Q(0)}}=\frac{\|A x - b\|_2}{\|b\|_2 }$$

measures the quality of the solution (adaptation).
"

# ╔═╡ 5c65704d-666f-4f15-bc8f-7741457f9af0
md"""
## Example

Let us solve the system


$$\begin{aligned}
x+y&=0\\
y+z&=1\\
x+z&=0\\
-x+y+z&=1\\
-x-z&=0
\end{aligned}$$

in the sense of least squares.
"""

# ╔═╡ 365fc919-988d-4a1b-b42c-b8ab6931f860
A=[1//1 1 0;0 1 1;1 0 1;-1 1 1;-1 0 -1]

# ╔═╡ dea071d2-66b7-44f0-a75c-e0c67e574561
b=[0//1,1,0,1,0]

# ╔═╡ d6e05d5d-6878-4865-934a-6d8846f1d157
x=(A'*A)\(A'*b)

# ╔═╡ 49874317-1bb6-4882-ae83-4644094bf87e
# Relative residual	
sqrt(norm(A*x-b)/norm(b))

# ╔═╡ 1052a0b7-748d-46d5-a5fb-9d3b1ba2b65e
md"
If the system is overdetermined, the standard command `/` computes 
the least squares solution using QR factorization:
"

# ╔═╡ de2b7ecb-df68-4280-a23b-d11bff24aa78
float(A)\float(b)

# ╔═╡ d8458b22-f8b3-4fda-8351-3abe5af8dc46
float(x)

# ╔═╡ c2f32713-6b34-4c34-9180-d759039891c5
md"""
## Example
"""

# ╔═╡ d92d5bed-8689-481f-b56f-d46fa4f835c1
begin
	Random.seed!(123)
	A₁=rand(20,10)
	b₁=rand(20);
end

# ╔═╡ 74feeef4-02b4-454c-9d86-c29d775c89c0
x₁=A₁\b₁

# ╔═╡ 73d14810-4471-40ad-b1ca-f1d210ad2eb2
q₁=sqrt(norm(A₁*x₁-b₁)/norm(b₁))

# ╔═╡ f7f6d6b8-44cc-4337-a601-3c4c5ef3fb77
md"""
## Perturbation theory

__Sensitivity of the least squares problem__ is given by following bounds (see [Matrix Computations, Section 5](https://books.google.hr/books?id=X5YfsuCWpxMC&printsec=frontcover&hl=hr#v=onepage&q&f=false)).

__Condition number__ of a general matrix $A$ is:

$$
\kappa_2(A)=\sqrt{\kappa(A^TA)}=\|A\|_2 \|(A^TA)^{-1} A^T\|_2.$$

Let $x$ and $\hat x$ be the least squares solutions of the systems $Ax=b$ and 
$(A+\delta A)\hat x=b+\delta b$, respectively. The __residuals__ are defined by

$$
\begin{aligned}
r&=Ax-b\\
\hat r&=(A+\delta A)\hat x-(b+\delta b).
\end{aligned}$$

Let

$$
\epsilon=\max \bigg\{ \frac{\|\delta A\|_2}{\|A\|_2},\frac{\|\delta b\|_2}{\|b\|_2}\bigg\}$$

and

$$
q=\frac {\|r\|_2}{\|b\|_2}\equiv\sin\theta <1.$$

Then,

$$
\begin{aligned}
\frac{\|\hat x-x\|_2}{\|x\|_2}&\leq \epsilon \bigg[\frac{2\,\kappa_2(A)}{\cos \theta} +\tan\theta \,\kappa_2^2(A)\bigg]+O(\epsilon^2),\\
\frac{\|\hat r-r\|_2}{\|b\|_2}&\leq \epsilon\,[1+ 2\,\kappa_2(A)](m-n)+O(\epsilon^2).
\end{aligned}$$

We see that the residual itself is less sensitive than the position where it is attained.
"""

# ╔═╡ fbcad473-e073-4edf-a46f-f1d28a3b753d
cond(A₁)

# ╔═╡ 4dcccfd8-fdf6-4cdc-8314-c2a0d002baff
δA₁=1e-4*(rand(20,10).-0.5)

# ╔═╡ bb45e359-ca50-49cd-a23c-16552fac65ed
xp₁=(A₁+δA₁)\b₁

# ╔═╡ 4bd4192c-4d5e-4869-b598-9edd04128c7a
begin
	r₁=A₁*x₁-b₁
	rp₁=(A₁+δA₁)*xp₁-b₁
end

# ╔═╡ d7a8a461-17f6-4a9e-9343-f58314e6a9ee
norm(xp₁-x₁)/norm(x₁), norm(rp₁-r₁)/norm(b₁)

# ╔═╡ 314c1246-1772-415f-9f9a-38b221b6eb96
md"""
## Error analysis and accuracy

If $\mathop{\mathrm{rang}}A =n$, the matrix $A^TA$ is symmetric and positive definite, so the system 
(*) can be solved using Cholesky factorization.

The computed solution $\hat x$ satisfies

$$
(A^TA +E)\hat x=A^Tb,$$

where

$$
\|E\|_2\approx \varepsilon \| A^TA\|_2,$$

so the bound for the relative error is

$$
\frac{\|\hat x -x\|_2}{\|x\|_2}\approx \varepsilon \kappa_2(A^TA) =\varepsilon \kappa^2_2(A).$$


Therefore, the relative error of the solution obtained using normal equation depends upon the __square of the condition number__, so it is better to use QR factorization.
"""

# ╔═╡ Cell order:
# ╟─fcc44b72-e162-4351-8601-f7402e2ed694
# ╟─a062b872-1eaa-11eb-005f-9d66fad5ee28
# ╟─a6fc5380-1eaa-11eb-11a9-3544f456255c
# ╟─b4ddcb00-1eaa-11eb-23c6-d15643cd207a
# ╠═bb659710-2421-11eb-25d4-af6542e470eb
# ╟─5c65704d-666f-4f15-bc8f-7741457f9af0
# ╠═365fc919-988d-4a1b-b42c-b8ab6931f860
# ╠═dea071d2-66b7-44f0-a75c-e0c67e574561
# ╠═d6e05d5d-6878-4865-934a-6d8846f1d157
# ╠═49874317-1bb6-4882-ae83-4644094bf87e
# ╟─1052a0b7-748d-46d5-a5fb-9d3b1ba2b65e
# ╠═de2b7ecb-df68-4280-a23b-d11bff24aa78
# ╠═d8458b22-f8b3-4fda-8351-3abe5af8dc46
# ╟─c2f32713-6b34-4c34-9180-d759039891c5
# ╠═d92d5bed-8689-481f-b56f-d46fa4f835c1
# ╠═74feeef4-02b4-454c-9d86-c29d775c89c0
# ╠═73d14810-4471-40ad-b1ca-f1d210ad2eb2
# ╟─f7f6d6b8-44cc-4337-a601-3c4c5ef3fb77
# ╠═fbcad473-e073-4edf-a46f-f1d28a3b753d
# ╠═4dcccfd8-fdf6-4cdc-8314-c2a0d002baff
# ╠═bb45e359-ca50-49cd-a23c-16552fac65ed
# ╠═4bd4192c-4d5e-4869-b598-9edd04128c7a
# ╠═d7a8a461-17f6-4a9e-9343-f58314e6a9ee
# ╟─314c1246-1772-415f-9f9a-38b221b6eb96
