### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# ╔═╡ b3adb110-f081-485e-8df0-8059a59d77c3
begin
	import Pkg
	Pkg.activate(mktempdir())
	Pkg.add("PlutoUI")
end

# ╔═╡ b897b490-2422-11eb-29f2-c1b07eae9d6e
begin
	using PlutoUI, Random, LinearAlgebra
	Random.seed!(123)
end

# ╔═╡ 7d696568-5099-453a-bc62-b5cc9da014cd
TableOfContents(title="📚 Table of Contents", aside=true)

# ╔═╡ fb66b990-4511-476f-8e77-87aa9c265041
md"""
# QR Factorization


__QR factorization__ of a $m\times n$ matrix $A$, where $m\geq n$, is

$$
A=QR,$$

where $Q$ is an __ortonormal matrix__ of size $m\times m$, or

$$
Q^TQ=Q Q^T=I,$$

and $R$ is a $m\times n$ upper triangular matrix.

We call an orthonormal matrix also __orthogonal matrix__.

For example,

$$\begin{aligned}
\begin{bmatrix} a_{11} & a_{12} & a_{13} \\
a_{21} & a_{22} & a_{23} \\
a_{31} & a_{32} & a_{33} \\
a_{41} & a_{42} & a_{43} \\
a_{51} & a_{52} & a_{53}
\end{bmatrix}=
\begin{bmatrix}
q_{11} & q_{12} & q_{13} & q_{14} & q_{15} \\
q_{21} & q_{22} & q_{23} & q_{24} & q_{25} \\
q_{31} & q_{32} & q_{33} & q_{34} & q_{35} \\
q_{41} & q_{42} & q_{43} & q_{44} & q_{45} \\
q_{51} & q_{52} & q_{53} & q_{54} & q_{55}
\end{bmatrix}
\begin{bmatrix}
r_{11} & r_{12} & r_{13} \\
0 & r_{22} & r_{23} \\
0 & 0 & r_{33} \\
0 & 0 & 0 \\
0 & 0 & 0 
\end{bmatrix}. 
\end{aligned}$$

The above relation also defines an  __economical QR factorization__

$$\begin{aligned}
\begin{bmatrix} a_{11} & a_{12} & a_{13} \\
a_{21} & a_{22} & a_{23} \\
a_{31} & a_{32} & a_{33} \\
a_{41} & a_{42} & a_{43} \\
a_{51} & a_{52} & a_{53}
\end{bmatrix}=
\begin{bmatrix}
q_{11} & q_{12} & q_{13} \\
q_{21} & q_{22} & q_{23} \\
q_{31} & q_{32} & q_{33} \\
q_{41} & q_{42} & q_{43} \\
q_{51} & q_{52} & q_{53}
\end{bmatrix}
\begin{bmatrix}
r_{11} & r_{12} & r_{13} \\
0 & r_{22} & r_{23} \\
0 & 0 & r_{33}
\end{bmatrix}.
\end{aligned}$$
"""

# ╔═╡ af380e90-aa7c-411a-8b4d-ab16d5df6312
md"""
# Gram-Schmidt orthogonalization

Equating columns starting from the first one, gives:

$$
\begin{aligned}
t&=a_{:1}\\
r_{11}&=\|t\|_2 \\
q_{:1}&=t\frac{1}{r_{11}}\\
r_{12}&= q_{:1}^Ta_{:2} \\
t&=a_{:2}-q_{:1}r_{12} \\
r_{22}&=\|t\|_2 \\
q_{:2}&=t\frac{1}{r_{22}} \\
r_{13}&=q_{:1}^Ta_{:3} \\
r_{23}&=q_{:2}^Ta_{:3} \\
t&=a_{:3}-q_{:1}r_{13}-q_{:2}r_{23}\\
r_{33}&=\|t\|_2 \\
q_{:3}&=t\frac{1}{r_{33}}.
\end{aligned}$$

Induction yields __Gram-Schmidt orthogonalization procedure__.
"""

# ╔═╡ 24c9f0ba-6e18-4362-840b-1a3581032ff5
function GramSchmidtQR(A::Array)
    m,n=size(A)
    R=zeros(n,n)
    Q=Array{Float64}(undef,m,n)
    R[1,1]=norm(A[:,1])
    Q[:,1]=A[:,1]/R[1,1]
    for k=2:n
        for i=1:k-1
            R[i,k]=Q[:,i]⋅A[:,k]
        end
        t=A[:,k]-sum([R[i,k]*Q[:,i] for i=1:k-1])
        R[k,k]=norm(t)
        Q[:,k]=t/R[k,k]
    end
    return Q,R
end 

# ╔═╡ d98a860c-b31c-40df-bb9e-0993c9e69dd6
md"
## Example
"

# ╔═╡ b6ee2c00-78e7-4d22-a851-d22c4a4946bc
A=rand(8,5)

# ╔═╡ b8994eec-d651-4921-b1a3-184a029e83b7
Q,R=GramSchmidtQR(A)

# ╔═╡ 19a0f0b8-dbab-4609-9f7f-17339af6bf63
Q

# ╔═╡ 54f456e7-133d-4f96-901f-5bcda000178f
Q'*Q

# ╔═╡ 1d11a929-212b-42d0-ab52-c793ed8905d2
R

# ╔═╡ 9038e960-55c0-437a-a9fd-73f48860a36a
# Residual
A-Q*R

# ╔═╡ b9bf7235-dcb8-47bf-8c3f-1179ba9e08e0
md"""
Algorithm `GramSchmidtQR()` is numerically unstable, so it is better to use __modified Gram-Schmidt algorithm__ or __Householder reflectors__ or __Givens rotations__  (see [Matrix Computations, Section 5](https://books.google.hr/books?id=X5YfsuCWpxMC&printsec=frontcover&hl=hr#v=onepage&q&f=false)).
"""

# ╔═╡ 3fd60d80-25d7-45e2-9776-7250841e8f74
md"""
# Householder reflectors

__QR factorization of vector__ $x$ is

$$
H \begin{bmatrix} x_1 \\ x_2 \\ \vdots \\ x_m 
\end{bmatrix}  =r,$$

where

$$
H=I - \frac{2}{v^Tv}v v^T, \qquad  
v=\begin{bmatrix}
x_1\pm \|x\|_2 \\ x_2 \\ x_3 \\ \vdots \\ x_m
\end{bmatrix}.$$ 

__Householder reflector__ $H$ is __symmetric__ and __orthogonal__ matrix (__prove it!__). Depending on the choice of sign in the definition of the vector $v$, we have

$$
r=\begin{bmatrix} \mp \|x\| \\ 0 \\ \vdots \\ 0
\end{bmatrix}.$$

For the sake of numerical stability, the standard choice is

$$
v_1=x_1+\mathop{\mathrm{sign}} (x_1) \|x\|_2.$$

Matrix $H$ is __not explicitly computed__, but the product $Hx$ is computed using formula

$$
Hx=x-\frac{2(v^Tx)}{v^Tv}v=x-\frac{2 (v\cdot x)}{v\cdot v}v,$$

which requires $O(6m)$ operations.
"""

# ╔═╡ 4a371793-1d5c-413f-8e89-dce986c972ef
function HouseholderVector(x::Array)
    # Computes v
    v=copy(x)
    v[1]=x[1]+sign(x[1])*norm(x)
    v
end

# ╔═╡ 814c6dd6-faa1-4211-8f64-bf349624dc37
begin
	x=rand(8)
	v=HouseholderVector(x)
	β=(2/(v⋅v))*(v⋅x)
	x-β*v
end

# ╔═╡ e282bbd3-6505-4afc-b304-3c327fa5db62
x

# ╔═╡ eb8f508a-9018-4896-bd82-e707e01157bc
norm(x)

# ╔═╡ df0bb186-4b2c-429a-b279-5843f97cd3f0
md"""
QR factorization of a matrix is computed by recursively applying QR factorization of a vector to its columns:
"""

# ╔═╡ 3bb4d086-a79a-44b1-92c3-5e29be4c36de
function HouseholderQR(A₁::Array)
    # Computes Q and R
    A=copy(A₁)
    m,n=size(A)
    Q=Matrix{Float64}(I,m,m) # eye
    for k=1:n
        v=HouseholderVector(A[k:m,k])
        β=(2/(v⋅v))*v
        A[k:m,k:n]=A[k:m,k:n]-β*(v'*A[k:m,k:n])
        Q[k:m,:]=Q[k:m,:]-β*(v'*Q[k:m,:])
    end
    R=triu(A)
    Q',R
end
    

# ╔═╡ d7ddbafd-8923-45c1-a0f2-55aa1dc86185
A

# ╔═╡ cb302477-773e-4fd9-a9d6-a9df50e45b18
Qₕ,Rₕ=HouseholderQR(A)

# ╔═╡ 12069de5-8064-437b-b96c-222baa7b4c7d
Qₕ'*A

# ╔═╡ 3ceae442-3a07-4400-933c-a04bdc8ee742
Rₕ

# ╔═╡ ef4702ee-44cc-49d8-85bb-39249c431036
md"""
The function `HouseholderQR()` is  illustrative. Professional programs have following properties:

* computing with block matrices (usual size of a block is 32 or 64),
* the Householder vector is scaled as $\hat v=v/v_1$. Thus,  $\hat v_1=1$, while the rest of the elements of $\hat v$'s are stored in the strict lower triangle of $A$,
* if the matrix $Q$ is required, the accumulation is done backwards using stored elements $\hat v$ (this reduces operation count),
* there is an option of returning economical factorization,
* there is an option of using __pivoting__ - in each step, the columnn of the current submatrix with largest norm is brought to the pivoting position, so

$$AP=QR,\quad |R_{kk}|\geq |R_{k+1,k+1}|$$

which can be used to determine __numerical rank__ of the matrix $A$.
"""

# ╔═╡ 4587b182-80fd-4bd7-87af-17074ea7c9d9
# ?qr # See the description

# ╔═╡ f944cdbe-794c-405e-a215-687cbb12b22b
# Return the QR object
F=qr(A)

# ╔═╡ 64e052c2-b798-4fb8-80fd-c018df2c8625
F.Q'*A

# ╔═╡ ef2ebf49-e7b2-4300-9a5c-3a2e42251e74
F.Q*F.R

# ╔═╡ 6bf029d7-68e2-47a4-87b9-8f66154933d4
Fₚ=qr(A,Val(true))

# ╔═╡ 68d0e359-9fb9-4a9b-9f13-e8039597e948
# Pivoting vector
Fₚ.p

# ╔═╡ dcc30d6f-3170-4e6d-ad51-741bce871a42
# Permutation matrix
Fₚ.P

# ╔═╡ ad9107ee-2abd-4dd7-b821-2e0a4bae4374
# Residual using permutation matrix
Fₚ.Q*Fₚ.R-A*Fₚ.P

# ╔═╡ 17cebd37-73b8-4699-a119-4237ab85d084
# Residual using pivot vektor
Fₚ.Q*Fₚ.R-A[:,Fₚ.p]

# ╔═╡ 60166ed0-dfe1-474a-91df-64c433867489
md"""
## Speed

The number of floating point operations needed to compute the QR factorization of a $n\times n$ matrix is $O\big(\frac{4}{3}n^3\big)$ to compute $R$ and $O\big(\frac{4}{3}n^3\big)$ to compute $Q$. 
"""

# ╔═╡ 7a2fb72c-2c5c-4e1a-895a-fb0eb5608876
begin
	n=512
	A₁=rand(n,n);
end

# ╔═╡ c7f654da-3ad0-4781-a1b0-f7f31a19bf74
# The output is in Julia terminal.
@time qr(A₁);

# ╔═╡ c0eaba40-c919-49af-96f8-700173e3a72a
@time qr(A₁,Val(true));

# ╔═╡ f6b568ec-20e5-4575-a62c-70d470589ef1
@time HouseholderQR(A₁);

# ╔═╡ 7e40ad34-0497-410d-9669-20ad8a9ca74b
md"""
## Accuracy

Matrices $\hat Q$ and $\hat R$ computed with the Householder method satisfy 

$$
\begin{aligned}
\hat Q^T\hat Q& =I+E, \qquad \|E \|_2\approx \varepsilon,\\ 
\| A-\hat Q\hat R\|_2& \approx \varepsilon\|A\|_2.
\end{aligned}$$

Also, there exists an orthogonal matrix $Q$ for which 

$$
\| A- Q\hat R\|_2\approx \varepsilon\|A\|_2.$$
"""

# ╔═╡ Cell order:
# ╠═b3adb110-f081-485e-8df0-8059a59d77c3
# ╠═b897b490-2422-11eb-29f2-c1b07eae9d6e
# ╠═7d696568-5099-453a-bc62-b5cc9da014cd
# ╟─fb66b990-4511-476f-8e77-87aa9c265041
# ╟─af380e90-aa7c-411a-8b4d-ab16d5df6312
# ╠═24c9f0ba-6e18-4362-840b-1a3581032ff5
# ╠═d98a860c-b31c-40df-bb9e-0993c9e69dd6
# ╠═b6ee2c00-78e7-4d22-a851-d22c4a4946bc
# ╠═b8994eec-d651-4921-b1a3-184a029e83b7
# ╠═19a0f0b8-dbab-4609-9f7f-17339af6bf63
# ╠═54f456e7-133d-4f96-901f-5bcda000178f
# ╠═1d11a929-212b-42d0-ab52-c793ed8905d2
# ╠═9038e960-55c0-437a-a9fd-73f48860a36a
# ╟─b9bf7235-dcb8-47bf-8c3f-1179ba9e08e0
# ╟─3fd60d80-25d7-45e2-9776-7250841e8f74
# ╠═4a371793-1d5c-413f-8e89-dce986c972ef
# ╠═814c6dd6-faa1-4211-8f64-bf349624dc37
# ╠═e282bbd3-6505-4afc-b304-3c327fa5db62
# ╠═eb8f508a-9018-4896-bd82-e707e01157bc
# ╟─df0bb186-4b2c-429a-b279-5843f97cd3f0
# ╠═3bb4d086-a79a-44b1-92c3-5e29be4c36de
# ╠═d7ddbafd-8923-45c1-a0f2-55aa1dc86185
# ╠═cb302477-773e-4fd9-a9d6-a9df50e45b18
# ╠═12069de5-8064-437b-b96c-222baa7b4c7d
# ╠═3ceae442-3a07-4400-933c-a04bdc8ee742
# ╟─ef4702ee-44cc-49d8-85bb-39249c431036
# ╠═4587b182-80fd-4bd7-87af-17074ea7c9d9
# ╠═f944cdbe-794c-405e-a215-687cbb12b22b
# ╠═64e052c2-b798-4fb8-80fd-c018df2c8625
# ╠═ef2ebf49-e7b2-4300-9a5c-3a2e42251e74
# ╠═6bf029d7-68e2-47a4-87b9-8f66154933d4
# ╠═68d0e359-9fb9-4a9b-9f13-e8039597e948
# ╠═dcc30d6f-3170-4e6d-ad51-741bce871a42
# ╠═ad9107ee-2abd-4dd7-b821-2e0a4bae4374
# ╠═17cebd37-73b8-4699-a119-4237ab85d084
# ╟─60166ed0-dfe1-474a-91df-64c433867489
# ╠═7a2fb72c-2c5c-4e1a-895a-fb0eb5608876
# ╠═c7f654da-3ad0-4781-a1b0-f7f31a19bf74
# ╠═c0eaba40-c919-49af-96f8-700173e3a72a
# ╠═f6b568ec-20e5-4575-a62c-70d470589ef1
# ╟─7e40ad34-0497-410d-9669-20ad8a9ca74b
