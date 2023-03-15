### A Pluto.jl notebook ###
# v0.17.1

using Markdown
using InteractiveUtils

# ╔═╡ a81eb2fe-bfb7-45d4-b983-ac3d43bdbb53
using PlutoUI, Random, LinearAlgebra

# ╔═╡ 66047586-f90b-4c82-bee1-d836dcdd064f
TableOfContents(title="📚 Table of Contents", aside=true)

# ╔═╡ 48167400-0d83-11eb-1c7d-359a2574c8b1
md"
# Gaussian Elimination


# Basic idea

The system of linar equations $Ax=b$
is solved in three steps (__without pivoting__):

1.  $A=LU\qquad$ (LU factorization, $O(\frac{2}{3}n^3)$ operations)
2.  $Ly=b\qquad$ (solving lower triangular system, $n^2$ operations)
3.  $Ux=y\qquad$ (solving upper triangular system, $n^2$ operations).

With pivoting we have

1.  $PA=LU$
2.  $Ly=P^T b$
3.  $Ux=y$
"

# ╔═╡ 221d2474-de59-4042-918f-534305d8708f
md"""

## Examples

The following class of problems gives rise to two separate types of issues, one we have already discussed, one we have not. Below  $\epsilon$ is the value generated by the function `eps()`.

Consider the system of linear equations


$\begin{aligned}\frac{\epsilon}{10} x_1 + x_2  = 1 \\ x_1 + x_2 = 2\end{aligned}$

A good approximate answer is $x_1 = x_2 =1$. Use the augmented system approach:

$$\begin{aligned}
& \left(\begin{array}{cc|c} \displaystyle\displaystyle\frac{\epsilon}{10} & 1 & 1 \cr 1 & 1 & 2 \end{array}\right) \sim
\left(\begin{array}{cc|c} \displaystyle \displaystyle\frac{\epsilon}{10} & 1 & 1 \cr 0 & 1-\displaystyle\displaystyle\frac{10}{\epsilon} & 2-\displaystyle\displaystyle\frac{10}{\epsilon}\end{array}\right)
\approx \left(\begin{array}{cc|c}  \displaystyle\displaystyle\frac{\epsilon}{10} & 1 & 1 \cr 0 & -\displaystyle\displaystyle\frac{10}{\epsilon} & -\displaystyle\displaystyle\frac{10}{\epsilon}\end{array}\right).
\end{aligned}$$

The last transformation is rounding to the machine precision $\epsilon$. The very significant "1" and "2" are __rounded away__ in the last line!

Back solve to get $x_1 = 0$, $x_2 = 1$. If you put these values back in the original system, note that $x_1+x_2 = 1$, so this is "way off."
"""

# ╔═╡ 89291e44-4ea6-4e74-99bd-a30e8ee4d895
[eps()/10 1;0 1-10/eps()]\[1; 2-10/eps()]

# ╔═╡ ed629240-6b9a-4a85-b412-97a06565a9cf
md"""
Again there is a "fix", it is called __partial pivoting__. Put largest uneliminated entry in the column in pivot or diagonal position :

$$\begin{aligned}
&\left(\begin{array}{cc|c}     1 & 1 & 2 \cr \displaystyle\frac{\epsilon}{10} & 1 & 1 \end{array}\right) \sim
\left(\begin{array}{cc|c}     1 & 1 & 2 \cr 0                   & 1-\displaystyle\frac{\epsilon}{10}&1-\displaystyle\frac{\epsilon}{10}\end{array}\right)
\approx   \left(\begin{array}{cc|c}     1 & 1 & 2 \cr 0                   & 1                    &1                    \end{array}\right)
\end{aligned}$$
"""

# ╔═╡ b8b17cd3-c828-438d-a2c0-13b1965ed778
[1 1;0 1-eps()/10]\[2; 1-eps()/10]

# ╔═╡ eb5aab00-0d73-11eb-2f45-771b2e23a5e3
md"""
This is the correct solution to machine precision.

Sometimes changing the algorithm does _no good at all_ ! Consider the system
$\begin{aligned}
(1+2\epsilon)x_1 + (1+2\epsilon)x_2 &= 2 \cr
(1+\epsilon)x_1 + x_2 &=2
\end{aligned}$

Using the augmented matrix approach

$$\left(\begin{array}{cc|c}
(1+2\epsilon)&     (1+2\epsilon )&     2 \cr
(1+\epsilon)&     1   &2 \end{array}\right).$$

Mulitply the first row by $\alpha = (1+\epsilon)/(1+2\epsilon)= 1-\epsilon + O(\epsilon^2)$
and add to the second row:

$$\left(\begin{array}{cc|c}
(1+2\epsilon)&     (1+2\epsilon )&     2 \cr
0           &     -\epsilon & 2\epsilon \end{array}\right).$$

The solution is $x_1 = 4$ and $x_2 =-2$. This is correct to machine precision.

A small change in the right hand side yields

$\left(\begin{array}{cc|c}
(1+2\epsilon)&     (1+2\epsilon )&     2+4\epsilon \cr
(1+\epsilon)&     1   &2 +\epsilon \end{array}\right).$

The correct answer is $x_1=x_2=1$, but with rounding we get  $x_1 =0$, $x_2 =2$. __Explain!__
Every trick I know except increasing the precision, yields similar wrong answers.
"""

# ╔═╡ c471e93a-5c8a-433f-880a-9e26788fa601
[1+2*eps() 1+2*eps(); 1+eps() 1]\[2+4*eps(); 2+eps()]

# ╔═╡ 9983cc49-9398-4772-9d04-3f7bbf2b47a1
[BigFloat(1)+2*eps() 1+2*eps(); 1+eps() 1]\[BigFloat(2)+4*eps(); 2+eps()]

# ╔═╡ de3e2152-7ec2-4366-a6fa-dc1476d19480
md"""
__Reason.__   IEEE Arithmetic rounds this system to


$$\left(\begin{array}{cc|c}
(1+2\epsilon)&     (1+2\epsilon )&     2+4\epsilon \cr
(1+\epsilon)&     1   &2           \end{array}\right) $$

which has the solution $x_1=0$ and $x_2=2$. This problem is very close to the singular system

$$\left(\begin{array}{cc|c} 1 & 1 & 2 \cr 1 & 1 & 2\end{array}\right)$$

which has the parametric solutions

$$\mathbf{x} =\begin{pmatrix} x_1 \\ x_2\end{pmatrix} = \begin{pmatrix} 1\\1\end{pmatrix}+
\beta \begin{pmatrix}-1 \\1\end{pmatrix}, \quad \beta \in \mathbb{R}.$$

Notice that $\begin{pmatrix} x_1 \\ x_2\end{pmatrix}= \begin{pmatrix} 1\\1\end{pmatrix}$ i $\begin{pmatrix} x_1 \\ x_2\end{pmatrix}=\begin{pmatrix} 0\\2\end{pmatrix}$ are two of those solutions.

__Question.__ What is the geometric interpretation of this system?
"""

# ╔═╡ c9bc269a-a306-4a35-acd8-aad3de58f56a
md"""
# LU factorization

$$A=\begin{pmatrix}\alpha & a^T \cr b  & B \end{pmatrix}=
\begin{pmatrix} 1 & 0 \cr l & I \end{pmatrix}
\begin{pmatrix} \beta & u^T \cr 0 & C \end{pmatrix}
=LU=\begin{pmatrix} \beta & u^T \cr l\beta & lu^T+ C\end{pmatrix}$$

implies

$$\beta=\alpha,\quad u=a,\quad l=b\beta^{-1},\quad C=B-lu^T=B-b\beta^{-1}a^T.$$

Induction yields the following algorithm:
"""

# ╔═╡ 66fcc372-4f27-4092-9552-8eb6e863bd4a
function mylu(A₁::Array{T}) where T # Strang, p. 100
    A=copy(A₁)
    n,m=size(A)
    # This accepts numbers and block matrices
    U=map(Float64,[zero(A[1,1]) for i=1:n, j=1:n])
    L=map(Float64,[zero(A[1,1]) for i=1:n, j=1:n])
    for k=1:n
        L[k,k]=one(A[1,1])
        for i=k+1:n
            L[i,k]=A[i,k]/A[k,k]
            for j=k+1:n
                A[i,j]=A[i,j]-L[i,k]*A[k,j]
            end
        end
        for j=k:n
            U[k,j]=A[k,j]
        end
    end
    return L,U
end

# ╔═╡ b62492f0-0d7f-11eb-0e42-d377410cec70
Random.seed!(123);

# ╔═╡ 19567c60-0d82-11eb-3405-8d0312a34b5f
begin
	n=6
	A=rand(n,n)
	b=rand(n)
end

# ╔═╡ 1e51f7d0-0d82-11eb-238b-f1179d7f9a30
A

# ╔═╡ 251a1cf0-0d82-11eb-3747-cf84a824570f
L,U=mylu(A)

# ╔═╡ 2a725e60-0d82-11eb-18ce-cd017246fc46
L

# ╔═╡ 2d7fc570-0d82-11eb-06d0-b115bbbd897c
U

# ╔═╡ ddbe2df0-0d82-11eb-0a77-e9de7611ec9b
L*U-A

# ╔═╡ 946eb7cd-8a97-4aa8-880c-bfba8e6efae1
md"""
# Triangular systems
"""

# ╔═╡ 29d27f7f-d2e5-4d1f-a667-39fc266ffa17
begin
	function myU(U::Array{T},b₁::Array{T}) where T
	    b=copy(b₁)
	    n=length(b)
	    for i=n:-1:1
	       for j=n:-1:i+1
	            b[i]=b[i]-U[i,j]*b[j]
	       end
	        b[i]=b[i]/U[i,i]
	    end
	    b
	end

	function myL(L::Array{T},b₁::Array{T}) where T
	    b=copy(b₁)
	    n=length(b)
	    for i=1:n
	        for j=1:i-1
	            b[i]=b[i]-L[i,j]*b[j]
	        end
	        b[i]=b[i]/L[i,i]
	    end
	    b
	end
end

# ╔═╡ 6a7341a2-0d82-11eb-0831-05182f30ffe3
# Solve the system using the built-in function
x=A\b

# ╔═╡ f3732010-0d82-11eb-3484-c38c2d6a1f31
# Solve the system with our functions
y=myL(L,b)

# ╔═╡ f6da8e00-0d82-11eb-2179-b3cfdddf689a
x₁=myU(U,y)

# ╔═╡ 04064a10-0d83-11eb-2efe-7d7ca4ef35e7
# Compare the solutions
x-x₁

# ╔═╡ 8a3a81c0-0df9-11eb-087a-c9298b2fd265
@which lu(A)

# ╔═╡ 07fae3c0-0dfa-11eb-02a3-b7efa4490b1c
# lu

# ╔═╡ 4225c750-b668-4331-b6b8-0509635e69c6
md"""
# Speed

The function `mylu()` is slow. Among other reasons, it allocates unnecessarily three matrices and it does not work with block matrices.

The function can be reformulated such that $L$ and $U$ are stored in the array $A$, where the diagonal  of $L$ is not stored since all elements are 1 (see [Gilbert Strang, 'Introduction to Linear Algebra', p. 100](https://books.google.hr/books?id=M19gPgAACAAJ&dq=strang%20introduction&hl=hr&source=gbs_book_other_versions)):

"""

# ╔═╡ 17427400-0d83-11eb-14e2-f5d29e1650e4
function mylu₁(A₁::Array{T}) where T # Strang, p. 100
    A=copy(A₁)
    n,m=size(A)
    for k=1:n-1
        ρ=k+1:n
        A[ρ,k]=A[ρ,k]/A[k,k]
        A[ρ,ρ]=A[ρ,ρ]-A[ρ,k]*A[k,ρ]'
    end
    A
end

# ╔═╡ 0b44bb30-0d84-11eb-1d3a-dfc67cf3cf20
mylu₁(A)

# ╔═╡ 13a09fb2-0d84-11eb-15d6-e1d3b4ba658e
L

# ╔═╡ 1ccfd9c0-0d84-11eb-097b-2bfde3a9790f
U

# ╔═╡ 6f3f3257-8dd9-4d3c-b18e-cdbb37d52e2a
md"""
Compare execution times of LAPACK-based function `lu()` and our naïve function `mylu()` on larger matrix.

Execute several times for more accurate timings.
"""

# ╔═╡ 5da85850-0d84-11eb-091b-df4a89e4d052
A₁=rand(512,512);

# ╔═╡ 46868fc0-0d84-11eb-0bea-f9ee72af7795
lu(A₁)

# ╔═╡ 9ce0bb70-0d84-11eb-14ac-5335e9985dbf
mylu₁(A₁);

# ╔═╡ 09b70bcd-43ad-46b2-9664-2809351f9f70
md"""
## Block variant

`mylu()` and $\mathsf{mylu}_1()$ are tens of times slower than `lu()`.

Let us rewrite $\mathsf{mylu}_1()$ to work with blocks (still there is no pivoting!):
"""

# ╔═╡ e80e7d80-0d84-11eb-2423-23867085be67
function mylu₂(A₁::Array{T}) where T # Strang, page 100
    A=copy(A₁)
    n,m=size(A)
    for k=1:n-1
        for ρ=k+1:n
            A[ρ,k]=A[ρ,k]/A[k,k]
            for l=k+1:n
                A[ρ,l]=A[ρ,l]-A[ρ,k]*A[k,l]
            end
        end
    end
    A
end

# ╔═╡ 9d740956-6b11-4c0f-bca9-58fcfa852a62
md"""
First a small test, $k=2$, $l=4$:
"""

# ╔═╡ 090a3f10-0d85-11eb-0181-fdc5aa091df7
begin
	k,l=2,4   # 32,16
	Ab=[rand(k,k) for i=1:l, j=1:l]
end

# ╔═╡ 21e796b0-0dfb-11eb-2493-6fe0d3c70180
Ab[2,1]

# ╔═╡ 24c23190-0d85-11eb-382d-2536a040dfc8
A₀=mylu₂(Ab);

# ╔═╡ 3d486c20-0d85-11eb-2c81-a7d99c190907
begin
	# Provjera
	U₀=triu(A₀)
	L₀=tril(A₀)
	for i=1:maximum(size(L₀))
		L₀[i,i]=Matrix{Float64}(I,size(L₀[1,1])) # eye(L[1,1])
	end
end

# ╔═╡ 642b68b0-0d85-11eb-01d7-734195427bd9
L₀

# ╔═╡ 46ecd28e-0dfb-11eb-1e7a-2d55904f621e
L₀[1,1]

# ╔═╡ 9615e1c0-0d85-11eb-161d-39a7eff4046f
Residual=L₀*U₀-Ab

# ╔═╡ aa933e40-0d85-11eb-2141-d36ff3fc471b
# Converting block matrix into a standard one
unblock(A) = mapreduce(identity, hcat, [mapreduce(identity, vcat, A[:,i]) for i = 1:size(A,2)])

# ╔═╡ b307b3d0-0d85-11eb-25dc-83e2d3fbcb4f
norm(unblock(Residual))

# ╔═╡ 86f3ce48-73ce-4626-914f-478cf3ad1154
md"""
Try larger dimensions ($n=k\cdot l$).
"""

# ╔═╡ 83801210-b142-4d11-8eac-5eab72a181b3
md"""
We see that $\mathsf{mylu}_2()$ is nearly as fast as `lu()` (on one core), with the remark that $\mathsf{mylu}_2()$ uses no pivoting.
The function is still not optimal since it allocates to much memory.
"""

# ╔═╡ 9a2ef752-2231-4282-aa5c-275894c21de5
md"""
# Pivoting

Standard implementations always compute Gaussian elimination using _partial pivoting_ :

In each step, the rows are reordered (interchanged) such that the pivot element has largest absolute value in the active part of the given column. As a consequence

$$|L_{ij}| \leq 1,$$

which sufficiently reduces element growth (in practice).
"""

# ╔═╡ 535aa570-77c2-4962-97c2-9661884a21c2
A₂=[0.00003 1;2 3]

# ╔═╡ 4e739aee-0d86-11eb-056f-8589740ddc96
L₂,U₂=mylu(A₂)

# ╔═╡ 64384480-0d86-11eb-3ac7-7f602aeaa6d8
begin
	# With pivoting
	P=[0 1;1 0]
	mylu(P*A₂)
end

# ╔═╡ 416e1cd9-0fe8-4288-a1cb-b60f60139fa5
md"""
Now the built-in function. Use it precisely.
"""

# ╔═╡ e86a4730-0d86-11eb-266b-5b41924f61a8
F=lu(A₂)

# ╔═╡ fb2a1530-0d86-11eb-0d42-77246926cb7f
F.L*F.U == A₂[F.p, :]

# ╔═╡ 5b7d3a80-0dfe-11eb-05b2-2b3501e351fb
F.L*F.U==F.P*A₂

# ╔═╡ 68c72cf0-0dfe-11eb-37fb-d5ac2e1142e7
F.P

# ╔═╡ 4483e520-0d88-11eb-1eb4-25f4eaf1a88d
# We try the previous matrix A
A

# ╔═╡ 573618c0-0def-11eb-0101-e1b09c384512
F₄=lu(A)

# ╔═╡ 74b55dc0-0def-11eb-37e2-71ba59aa8295
F₄.L*F₄.U==A[F₄.p,:]

# ╔═╡ 90e4a320-0def-11eb-189e-7fc9c7b53942
# There are rounding errors
F₄.L*F₄.U-A[F₄.p,:]

# ╔═╡ 630a82aa-6998-4325-a0c9-d44f60c0df31
md"""
## Complete pivoting

The following function computes Gaussian elimination using __complete pivoting__ - in each step, rows and columns are interchanged such that the element with the largest absolute value in the current submatrix is brought to the pivot position. Theoretically, element growth is bounded, but the bound is $O(2^n)$ which is not useful.
"""

# ╔═╡ 18ad03b0-0d87-11eb-06f9-45ac1b7e3b04
function gecp(A₁::Array{T}) where T
    # Gaussian elimination with complete pivoting
    # On exit, either Pr*L*U*Pc'=A or Pr'*A*Pc=L*U
    A=copy(A₁)
    n,m=size(A)
    Pr=Matrix{Float64}(I,n,n)
    Pc=Matrix{Float64}(I,n,n)
    D=zeros(n)
    for i=1:n-1
        amax,indm=findmax(abs.(A[i:n,i:n]))
        imax=indm[1]+i-1
        jmax=indm[2]+i-1
        # Interchanging rows
        if (imax != i)
            temp = Pr[:,i]
            Pr[:,i] = Pr[:,imax]
            Pr[:,imax] = temp
            temp = A[i,:]
            A[i,:] = A[imax,:]
            A[imax,:] = temp
        end
        # Interchanging columns
        if (jmax != i)
            temp = Pc[:,i]
            Pc[:,i] = Pc[:,jmax]
            Pc[:,jmax] = temp
            temp = A[:,i]
            A[:,i] = A[:,jmax]
            A[:,jmax] = temp
        end
        # Elimination
        D[i]=A[i,i]
        A[i+1:n,i] = A[i+1:n,i]/D[i]
        A[i+1:n,i+1:n] = A[i+1:n,i+1:n] - A[i+1:n,i]*A[i,i+1:n]'
        A[i,i+1:n]=A[i,i+1:n]/D[i]
    end
    D[n]=A[n,n]
    L=I+tril(A,-1)
    U=I+triu(A,1)
    U=Diagonal(D)*U
    return L,U,Pr,Pc
end

# ╔═╡ 1f3a42b0-0d87-11eb-1fef-8f0a35eb3cce
Lₚ,Uₚ,Pr,Pc=gecp(A)

# ╔═╡ 3d328840-0d87-11eb-3c28-5dbb05be31f8
norm(Pr*Lₚ*Uₚ*Pc'-A)

# ╔═╡ 1ea5da30-0df0-11eb-2bb6-3bcf58c68adb
yₚ=myL(Lₚ,Pr'*b)

# ╔═╡ 284e14d0-0df0-11eb-2255-7b26982e1bbf
zₚ=myU(Uₚ,yₚ)

# ╔═╡ 355977a0-0df0-11eb-0e2b-0b5161d7979e
xₚ=Pc*zₚ

# ╔═╡ 3f504770-0df0-11eb-3049-ddac0626728f
# Residual
A*xₚ-b

# ╔═╡ 55d12cd0-0d87-11eb-10cc-edca8db298a1
md"""
# Accuracy

Consider the system $Ax=b$, where $A$ is nonsingular.

In order to apply concepts from the notebook
[L04 Backward Error and Stable Algorithms](L04%20Backward%20Error%20and%20Stable%20Algorithms.ipynb), we need to:

1. develop perturbation theory for the given problem, and
2. analyse errors in the algorithm (Gaussian elimination)

## Perturbation theory

Let

$$(A+\delta A)\hat x=(b+\delta b)$$

for some $\hat x=x+\delta x$.

We want to estimate

$$\frac{\| \hat x - x \|}{\| x\|} \equiv \frac{\| \delta x\|}{\| x\|}.$$

We introduce some notation (see [Matrix Computations, section 2.6.2](https://books.google.hr/books?id=X5YfsuCWpxMC&printsec=frontcover&hl=hr#v=onepage&q&f=false)):

$$\delta A=\varepsilon F, \quad \delta b=\varepsilon f, \qquad \hat x=x(\varepsilon)$$

which yields one-dimensional problem

$$(A+\varepsilon F)\,x(\varepsilon)=b+\varepsilon f$$

for some (unknown) matrix $F$ and vector $f$.

Differentiating with respect to $\varepsilon$ gives

$$Fx(\varepsilon)+(A+\varepsilon F)\, \dot x(\varepsilon)=f.$$

Setting $\varepsilon=0$ gives

$$F x+A\dot x(0)=f,$$

or

$$\dot x(0)=A^{-1}(f-Fx).$$

Taylor expansion arround $\varepsilon=0$ gives

$$x(\varepsilon)=x(0)+\varepsilon \dot x(0) +O(\varepsilon^2),$$

that is, by neglecting $O(\varepsilon^2)$ term,

$$\hat x-x=\varepsilon A^{-1}(f-Fx)=A^{-1} (\varepsilon f - \varepsilon F x) = A^{-1} (\delta b - \delta A x).$$

Properties of norm imply

$$\| \hat x-x\|\leq \| A^{-1} \| (\| \delta b \|  + \| \delta A \| \cdot \|  x\| ).$$

Finally, since $\| b\| \leq \| A\| \| x\|$, we have

$$\frac{\| \hat x-x\|}{\| x\|}\leq \| A\|  \cdot \| A^{-1} \| \bigg(\frac{\| \delta b \|}{\|b\|}  + \frac{\| \delta A \|}{ \|  A\|} \bigg). \tag{1}$$

The number

$$\kappa(A)\equiv \| A\|  \cdot \| A^{-1} \|$$

is __condition number__ (or __condition__) the matrix $A$, and it tells us
how much are relative changes in the input data (matrix $A$ and vector $b$)
relatively amplified in the solution.

Consider the following example:

"""

# ╔═╡ 06273c00-0d88-11eb-2259-230e34f04417
A₃= [0.234 0.458; 0.383 0.750]

# ╔═╡ 542b76a0-0d88-11eb-0672-c95813a3ccdc
b₃=[0.224;0.367]

# ╔═╡ 29a64d0e-0d88-11eb-0f2b-dfd116b214c4
mylu₁(A₃)

# ╔═╡ 3b73569e-0d88-11eb-271c-b983eb9cb3f5
F₃=lu(A₃)

# ╔═╡ f2d157d0-0df0-11eb-2ab9-4d55d1b5e307
x₃=A₃\b₃

# ╔═╡ fd0b7230-0df0-11eb-1443-67ea60bb2b7f
x₃[1]

# ╔═╡ 1e4c2c00-0df1-11eb-2ca2-3b4f08d92e9a
x₃[2]

# ╔═╡ 2bff8ea0-0df1-11eb-08a3-375007f3f276
begin
	δb₃=[0.00009; 0.000005]
	A₃\(b₃+δb₃)
end

# ╔═╡ 5179d372-0df1-11eb-0183-8bcc73149584
begin
	δA₃=[-0.001 0;0 0]
    x₄=(A₃+δA₃)\b₃
end

# ╔═╡ 69b8cbd0-0df1-11eb-3fcd-a9f0865efdce
cond(A₃), norm(δA₃)/norm(A₃), norm(x₄-x₃)/norm(x₃)

# ╔═╡ a47802e0-0df1-11eb-3f9f-2fe1ebc781fd
md"""
## Errors in Gaussian elimination

According to [Matrix Computations, section 3.3](https://books.google.hr/books?id=X5YfsuCWpxMC&printsec=frontcover&hl=hr#v=onepage&q&f=false), the COMPUTED factors
$\hat L$ and $\hat U$ satisfy

$$
\hat L\cdot \hat U = A+\delta A,$$

where (the inequality is interpreted elementwise, $\varepsilon$ is the machine precision)

$$
| \delta A|\leq 3(n-1) \varepsilon (|A|+|\hat L| \cdot |\hat U|) +O(\varepsilon^2).$$

Neglecting the $O(\varepsilon^2)$ term and taking norms yields

$$
\|\delta A \| \lesssim O(n)\varepsilon (\| A\| + \| \hat L\| \cdot \| \hat U\|),$$

so

$$
 \frac{\|\delta A \|}{\|A\|} \lesssim O(n)\varepsilon \bigg(1+\frac{\| \hat L\| \cdot \| \hat U\|}{\|A\|}\bigg).$$

If Gaussian eleimination is computed using row pivoting, then, most probably, the last quotient will be small ($\approx 1$). Further, the error in solving tirangular systems is not larger than this one, so inserting it into (1) it follows that the error in the computed solution satisfies

$$
\frac{\| \hat x-x\|}{\| x\|}\leq \kappa(A) O(n\varepsilon).$$

To conclude:

> _If the condition number is large, the solution may be inaccurate._

###  Vandermonde matrix
"""

# ╔═╡ f1c74560-0df1-11eb-19a7-c9ad6aed7410
begin
	nᵥ=10
	v=rand(nᵥ)
end

# ╔═╡ 026d0d00-0df2-11eb-26fd-cbc13048e56c
begin
	# Vandermonde matrices are notoriously ill-conditioned.
	V=Array{Float64}(undef,nᵥ,nᵥ)
	for i=1:nᵥ
	    V[:,i]=v.^(i-1)
	end
	V=V'
end

# ╔═╡ 2f58af40-0df2-11eb-28e4-5b1911f53b83
bᵥ=rand(nᵥ)

# ╔═╡ 3da4f680-0df2-11eb-23d4-f3fb28cdc8e7
xᵥ=V\bᵥ

# ╔═╡ 439224f0-0df2-11eb-2e57-539a3470de32
cond(V)

# ╔═╡ 50da42a0-0df2-11eb-2ded-c52a76acc155
begin
	Vbig=map(BigFloat,V)
	bbig=map(BigFloat,bᵥ)
	xbig=Vbig\bbig;
end

# ╔═╡ 56e4bd10-0df2-11eb-0152-556ef692a70e
map(Float64,norm(xbig-xᵥ)/norm(xbig))

# ╔═╡ 50ed08c1-391f-4baa-8cfa-54db04038fb1
md"""
## Artificial ill-conditioning
"""

# ╔═╡ a496c05e-0def-11eb-0ae1-83f3cdccf36e
Aᵤ=[1 1; 1 2]

# ╔═╡ aeee4670-0df2-11eb-0ef9-0bb353d8ebfe
bᵤ=[1;3]

# ╔═╡ aef0b770-0df2-11eb-3f66-8d09a5970a49
xᵤ=Aᵤ\bᵤ

# ╔═╡ aef1a1d0-0df2-11eb-0aeb-c5f670c48b32
xᵤ,cond(Aᵤ)

# ╔═╡ af144500-0df2-11eb-37e0-af9d5cd06c65
A₅=[1e-4 1e-4;1 2]

# ╔═╡ af16b600-0df2-11eb-2852-7f1e61314674
b₅=[1e-4;3]

# ╔═╡ af2b7680-0df2-11eb-2ff2-6b188cbd6b7f
x₅=A₅\b₅

# ╔═╡ af2dc072-0df2-11eb-2958-19fc8d01e81d
x₅,cond(A₅),xᵤ-x₅

# ╔═╡ bfc0ea15-556f-4109-b626-cb724ee14bfd
md"""
## Condition estimation

Computing the condition number according to the definition $\kappa(A)=\|A\| \cdot \|A^{-1}\|$ requires the inverse matrix, which requires $O(n^3)$ operacija. That is the same order of magnitude of operations needed to solve the entire system. However, when solving the system, the triangular factors $L$ and $U$ are already at our disposal, which can be used to approximate condition number using just $O(n^2)$ operation.
Details of this approach can be found in [Matrix Computations, section 3.5.4](https://books.google.hr/books/about/Matrix_Computations.html?id=X5YfsuCWpxMC&redir_esc=y).
LAPACK routine
[dtrcon.f](http://www.netlib.org/lapack/explore-html/d9/d84/dtrcon_8f_source.html) computes approximate reciprocal of the condition number of a triangular matrix.

Let us estimate the condition number of the Vandermonde matrix from the previous example.

"""

# ╔═╡ 8707357c-84c9-4d30-aa81-1172d7ac715e
#?LAPACK.trcon!

# ╔═╡ e38b8370-0df2-11eb-0ed5-ab750f73de17
begin
	Fᵥ=lu(V)
	cond(V,1), cond(Fᵥ.L,1), cond(Fᵥ.U,1)
end

# ╔═╡ 11add0f0-0df3-11eb-2f01-5b985574b265
1 ./LAPACK.trcon!('O','L','U',Fᵥ.L),1 ./LAPACK.trcon!('O','U','N',Fᵥ.U)

# ╔═╡ 24dc3f3e-0df3-11eb-04f6-a58c63e5ba58
md"""
## Residual


The computed solution $\hat x$ of the system $Ax=b$ is the exact solution of a nearby system (see [Afternotes on Numerical Analysis, str. 128](https://books.google.hr/books?id=w-2PWh01kWcC&printsec=frontcover&hl=hr#v=onepage&q&f=false)):


$$(A+\delta A)\,\hat x=b. \tag{1}$$

__Residual__ is defined as

$$r=b-A\hat x.$$

Then,

$$0=b-(A+\delta A)\,\hat x=r- \delta A\,\hat x.$$

Therefore,

$$\| r\| = \| \delta A\,\hat x \| \leq \| \delta A\| \cdot \|\hat x \|,$$

or

$$\frac{\|  \delta A\|}{\|A \|} \geq \frac{\|r\|}{\| A\| \cdot \|\hat x \|}.$$

Thus, if the   __relative rezidual__

$$\frac{r}{\| A\| \cdot \|\hat x \|}$$

has large norm, then  _the solution is not computed stably._

On the other side, if the relative residual is small in norm, then _the solution is computed stably_. Indeed, for (here we are using 2-norm)

$$\delta A=\frac{r\hat x^T}{\|\hat x\|^2_2}$$

(1) holds:

$$b-(A+\delta A)\hat x=(b-A\hat x)-\delta A \hat x = r-\frac{r\hat x^T \hat x}{\|\hat x\|^2_2}
= r-\frac{r \|\hat x^T \hat x\|_2}{\|\hat x\|^2}=r-r=0.$$

Also,

$$\frac{\|  \delta A\|_2}{\|A \|}  \leq  \frac{\|r\|_2\|\hat x \|_2}{\| A\| \cdot \|\hat x \|^2_2}=
\frac{\|r\|_2}{\| A\| \cdot \|\hat x \|_2}.$$

Let us compute residuals for the previous exmple of dimension $2$:

"""

# ╔═╡ 88f2801e-0df3-11eb-35ac-c32cb53aef8a
rᵤ=bᵤ-Aᵤ*xᵤ

# ╔═╡ e5108050-0df3-11eb-2d96-fddf9e91ef9e
norm(rᵤ)/(norm(Aᵤ)*norm(xᵤ))

# ╔═╡ f77005e0-0df3-11eb-0d2b-3b00dddef4f3
r₅=b₅-A₅*x₅

# ╔═╡ 0776ce62-0df4-11eb-1f95-3900b12d5087
norm(r₅)/(norm(A₅)*norm(x₅))

# ╔═╡ 1e40da00-0df4-11eb-3d74-03fb919b4781
md"
Residual for the Vandermonde system:
"

# ╔═╡ 2ec7cf00-0df4-11eb-03c6-03c37219650d
rᵥ=bᵥ-V*xᵥ

# ╔═╡ 3a45b400-0df4-11eb-3e3c-41fed6f7a499
norm(rᵥ)/(norm(V)*norm(xᵥ))

# ╔═╡ 40c1dc00-0df4-11eb-10a6-bf598057f7fb
md"
We conclude that the solution $x_v$ is computed stably, that is, with small backward error in the input data. This still does not mean that the relative error in the computed solution is small.
"

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[compat]
PlutoUI = "~0.7.9"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "8076680b162ada2a031f707ac7b4953e30667a37"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.2"

[[Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[LinearAlgebra]]
deps = ["Libdl"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[Parsers]]
deps = ["Dates"]
git-tree-sha1 = "438d35d2d95ae2c5e8780b330592b6de8494e779"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.0.3"

[[PlutoUI]]
deps = ["Base64", "Dates", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "Suppressor"]
git-tree-sha1 = "44e225d5837e2a2345e69a1d1e01ac2443ff9fcb"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.9"

[[Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[Random]]
deps = ["Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[Suppressor]]
git-tree-sha1 = "a819d77f31f83e5792a76081eee1ea6342ab8787"
uuid = "fd094767-a336-5f1f-9728-57cf17d0bbfb"
version = "0.2.0"

[[Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
"""

# ╔═╡ Cell order:
# ╠═a81eb2fe-bfb7-45d4-b983-ac3d43bdbb53
# ╠═66047586-f90b-4c82-bee1-d836dcdd064f
# ╟─48167400-0d83-11eb-1c7d-359a2574c8b1
# ╟─221d2474-de59-4042-918f-534305d8708f
# ╠═89291e44-4ea6-4e74-99bd-a30e8ee4d895
# ╟─ed629240-6b9a-4a85-b412-97a06565a9cf
# ╠═b8b17cd3-c828-438d-a2c0-13b1965ed778
# ╟─eb5aab00-0d73-11eb-2f45-771b2e23a5e3
# ╠═c471e93a-5c8a-433f-880a-9e26788fa601
# ╠═9983cc49-9398-4772-9d04-3f7bbf2b47a1
# ╟─de3e2152-7ec2-4366-a6fa-dc1476d19480
# ╟─c9bc269a-a306-4a35-acd8-aad3de58f56a
# ╠═66fcc372-4f27-4092-9552-8eb6e863bd4a
# ╠═b62492f0-0d7f-11eb-0e42-d377410cec70
# ╠═19567c60-0d82-11eb-3405-8d0312a34b5f
# ╠═1e51f7d0-0d82-11eb-238b-f1179d7f9a30
# ╠═251a1cf0-0d82-11eb-3747-cf84a824570f
# ╠═2a725e60-0d82-11eb-18ce-cd017246fc46
# ╠═2d7fc570-0d82-11eb-06d0-b115bbbd897c
# ╠═ddbe2df0-0d82-11eb-0a77-e9de7611ec9b
# ╟─946eb7cd-8a97-4aa8-880c-bfba8e6efae1
# ╠═29d27f7f-d2e5-4d1f-a667-39fc266ffa17
# ╠═6a7341a2-0d82-11eb-0831-05182f30ffe3
# ╠═f3732010-0d82-11eb-3484-c38c2d6a1f31
# ╠═f6da8e00-0d82-11eb-2179-b3cfdddf689a
# ╠═04064a10-0d83-11eb-2efe-7d7ca4ef35e7
# ╠═8a3a81c0-0df9-11eb-087a-c9298b2fd265
# ╠═07fae3c0-0dfa-11eb-02a3-b7efa4490b1c
# ╟─4225c750-b668-4331-b6b8-0509635e69c6
# ╠═17427400-0d83-11eb-14e2-f5d29e1650e4
# ╠═0b44bb30-0d84-11eb-1d3a-dfc67cf3cf20
# ╠═13a09fb2-0d84-11eb-15d6-e1d3b4ba658e
# ╠═1ccfd9c0-0d84-11eb-097b-2bfde3a9790f
# ╟─6f3f3257-8dd9-4d3c-b18e-cdbb37d52e2a
# ╠═5da85850-0d84-11eb-091b-df4a89e4d052
# ╠═46868fc0-0d84-11eb-0bea-f9ee72af7795
# ╠═9ce0bb70-0d84-11eb-14ac-5335e9985dbf
# ╟─09b70bcd-43ad-46b2-9664-2809351f9f70
# ╠═e80e7d80-0d84-11eb-2423-23867085be67
# ╟─9d740956-6b11-4c0f-bca9-58fcfa852a62
# ╠═090a3f10-0d85-11eb-0181-fdc5aa091df7
# ╠═21e796b0-0dfb-11eb-2493-6fe0d3c70180
# ╠═24c23190-0d85-11eb-382d-2536a040dfc8
# ╠═3d486c20-0d85-11eb-2c81-a7d99c190907
# ╠═642b68b0-0d85-11eb-01d7-734195427bd9
# ╠═46ecd28e-0dfb-11eb-1e7a-2d55904f621e
# ╠═9615e1c0-0d85-11eb-161d-39a7eff4046f
# ╠═aa933e40-0d85-11eb-2141-d36ff3fc471b
# ╠═b307b3d0-0d85-11eb-25dc-83e2d3fbcb4f
# ╟─86f3ce48-73ce-4626-914f-478cf3ad1154
# ╟─83801210-b142-4d11-8eac-5eab72a181b3
# ╟─9a2ef752-2231-4282-aa5c-275894c21de5
# ╠═535aa570-77c2-4962-97c2-9661884a21c2
# ╠═4e739aee-0d86-11eb-056f-8589740ddc96
# ╠═64384480-0d86-11eb-3ac7-7f602aeaa6d8
# ╟─416e1cd9-0fe8-4288-a1cb-b60f60139fa5
# ╠═e86a4730-0d86-11eb-266b-5b41924f61a8
# ╠═fb2a1530-0d86-11eb-0d42-77246926cb7f
# ╠═5b7d3a80-0dfe-11eb-05b2-2b3501e351fb
# ╠═68c72cf0-0dfe-11eb-37fb-d5ac2e1142e7
# ╠═4483e520-0d88-11eb-1eb4-25f4eaf1a88d
# ╠═573618c0-0def-11eb-0101-e1b09c384512
# ╠═74b55dc0-0def-11eb-37e2-71ba59aa8295
# ╠═90e4a320-0def-11eb-189e-7fc9c7b53942
# ╟─630a82aa-6998-4325-a0c9-d44f60c0df31
# ╠═18ad03b0-0d87-11eb-06f9-45ac1b7e3b04
# ╠═1f3a42b0-0d87-11eb-1fef-8f0a35eb3cce
# ╠═3d328840-0d87-11eb-3c28-5dbb05be31f8
# ╠═1ea5da30-0df0-11eb-2bb6-3bcf58c68adb
# ╠═284e14d0-0df0-11eb-2255-7b26982e1bbf
# ╠═355977a0-0df0-11eb-0e2b-0b5161d7979e
# ╠═3f504770-0df0-11eb-3049-ddac0626728f
# ╟─55d12cd0-0d87-11eb-10cc-edca8db298a1
# ╠═06273c00-0d88-11eb-2259-230e34f04417
# ╠═542b76a0-0d88-11eb-0672-c95813a3ccdc
# ╠═29a64d0e-0d88-11eb-0f2b-dfd116b214c4
# ╠═3b73569e-0d88-11eb-271c-b983eb9cb3f5
# ╠═f2d157d0-0df0-11eb-2ab9-4d55d1b5e307
# ╠═fd0b7230-0df0-11eb-1443-67ea60bb2b7f
# ╠═1e4c2c00-0df1-11eb-2ca2-3b4f08d92e9a
# ╠═2bff8ea0-0df1-11eb-08a3-375007f3f276
# ╠═5179d372-0df1-11eb-0183-8bcc73149584
# ╠═69b8cbd0-0df1-11eb-3fcd-a9f0865efdce
# ╟─a47802e0-0df1-11eb-3f9f-2fe1ebc781fd
# ╠═f1c74560-0df1-11eb-19a7-c9ad6aed7410
# ╠═026d0d00-0df2-11eb-26fd-cbc13048e56c
# ╠═2f58af40-0df2-11eb-28e4-5b1911f53b83
# ╠═3da4f680-0df2-11eb-23d4-f3fb28cdc8e7
# ╠═439224f0-0df2-11eb-2e57-539a3470de32
# ╠═50da42a0-0df2-11eb-2ded-c52a76acc155
# ╠═56e4bd10-0df2-11eb-0152-556ef692a70e
# ╟─50ed08c1-391f-4baa-8cfa-54db04038fb1
# ╠═a496c05e-0def-11eb-0ae1-83f3cdccf36e
# ╠═aeee4670-0df2-11eb-0ef9-0bb353d8ebfe
# ╠═aef0b770-0df2-11eb-3f66-8d09a5970a49
# ╠═aef1a1d0-0df2-11eb-0aeb-c5f670c48b32
# ╠═af144500-0df2-11eb-37e0-af9d5cd06c65
# ╠═af16b600-0df2-11eb-2852-7f1e61314674
# ╠═af2b7680-0df2-11eb-2ff2-6b188cbd6b7f
# ╠═af2dc072-0df2-11eb-2958-19fc8d01e81d
# ╟─bfc0ea15-556f-4109-b626-cb724ee14bfd
# ╠═8707357c-84c9-4d30-aa81-1172d7ac715e
# ╠═e38b8370-0df2-11eb-0ed5-ab750f73de17
# ╠═11add0f0-0df3-11eb-2f01-5b985574b265
# ╟─24dc3f3e-0df3-11eb-04f6-a58c63e5ba58
# ╠═88f2801e-0df3-11eb-35ac-c32cb53aef8a
# ╠═e5108050-0df3-11eb-2d96-fddf9e91ef9e
# ╠═f77005e0-0df3-11eb-0d2b-3b00dddef4f3
# ╠═0776ce62-0df4-11eb-1f95-3900b12d5087
# ╟─1e40da00-0df4-11eb-3d74-03fb919b4781
# ╠═2ec7cf00-0df4-11eb-03c6-03c37219650d
# ╠═3a45b400-0df4-11eb-3e3c-41fed6f7a499
# ╟─40c1dc00-0df4-11eb-10a6-bf598057f7fb
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
