### A Pluto.jl notebook ###
# v0.12.4

using Markdown
using InteractiveUtils

# ╔═╡ c162d7c0-0efc-11eb-16ad-3d5ff641f6e0
using LinearAlgebra

# ╔═╡ 70942b12-0bd1-44e6-8e04-862d4948e8a0
md"""
# Positive Definite Matrices and Cholesky Factorization 


Matrix $A$ is _positive definite_ if it is symmetric and all its eigenvalues are non-negative.

Positive definite matrix can be factored (without pivoting) as

$$
A=L L^T$$

where $L$ is lower triangular matrix with non-negative diagonal elements. This factorization is called
__Cholesky factorization__ (see [Matrix Computations, section 4.2](https://books.google.hr/books?id=X5YfsuCWpxMC&printsec=frontcover&hl=hr#v=onepage&q&f=false)).

From

$$A=\begin{pmatrix}\alpha & a^T \cr a  & B \end{pmatrix}=
\begin{pmatrix} \beta & 0 \cr l & I \end{pmatrix}
\begin{pmatrix} \beta & l^T \cr 0 & C \end{pmatrix}
=LL^T=\begin{pmatrix} \beta^2 & \beta l^T \cr l\beta & ll^T+ C\end{pmatrix}$$

it follows

$$\beta=\sqrt{\alpha},\quad l=a\beta^{-1}=a\alpha^{-1/2},\quad C=B-ll^T=B-a\alpha^{-1}a^T.$$

Induction yields the following algorithm:

"""

# ╔═╡ a808c71a-97bb-4106-89fd-dcc0e902d2cb
function mychol(A₁::Matrix{T}) where T
    A=copy(A₁)
    n,m=size(A)
    for k=1:n
        A[k,k]=sqrt(A[k,k])
        for j=k+1:n
            A[k,j]=A[k,j]/A[k,k]
        end
        for j=k+1:n
            for i=k+1:n
                A[i,j]=A[i,j]-A[k,i]*A[k,j]
            end
        end
    end
    return triu(A)
end

# ╔═╡ ea4f0c1f-55e7-4423-a8df-b9de1889e511
begin
	A=rand(6,6)
	A=A*A'
end

# ╔═╡ dfea6e0a-a718-49ef-8e60-bd0ad6e204db
# Built-in function
C=cholesky(A)

# ╔═╡ 396ba935-a1fc-4978-a64e-990bb1d78f44
# Extract L from the structure
L=C.U

# ╔═╡ 1ff3e65c-8253-47d7-9caa-466a1359a434
# Residual 
L'*L-A

# ╔═╡ 286a74e5-6a3b-41ee-a227-cb5d88dd027e
# Our function
L₁=mychol(A)

# ╔═╡ 2e54b52f-2ba0-4de0-a7bf-73bf1c0d5b6f
# Residual
L₁'*L₁-A

# ╔═╡ 7decf367-ee06-401c-b742-b910e4b647a9


# ╔═╡ Cell order:
# ╟─70942b12-0bd1-44e6-8e04-862d4948e8a0
# ╠═a808c71a-97bb-4106-89fd-dcc0e902d2cb
# ╠═ea4f0c1f-55e7-4423-a8df-b9de1889e511
# ╠═c162d7c0-0efc-11eb-16ad-3d5ff641f6e0
# ╠═dfea6e0a-a718-49ef-8e60-bd0ad6e204db
# ╠═396ba935-a1fc-4978-a64e-990bb1d78f44
# ╠═1ff3e65c-8253-47d7-9caa-466a1359a434
# ╠═286a74e5-6a3b-41ee-a227-cb5d88dd027e
# ╠═2e54b52f-2ba0-4de0-a7bf-73bf1c0d5b6f
# ╠═7decf367-ee06-401c-b742-b910e4b647a9
