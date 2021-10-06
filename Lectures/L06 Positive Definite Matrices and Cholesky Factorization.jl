### A Pluto.jl notebook ###
# v0.16.0

using Markdown
using InteractiveUtils

# ╔═╡ c162d7c0-0efc-11eb-16ad-3d5ff641f6e0
using LinearAlgebra

# ╔═╡ 70942b12-0bd1-44e6-8e04-862d4948e8a0
md"""
# Positive Definite Matrices and Cholesky Factorization 


Matrix $A$ is __positive definite__ if it is symmetric and all its eigenvalues are non-negative.

Positive definite matrix can be factored (without pivoting) as

$$
A=L L^T$$

where $L$ is lower triangular matrix with non-negative diagonal elements. This factorization is called __Cholesky factorization__ (see [Matrix Computations, section 4.2](https://books.google.hr/books?id=X5YfsuCWpxMC&printsec=frontcover&hl=hr#v=onepage&q&f=false)).

From

$$A=\begin{pmatrix}\alpha & a^T \cr a  & B \end{pmatrix}=
\begin{pmatrix} \beta & 0 \cr l & I \end{pmatrix}
\begin{pmatrix} \beta & l^T \cr 0 & C \end{pmatrix}
=\begin{pmatrix} \beta^2 & \beta l^T \cr l\beta & ll^T+ C\end{pmatrix}$$

it follows

$$\beta=\sqrt{\alpha},\quad l=a\beta^{-1}=a\alpha^{-1/2},\quad C=B-ll^T=B-a\alpha^{-1}a^T.$$

We see that we must have $\alpha>0$. Also, since the matrix $B$ is symmetric, so is $C$.

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

# ╔═╡ 0f3ec098-91c5-42a2-a212-162b2e96d949
A==A'

# ╔═╡ bae32bb1-e9de-470c-8f62-b383612b095b
eigvals(A)

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

# ╔═╡ 4483ea9f-345a-4f24-8b1c-49252d6ae9d7
L-L₁

# ╔═╡ 9a8b96f8-4e81-4572-89b7-d3a50c3627ec
md"
## Cholesky factorization with pivoting

Due to positive definitesess, the element with largest modulus lies on the main diagonal. Therefore, the complete pivoting is achieved by chosing the largest element on the diagonal of each submatrix.
"

# ╔═╡ b7696bc8-64b2-44c8-88f5-1d3af2f1d614
Cₚ=cholesky(A,Val(true))

# ╔═╡ 9e7fd697-4e1f-4b40-b90d-1d5696d76c89
Cₚ.P

# ╔═╡ 29dd9784-a9d3-486a-8a5a-47e9c2898197
Lₚ=Cₚ.L

# ╔═╡ bbe91c0a-6cb2-4b89-bafa-a5ecfa9298a4
Lₚ*Lₚ'-Cₚ.P'*A*Cₚ.P

# ╔═╡ 22e42b00-f6c7-4dcf-bec7-9407045dc93f
Lₚ*Lₚ'-A[Cₚ.p,Cₚ.p]

# ╔═╡ 89f58a12-4cac-40c2-8ef4-ee77d575d1a6
md"
## Block-Cholesky
"

# ╔═╡ 78704bf8-d55b-45ff-ac33-976fde87ca49
function mycholb(A₁::Matrix{T}) where T
    A=copy(A₁)
    n,m=size(A)
    for k=1:n
		C=cholesky(A[k,k])
        A[k,k]=C.U
        for j=k+1:n
            A[k,j]=C.L\A[k,j] # Solving triangular system
        end
        for j=k+1:n
            for i=k+1:j
                A[i,j]-=transpose(A[k,i])*A[k,j]
			end
        end
    end
    return triu(A)
end

# ╔═╡ 6c1a68b1-c396-43d5-844c-5e980d04d3f1
# Generate block matrix
begin
	k,l=32,16   # 32,16
	Ab=[rand(k,k) for i=1:l, j=1:l]
	Ab=Ab'*Ab
end

# ╔═╡ bbe2bceb-9cbb-416e-bd60-d95e3417790f
Lb=mycholb(Ab)

# ╔═╡ 291e0ce8-58dc-40e2-a0b9-12d67bf09713
# Residual
norm(Lb'*Lb-Ab)

# ╔═╡ 6455e2e5-1bd3-4acb-9742-995e4330277e
# Converting block matrix into a standard one
unblock(A) = mapreduce(identity, hcat, [mapreduce(identity, vcat, A[:,i]) for i = 1:size(A,2)])

# ╔═╡ 138dfe1b-ff8e-4b30-922b-2c0a4a2a9e0a
Ab₀=unblock(Ab);

# ╔═╡ 78b16d3f-2982-4732-bbc5-b6a611acc2c4
cholesky(Ab₀);

# ╔═╡ 13956a4e-0c20-4a27-a467-90b66c51f560
md"
Execution times of our block algorithm `mycholb()` and the built-in function `cholesky()` are similar.
"

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[LinearAlgebra]]
deps = ["Libdl"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
"""

# ╔═╡ Cell order:
# ╠═c162d7c0-0efc-11eb-16ad-3d5ff641f6e0
# ╟─70942b12-0bd1-44e6-8e04-862d4948e8a0
# ╠═a808c71a-97bb-4106-89fd-dcc0e902d2cb
# ╠═ea4f0c1f-55e7-4423-a8df-b9de1889e511
# ╠═0f3ec098-91c5-42a2-a212-162b2e96d949
# ╠═bae32bb1-e9de-470c-8f62-b383612b095b
# ╠═dfea6e0a-a718-49ef-8e60-bd0ad6e204db
# ╠═396ba935-a1fc-4978-a64e-990bb1d78f44
# ╠═1ff3e65c-8253-47d7-9caa-466a1359a434
# ╠═286a74e5-6a3b-41ee-a227-cb5d88dd027e
# ╠═2e54b52f-2ba0-4de0-a7bf-73bf1c0d5b6f
# ╠═4483ea9f-345a-4f24-8b1c-49252d6ae9d7
# ╟─9a8b96f8-4e81-4572-89b7-d3a50c3627ec
# ╠═b7696bc8-64b2-44c8-88f5-1d3af2f1d614
# ╠═9e7fd697-4e1f-4b40-b90d-1d5696d76c89
# ╠═29dd9784-a9d3-486a-8a5a-47e9c2898197
# ╠═bbe91c0a-6cb2-4b89-bafa-a5ecfa9298a4
# ╠═22e42b00-f6c7-4dcf-bec7-9407045dc93f
# ╟─89f58a12-4cac-40c2-8ef4-ee77d575d1a6
# ╠═78704bf8-d55b-45ff-ac33-976fde87ca49
# ╠═6c1a68b1-c396-43d5-844c-5e980d04d3f1
# ╠═bbe2bceb-9cbb-416e-bd60-d95e3417790f
# ╠═291e0ce8-58dc-40e2-a0b9-12d67bf09713
# ╠═6455e2e5-1bd3-4acb-9742-995e4330277e
# ╠═138dfe1b-ff8e-4b30-922b-2c0a4a2a9e0a
# ╠═78b16d3f-2982-4732-bbc5-b6a611acc2c4
# ╟─13956a4e-0c20-4a27-a467-90b66c51f560
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
