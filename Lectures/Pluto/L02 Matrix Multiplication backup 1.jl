### A Pluto.jl notebook ###
# v0.10.0

using Markdown

# ╔═╡ 1db27bd5-d11b-4cee-8f5e-53427cf82786
md"""
# Matrix Multiplication



We can multiply matrices in  __three different ways__: 

\begin{align*}
\begin{bmatrix}    
1& 2& 3\\
  4 &5 &6\\
  7 &8& 9
\end{bmatrix}
\begin{bmatrix}
   1&  2&  0\\
   4&  3&  2\\
   1& -1&  1
 \end{bmatrix}&=
 \begin{bmatrix}
   (1\cdot 1 + 2\cdot 4+3\cdot 1) & 1\cdot 2+2\cdot 3+3\cdot(-1))&  
    (1\cdot 0 + 2\cdot 2 + 3\cdot 1)\\
   (4\cdot 1+5\cdot 4+6\cdot 1) & (4\cdot 2+5\cdot 3+6\cdot (-1)) &  
   (4\cdot 0+5\cdot 2+6\cdot 1)\\
   (7\cdot 1+8\cdot 5+9\cdot 1) & (7\cdot 2+8\cdot 3+9\cdot (-1)) &  
   (7\cdot 0+8\cdot 2+9\cdot 1)
 \end{bmatrix} \\ \\
\begin{bmatrix}    
  1& 2& 3\\
  4 &5 &6\\
  7 &8& 9
\end{bmatrix}
\begin{bmatrix}
   1&  2&  0\\
   4&  3&  2\\
   1& -1&  1
 \end{bmatrix}& =
 \begin{bmatrix} 1\\ 4\\ 7 \end{bmatrix}
 \begin{bmatrix} 1&2&0 \end{bmatrix} + 
 \begin{bmatrix} 2\\ 5\\ 8 \end{bmatrix}
 \begin{bmatrix} 4&3&2 \end{bmatrix} + 
 \begin{bmatrix} 3\\ 6\\ 9 \end{bmatrix}
 \begin{bmatrix} 1&-1&1 \end{bmatrix}\\ \\
\begin{bmatrix}    
  1& 2& 3\\
  4 &5 &6\\
  7 &8& 9
\end{bmatrix}
\begin{bmatrix}
   1&  2&  0\\
   4&  3&  2\\
   1& -1&  1
 \end{bmatrix}&=
 \begin{bmatrix}
   1 \begin{bmatrix} 1\\ 4\\ 7\end{bmatrix}+
   4 \begin{bmatrix} 2\\ 5\\ 8\end{bmatrix}+
   1 \begin{bmatrix} 3\\ 6\\ 9\end{bmatrix} &
   2 \begin{bmatrix} 1\\ 4\\ 7\end{bmatrix}+
   3 \begin{bmatrix} 2\\ 5\\ 8\end{bmatrix}+
   (-1) \begin{bmatrix} 3\\ 6\\ 9\end{bmatrix} &
   0 \begin{bmatrix} 1\\ 4\\ 7\end{bmatrix}+
   2 \begin{bmatrix} 2\\ 5\\ 8\end{bmatrix}+
   1 \begin{bmatrix} 3\\ 6\\ 9\end{bmatrix}
 \end{bmatrix}
\end{align*}

The formulas differ in the way they approach memory, and thus in speed, and in rounding errors.
"""

# ╔═╡ d4ca4508-4d02-41f3-b343-e7d9e3d9349c
n=5
A=ones(5,5)
B=ones(5,5)
C=zeros(n,n)
for i=1:n
    for j=1:n
        for k=1:n
            C[i,j]=C[i,j]+A[i,k]*A[k,j]
        end
        @show i,j,C[i,j]
    end
end
C

# ╔═╡ a05d46aa-b091-4551-9919-c97b404cf4fa
A

# ╔═╡ c60ff0fb-266e-43f2-926a-64913ad8cb3e
B

# ╔═╡ c56ebccc-e0c3-4926-a7f2-6363bd0cc563
# Inner loop: _dot
using LinearAlgebra
for i=1:n
    for j=1:n
        @show C[i,j]=(A[i,:]⋅B[:,j])[]
    end
end

# ╔═╡ 242cd1ff-49f5-476f-a9c1-12fc4bccea04
C=zeros(n,n)
for k=1:n
    for j=1:n
        for i=1:n
            C[i,j]=C[i,j]+A[i,k]*A[k,j]
            # @show i,j,C[i,j]
        end
    end
    @show C
end

# ╔═╡ 2ddf86e1-0751-4c0e-aee5-c3b1d4bed2fa
# Inner loop: _syrk
C=zeros(n,n)
for i=1:n
    C=C+A[:,i]*A[i,:]'
    println(C)
end

# ╔═╡ da5dab82-8aa2-4bf3-9273-69bb3675508f
C=zeros(n,n)
for j=1:n
    for k=1:n
        for i=1:n
            C[i,j]=C[i,j]+A[i,k]*A[k,j]
            @show i,j,C[i,j]
        end
    end
end

# ╔═╡ 7d8c987a-3a4e-48d9-8020-2170090f9934
# Inner loop: _axpy
C=zeros(n,n)
for i=1:n
    for k=1:n
        C[:,i]=C[:,i]+A[:,k]*B[k,i]
    end
    @show C
end

# ╔═╡ 8cef7e23-c539-4508-9176-37132955257d
md"""
`_axpy` is $y=\alpha x+y$
"""

# ╔═╡ 05adbae2-88a2-49ed-9ca6-129757f94266
md"""
## Basic Linear Algebra Subroutines - [BLAS](http://www.netlib.org/blas/)

For example, [`ddot.f`](http://www.netlib.org/lapack/explore-html/d5/df6/ddot_8f_source.html) - notice _loop unrolling_ :
"""

# ╔═╡ 6a085347-8f1d-4dc4-8742-774c4e279cfc
md"""
$x\cdot y=\sum_i x_i\cdot y_i=x_1y_1+x_2y_2+\cdots +x_ny_n$
"""

# ╔═╡ b3f1858f-b1ec-4f60-8343-b4b12d2c60c7
md"""
```
DOUBLE PRECISION FUNCTION DDOT(N,DX,INCX,DY,INCY)
*     .. Scalar Arguments ..
      INTEGER INCX,INCY,N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION DX(*),DY(*)
*     ..
*
*  Purpose
*  =======
*
*     forms the dot product of two vectors.
*     uses unrolled loops for increments equal to one.
*     jack dongarra, linpack, 3/11/78.
*     modified 12/3/93, array(1) declarations changed to array(*)
*
*
*     .. Local Scalars ..
      DOUBLE PRECISION DTEMP
      INTEGER I,IX,IY,M,MP1
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC MOD
*     ..
      DDOT = 0.0d0
      DTEMP = 0.0d0
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) GO TO 20
*
*        code for unequal increments or equal increments
*          not equal to 1
*
      IX = 1
      IY = 1
      IF (INCX.LT.0) IX = (-N+1)*INCX + 1
      IF (INCY.LT.0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
          DTEMP = DTEMP + DX(IX)*DY(IY)
          IX = IX + INCX
          IY = IY + INCY
   10 CONTINUE
      DDOT = DTEMP
      RETURN
*
*        code for both increments equal to 1
*
*
*        clean-up loop
*
   20 M = MOD(N,5)
      IF (M.EQ.0) GO TO 40
      DO 30 I = 1,M
          DTEMP = DTEMP + DX(I)*DY(I)
   30 CONTINUE
      IF (N.LT.5) GO TO 60
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
          DTEMP = DTEMP + DX(I)*DY(I) + DX(I+1)*DY(I+1) +
     +            DX(I+2)*DY(I+2) + DX(I+3)*DY(I+3) + DX(I+4)*DY(I+4)
   50 CONTINUE
   60 DDOT = DTEMP
      RETURN
      END
```

"""

# ╔═╡ 85135806-d982-4ca1-b852-fc53cfd0903e
md"""
## Speed of Computation
"""

# ╔═╡ 25d321c6-26bc-4d2d-ab74-8c048560e57c
function AB(A::Array{T},B::Array{T}) where T
    n=maximum(size(A))
    C=zeros(n,n)
    # C=[zeros(A[1,1]) for i=1:n, j=1:n]
    for i=1:n
        for j=1:n
            for k=1:n
                C[i,j]=C[i,j]+(A[i,k]*B[k,j])
            end
        end
    end
    C
end

# ╔═╡ 533c1833-8b02-4d45-b1fc-ad0255d0380a
import Random
Random.seed!(123)
n=512
A=rand(n,n)
B=rand(n,n)

# ╔═╡ 1783e655-bc2f-4442-a772-2b98f5ab2e96
# Run 2 times, the second measurement is relevant.
@time C=A*B

# ╔═╡ 698da56b-ed97-4e26-8456-403080e095e3
operations_in_second=(2*n^3)/0.003

# ╔═╡ ea827322-031d-42b5-aafd-10c2d12fa469
md"""
__Problem.__ Calculate largest $n$ for which three square $n\times n$ matrices of `Float64` numbers fit in your computer's RAM. How long will the multiplication last?  
"""

# ╔═╡ d1dccfbf-7878-4b48-a376-f0bf9b61ad59
# Run 2 times, the second measurement is relevant. 
# Our program is considerably slower!! 
@time AB(A,B);

# ╔═╡ 2015d3fc-37d2-43ee-99da-fe7029b73ee4
md"""
### Block variant

To speed-up our program by better usgae of cache memory, we need to compute with block-matrices (BLAS 3). 

BLAS level  |  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  operation &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  | memory | # ops
:---|:-----|:---|:---
1 | $y= \alpha x+y$  | $2n$ | $2n$
2 | $y=\alpha Ax+\beta y$ | $n^2$ | $2n^2$
3 | $C=\alpha A\cdot B+\beta C$| $3n^2$ | $2n^3$

__Problem.__ Study the subroutines `daxpy.f`, `dgemmv.f` i `dgemm.f`.
"""

# ╔═╡ be605463-3a97-4869-b9ba-426dc1822c96
function AB(A::Array{T},B::Array{T}) where T
    n=maximum(size(A))
    # C=zeros(n,n)
    C=[zero(A[1,1]) for i=1:n, j=1:n]
    for i=1:n
        for j=1:n
            for k=1:n
                C[i,j]+=(A[i,k]*B[k,j])
            end
        end
    end
    C
end

# ╔═╡ 6fc4fe72-e7ac-4c14-aa9c-37a5459d7004
# Try k,l=32,16 and k,l=64,8. Run twice.
k,l=64,8
Ab=[rand(k,k) for i=1:l, j=1:l]
Bb=[rand(k,k) for i=1:l, j=1:l];

# ╔═╡ 9bc2dc46-206d-46c5-b3c3-17d398d2cf6d
# This is considerably faster.
@time C=AB(Ab,Bb);

# ╔═╡ a7de05aa-0ed1-4e11-88de-de70d629d30b
md"""
For even better speed-up, we should use all cores, as `*` does. 

To run Julia in multi-threading environment in Linux insert the line

    export JULIA_NUM_THREADS=`nproc`

in  `.bashrc` file. For Windows, set the environment variable JULIA_NUM_THREADS to number 
of processors (run msinfo32).

"""

# ╔═╡ 209d4041-a0d9-455d-a28d-1a7cc632082f
Threads.nthreads()

# ╔═╡ 3b096734-bd85-4e21-ad81-6c1ed99e2f43
md"""
## Accuracy of computation

### Basic operations

For operations $\odot \in  \{+,*,/\}$ we have (for simplicity, we use $\epsilon$ for $\epsilon_M$)

$$
fl(a \odot b)=(1+\varepsilon_\odot)(a\odot b),\qquad |\varepsilon_\odot|\leq \varepsilon.
$$


### Addition

If we add fully accurately wto numbers which have (small) errors from previous computations, the equality 

$$ a(1+\varepsilon_a)+b(1+\varepsilon_b)= (1+\varepsilon_{ab})(a+b)
$$

gives

$$
\varepsilon_{ab}=\frac{a\varepsilon_a+b\varepsilon_a}{a+b},
$$

so  

$$
|\varepsilon_{ab}|\leq \varepsilon\, \frac{|a|+  |b|}{|a+b|}.
$$

If  $a$ and $b$ are large numbers of different signs and $a-b$ is tiny, the error can be huge (_catastrophic cancelation_ ).

### Scalar (dot) product

For vectors $x$ and $y$ from $\mathbb{R}^n$, the recursive application of the above fomula gives absolute error

$$|fl(x\cdot y)-x\cdot y|\leq O(n\varepsilon) |x|\cdot |y| \tag{1}$$

and relative error

$$\frac{|fl(x\cdot y)-x\cdot y|}{|x\cdot y|}\leq O(n\varepsilon) \frac{|x|\cdot |y|}{|x\cdot y|}$$

If vectors $x$ and $y$ are nearly perpendicular, relative error can be large.

### Matrix multiplication


From formula (1), for matrices $A$ and $B$ of order $n$, it follows

$$ |fl(A\cdot B) -A\cdot B| \leq O(n\varepsilon) |A|\cdot |B|.$$ 

"""

# ╔═╡ 1ff4efee-a106-46b4-855d-7b76265cb050
n=1_000_000
x=rand(n)
y=rand(n) .-0.5;

# ╔═╡ 1850009e-d407-46a7-be70-88fdb4a774be
d=x⋅y

# ╔═╡ b76c7b7c-2e9e-4985-a977-6ae52d0d4dcd
ab=abs.(x)⋅abs.(y)

# ╔═╡ 89027fca-0854-469b-947b-e68f62b3b8b7
# Check the solution using `BigFloat` numbers (70 decimal digits)
xbig=map(BigFloat,x)
ybig=map(BigFloat,y)
ybig[1]

# ╔═╡ c1387f1e-f3ce-49fd-8daa-d270b4714987
dbig=xbig⋅ybig

# ╔═╡ 815fca39-6bf8-44fc-ab5c-5146e444eb35
abserr=map(Float64,abs(d-dbig))

# ╔═╡ 24a79f79-ecc1-4dae-ab00-0e5be89d684f
relerr=map(Float64,abserr/d)

# ╔═╡ 23b3c19c-e8a0-4323-a501-f9b4fa53222b
n=256
A=rand(n,n)
B=rand(n,n).-0.5;

# ╔═╡ b4831348-c018-4133-907b-1a34987a9518
Ab=map(BigFloat,A)
Bb=map(BigFloat,B);

# ╔═╡ de747a60-3953-48ba-bf8c-4238beb07666
C=A*B
# This lasts little longer
Cb=Ab*Bb
abserr=abs.(C-map(Float64,Cb))

# ╔═╡ d43cea1b-5574-4341-aa91-c7ea2efc5e55
abs.(A)*abs.(B)*n*eps()

# ╔═╡ 2b7f29bf-a440-4f07-868b-ae42f3624b05
md"""
Here relative error is not easy to establish, so we need better measure. 
In the next lecture, we shall explain _vector and matrix norms_. 
"""

# ╔═╡ bf6a82b0-28a9-4df2-b967-7c87df059416


# ╔═╡ Cell order:
# ╟─1db27bd5-d11b-4cee-8f5e-53427cf82786
# ╠═d4ca4508-4d02-41f3-b343-e7d9e3d9349c
# ╠═a05d46aa-b091-4551-9919-c97b404cf4fa
# ╠═c60ff0fb-266e-43f2-926a-64913ad8cb3e
# ╠═c56ebccc-e0c3-4926-a7f2-6363bd0cc563
# ╠═242cd1ff-49f5-476f-a9c1-12fc4bccea04
# ╠═2ddf86e1-0751-4c0e-aee5-c3b1d4bed2fa
# ╠═da5dab82-8aa2-4bf3-9273-69bb3675508f
# ╠═7d8c987a-3a4e-48d9-8020-2170090f9934
# ╟─8cef7e23-c539-4508-9176-37132955257d
# ╟─05adbae2-88a2-49ed-9ca6-129757f94266
# ╟─6a085347-8f1d-4dc4-8742-774c4e279cfc
# ╟─b3f1858f-b1ec-4f60-8343-b4b12d2c60c7
# ╟─85135806-d982-4ca1-b852-fc53cfd0903e
# ╠═25d321c6-26bc-4d2d-ab74-8c048560e57c
# ╠═533c1833-8b02-4d45-b1fc-ad0255d0380a
# ╠═1783e655-bc2f-4442-a772-2b98f5ab2e96
# ╠═698da56b-ed97-4e26-8456-403080e095e3
# ╟─ea827322-031d-42b5-aafd-10c2d12fa469
# ╠═d1dccfbf-7878-4b48-a376-f0bf9b61ad59
# ╟─2015d3fc-37d2-43ee-99da-fe7029b73ee4
# ╠═be605463-3a97-4869-b9ba-426dc1822c96
# ╠═6fc4fe72-e7ac-4c14-aa9c-37a5459d7004
# ╠═9bc2dc46-206d-46c5-b3c3-17d398d2cf6d
# ╟─a7de05aa-0ed1-4e11-88de-de70d629d30b
# ╠═209d4041-a0d9-455d-a28d-1a7cc632082f
# ╟─3b096734-bd85-4e21-ad81-6c1ed99e2f43
# ╠═1ff4efee-a106-46b4-855d-7b76265cb050
# ╠═1850009e-d407-46a7-be70-88fdb4a774be
# ╠═b76c7b7c-2e9e-4985-a977-6ae52d0d4dcd
# ╠═89027fca-0854-469b-947b-e68f62b3b8b7
# ╠═c1387f1e-f3ce-49fd-8daa-d270b4714987
# ╠═815fca39-6bf8-44fc-ab5c-5146e444eb35
# ╠═24a79f79-ecc1-4dae-ab00-0e5be89d684f
# ╠═23b3c19c-e8a0-4323-a501-f9b4fa53222b
# ╠═b4831348-c018-4133-907b-1a34987a9518
# ╠═de747a60-3953-48ba-bf8c-4238beb07666
# ╠═d43cea1b-5574-4341-aa91-c7ea2efc5e55
# ╟─2b7f29bf-a440-4f07-868b-ae42f3624b05
# ╠═bf6a82b0-28a9-4df2-b967-7c87df059416
