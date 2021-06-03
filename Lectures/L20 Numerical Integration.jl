### A Pluto.jl notebook ###
# v0.12.18

using Markdown
using InteractiveUtils

# ╔═╡ 7d2eceff-5dc8-4088-8504-8fa469948983
using QuadGK, FastGaussQuadrature, LinearAlgebra

# ╔═╡ 456ae14f-d92c-4561-bcfd-d1b9849b8244
using FFTW

# ╔═╡ c42b0321-1003-4935-857e-a8422fabd485
md"""

# Numerical Integration

## Newton-Cotes formulas

Function $f(x):[a,b]\to\mathbb{R}$ is interpolated by the polynomial of degree $n$ through $n+1$ equally spaced points, and the integral is approximated by the integral of the interpolation polynomial. The polynomial $P_n(x)$ can be computed in Lagrange form (see the notebook [L10 Interpolating Functions.jl](https://ivanslapnicar.github.io/NumericalMathematics/L10%20Interpolating%20Functions.jl.html)): let 

$$
L_k(x)=\prod_{{i=0}\atop {i\neq k}}^n \frac{x-x_i}{x_k-x_i},$$

Then

$$
f(x)\approx P_n(x)=\sum_{k=0}^n f(x_k) L_k(x),$$

so

$$
\int_a^b f(x)\, dx\approx \int_a^b P_n(x) \, dx=\sum_{k=0}^n f(x_k) \int_a^b L_k(x)\, dx =(b-a)\sum_{k=0}^n \omega_k f(x_k). \qquad (1)$$

We have

$$x_i=a+\displaystyle\frac{b-a}{n}\, i.$$

The substitution $x=a+(b-a)\,t$ gives

$$
\frac{x-x_i}{x_k-x_i}=\frac{nt-i}{k-i},$$

so the  __weights__ $\omega_k$ are equal to

$$
\omega_k=\frac{1}{b-a}\int_a^b L_k(x)\, dx = \int_0^1 \prod_{{i=0}\atop {i\neq k}}^n \frac{nt-i}{k-i} \, dt. \qquad(2)$$
"""

# ╔═╡ 5071e7f3-d353-4429-8885-ca9431edf5a9
md"""
### Trapezoidal rule

For $n=1$,  formula (2) gives 

$$\omega_0=\omega_1=\frac{1}{2}.$$

The Newton-Cotesova formula (1) gives

$$\int_a^b f(x)\, dx\approx \int_a^b P_1(x) \, dx=\frac{b-a}{2}(f(a)+f(b)),$$

which is the area of the trapeziod with vertices $(a,0)$, $(a,f(a))$, $b,f(b))$, and $(b,0)$.

The precise meaning of $\approx$ is given by the following theorem:

__Theorem.__ If the function $f$ is continuous on the interval $[a,b]$, then there exists $c\in(a,b)$ such that

$$\int_a^b f(x)\, dx= \frac{b-a}{2}(f(a)+f(b))-\frac{(b-a)^3}{12}f''(c).$$

_Proof:_ According to Theorem from the notebook  [L10 Interpolating Functions.jl](https://ivanslapnicar.github.io/NumericalMathematics/L10%20Interpolating%20Functions.jl.html), there exists $c\in(a,b)$ such that

$$
E=\int_a^b f(x)\, dx -\int_a^b P_1(x) \, dx = \int_a^b (f(x)- P_1(x)) \, dx = \int_a^b \frac{f''(c)}{2}(x-a)(x-b)\, dx.$$

Using the substitution $x=a+(b-a)\,t$, we have

$$
E=\frac{f''(c)}{2} (b-a)^3 \int_0^1 t(t-1)\, dt= -\frac{(b-a)^3}{12}f''(c),$$

which completes the proof.

Let us divide the interval $[a,b]$ into $n$ sub-intervals of equal length,

$$[x_{i-1},x_{i}],\quad  i=1,2,\ldots,n,$$ 

and set

$$
\Delta x=\frac{b-a}{n}, \quad y_i=f(x_i).$$

On each sub-interval we have 

$$\int_{x_{i-1}}^{x_i} f(x)\, dx= \frac{\Delta x}{2}(y_{i-1}+y_i)
-\frac{(\Delta x)^3}{12}f''(c_i),\quad c_i\in[x_{i-1},x_i].$$

Summation yields

$$
\int_a^b f(x)\, dx=I_n+E_n,$$ 

where $I_n$ is the __Trapezoidal rule__,

$$
I_n=\Delta x\bigg( \frac{y_0}{2} +y_1+y_2+\cdots +y_{n-1}+\frac{y_n}{2}\bigg),$$

and the __error__ $E_n$ is

$$
E_n =-\frac{(\Delta x)^3}{12}\sum_{i=1}^n f''(c_i)=
-\frac{b-a}{12}(\Delta x)^2 f''(c), \quad c\in[a,b].$$

Here we have used the fact that the continuity of $f''$ on the interval $(a,b)$ implies the existence of the point $c\in[a,b]$ for which 

$$
\frac{1}{n}\sum_{i=1}^n f''(c_i)=f''(c).$$

Also, 

$$
|E_n|\leq \frac{b-a}{12}(\Delta x)^2 \max_{x\in(a,b)} |f''(x)|. \tag{3}$$


The derivation of the Trapeziodal rule and the error estimate can be found in the textbook [Numerical Mathematics and Computing, Section 5.2](https://web.ma.utexas.edu/CNA/NMC6/).

"""

# ╔═╡ 04f1b0f8-10bc-4aaf-bbf2-c411d740db0e
md"""
### Simpson's formula

For $n=2$, formula (2) gives 

$$\omega_0=\frac{1}{6},\quad \omega_1=\frac{2}{3},\quad \omega_2=\frac{1}{6}.$$

We divide $[a,b]$ into even number $n$ sub-intervals of equal length,

$$[x_{2i-1},x_{2i+1}],\quad i=1,2,\ldots,\frac{n}{2}.$$

We then apply Newton-Cotes formula (1) to each sub-interval, and sum, thus obtaining __Simpson's formula__:

$$
I_n=\frac{\Delta x}{3}\big(y_0 +4(y_1+y_3\cdots +y_{n-1})+2(y_2+y_4+\cdots+y_{n-2})+y_n\big).$$

It holds

$$
\int_a^b f(x)\, dx =I_n+E_n,$$

where the __error__ $E_n$ is bounded by

$$
|E_n|\leq \frac{b-a}{180}(\Delta x)^4 \max_{x\in(a,b)} |f^{(4)}(x)|.$$

Details can be found in the textbook [Numerical Mathematics and Computing, Section 6.1](https://web.ma.utexas.edu/CNA/NMC6/).

"""

# ╔═╡ 7a11ccd0-68c1-46d3-9249-3cc9f867ef65
md"""
### Richardson's extrapolation

Estimation of errors using formulas (3) and (4) can be complex. Provided certain conditions are met, the __Richardson's extrapolation__ enables us to estimate the error using approximation of the integral with $n/2$ points. If the error bound contains the term $(\Delta x)^m$, then the error is approximated by

$$
E=\frac{\big(\frac{n}{2}\big)^m}{n^m-\big(\frac{n}{2}\big)^m}(I_n-I_{n/2}).$$

Here the sign of $E$ is also the sign of the error, that is, if $E>0$, then (approximately)

$$
\int_a^b f(x)\, dx\in[I_n,I_n+E],$$

and if $E\leq 0$, then (approximately)

$$
\int_a^b f(x)\, dx\in[I_n+E,I_n].$$

To prove the above formulas, we assume that there exists $\omega$ such that 

$$
\displaystyle I=I_n+E_n, \qquad E_n=\omega\,  (\Delta x)^m \tag{5}$$

for every $ \Delta x$. This assumption is not always exactly fulfilled, but in many cases we can assume it holds. 
For example, in Trapezoidal rule we can set

$$\omega=-\frac{b-a}{12}\max_{x\in(a,b)}f''(x),$$

so the assumption (5) means that we can take (approximately) the same $\omega$ for different values of $\Delta x$. Using the assumption, we have

$$
\begin{aligned}
E_{n}&=I-I_n =\omega (\Delta x)^m =\omega \, \left(\frac{b-a}{n}\right)^m,\cr
E_{n/2}&=I-I_{n/2} =\omega (2\Delta x)^m =\omega\,  \left(\frac{b-a}{\frac{n}{2}}\right)^m.
\end{aligned}$$

Thus,

$$
\displaystyle I_{n}-I_{n/2}=\omega \, (b-a)^m \left(\frac{1}{(n/2)^m}-\frac{1}{n^m}\right),$$

so

$$
\displaystyle \omega = \frac{n^m (n/2)^m}{(b-a)^m} \frac{I_{n}-I_{n/2}}{n^m-(n/2)^m},$$

which finally yields

$$
\displaystyle E_{n}\approx \frac{(n/2)^m} {n^m-(n/2)^m} (I_{n}-I_{n/2})=E.$$

"""

# ╔═╡ e9054e4d-0e5a-4b26-a323-895ec36ddf73
function Trapezoid(f::Function,a::Number,b::Number,n::Int64)
    # n is the number of intervals
    m=2
    X=range(a,stop=b,length=n+1)
    Y=map(f,X)
    Δx=(b-a)/n
    I=Δx*(Y[1]/2+sum(Y[2:end-1])+Y[end]/2)
    # Richardson's extrapolation
    I₂=2*Δx*(Y[1]/2+sum(Y[3:2:end-2])+Y[end]/2)
    E=(n/2)^m*(I-I₂)/(n^m-(n/2)^m)
    I,E
end 

# ╔═╡ 0283137b-539e-43c2-b2e4-34f1890d78b9
function Simpson(f::Function,a::Number,b::Number,n::Int64)
    # n is the number of intervals, divisible by 4 
    m=4
    X=range(a,stop=b,length=n+1)
    Y=map(f,X)
    Δx=(b-a)/n
    I=Δx/3*(Y[1]+4*sum(Y[2:2:end-1])+2*sum(Y[3:2:end-2])+Y[end])
    # Richardson's extrapolation
    I₂=2*Δx/3*(Y[1]+4*sum(Y[3:4:end-2])+2*sum(Y[5:4:end-4])+Y[end])
    E=(n/2)^m*(I-I₂)/(n^m-(n/2)^m)
    I,E
end 

# ╔═╡ 4baefc88-ba8b-43f7-8cec-48bf4d80c314
md"""
### Elliptic integral

Let us compute the circumference of the ellipse with semi-axes $2$ and $1$. Parametric equations of the ellipse are

$$
x=2\cos t,\quad y=\sin t,\quad t\in[0,\pi/2],$$

so the quarter of the circumference is

$$
\frac{O}{4}\int\limits_0^{\pi/2} \sqrt{(-2\sin t)^2+(\cos t)^2}\, dt
=2\int\limits_0^{\pi/2} \sqrt{1-\frac{3}{4}(\cos t)^2}\, dt.$$

The integral on the right-hand side is the __elliptic integral of the second kind__ which is not elementary solvable, but can be found 
[tabulated](http://nvlpubs.nist.gov/nistpubs/jres/50/jresv50n1p43_A1b.pdf). We see that 

$$
C=\int\limits_0^{\pi/2} \sqrt{1-\frac{3}{4}(\cos t)^2}\, dt \approx 8\cdot 1.21125.$$

"""

# ╔═╡ d0135ecd-9e39-4233-ab0d-74ac85ad4e2f
begin
	f₁(x)=sqrt(1-(3.0)/4*cos(x)^2)
	Trapezoid(f₁,0,π/2,4)
end

# ╔═╡ c0ba2cf8-7751-4bf2-8fd0-c4bb69e95504
Trapezoid(f₁,0,pi/2,10)

# ╔═╡ e370021d-2f28-4ec2-85c6-6bb094d242d1
Trapezoid(f₁,0,pi/2,24)

# ╔═╡ faab641c-5a27-4e36-a4dd-c8e82272eae9
Simpson(f₁,0,π/2,4)

# ╔═╡ c58f0660-6dd6-4045-bb33-a74351c48f73
Simpson(f₁,0,π/2,16)

# ╔═╡ a7cb847f-bc7d-4299-9b13-37a698a1f28b
Simpson(f₁,0,π/2,24)

# ╔═╡ aa07f01e-3bb9-4b47-b197-aba3c1310a2c
md"""
### The number $\pi$

It holds

$$
\int_0^1 \frac{4}{1+x^2}\, dx=\pi.$$

Let us approximate $\pi$ using numerical integration and check the error.
Using the Trapezoid rule we can obtain at most five significant digits. The Simpson's formula is more accurate, but slower. 
"""

# ╔═╡ 311b6a2a-021b-4fc2-b9b0-5bafa0b1d4a6
# For comparison
BigFloat(π)

# ╔═╡ 1c13d010-afbf-45b4-8210-a7cfea8641a4
begin
	f₂(x)=4/(1+x^2)
	myπ=Trapezoid(f₂,0,1,10)
end

# ╔═╡ 94267bb3-624e-49b4-bbc8-2c24aacbd095
myπ[1]-BigFloat(π)

# ╔═╡ cb4df66a-1b00-403f-b49b-7a275dc2929f
myπ₁=Trapezoid(f₂,0,1,100)

# ╔═╡ 5cb0017b-1722-4f38-9f53-76ff788fce9e
myπ₁[1]-BigFloat(π)

# ╔═╡ d1b0fca4-b727-4670-ad25-3545dfb62629
myπ₂=Simpson(f₂,0,1,16)

# ╔═╡ 29e800e1-b38e-4dbf-a1f1-6084b35ef762
myπ₂[1]-BigFloat(π)

# ╔═╡ ad0d05f0-8e66-4410-956b-bf6802ada45a
myπ₃=Simpson(f₂,0,1,64)

# ╔═╡ 11c3d407-0923-4d65-85a2-b30911fd75d2
myπ₃[1]-BigFloat(π)

# ╔═╡ e7b2a060-7ae9-4418-a434-1693f73d04ac
md"""
## Gaussian quadrature

Similarly to formula (1), the integral is approximated by a sum of products of function values and corresponding weights:

$$
\int_{a}^b \omega(x) f(x)\, dx=\sum_{k=1}^n \omega_k f(x_k),$$

where $\omega(x)$ is the __weight function__.

The points $x_k$ are the zeros of the corresponding orthogonal polynomial $P_{n}(x)$ of degree $n$, for example, __Legendre polynomial__ for $[a,b]=[-1,1]$ and $\omega(x)=1$, or __Chebyshev polynomial__ for 
$[a,b]=[-1,1]$ and $\omega(x)=\displaystyle\frac{1}{\sqrt{1-x^2}}$ (see the notebook  [L12 Orthogonal Polynomials.ipynb](https://nbviewer.jupyter.org/github/ivanslapnicar/NumericalMathematics/blob/master/Lectures/Jupyter/L12%20Orthogonal%20Polynomials.ipynb)).

The __weights__ are

$$
\omega_k=\int_a^b \omega(x) \prod_{{i=1}\atop {i\neq k}}^n\frac{x-x_i}{x_k-x_i} \, dx.$$

The __error__ is

$$
E=\frac{f^{(2n)}(\xi)}{(2n)!}\int_a^b \omega(x) P_n^2(x)\, dx.$$

__Remark__. Legendre and Chebyshev polynomials are defined on the interval $[-1,1]$, so for other intervals we use the substitution

$$
\int_{a}^b \omega(x) f(x)\, dx = \frac{b-a}{2} \int_{-1}^1 \omega\bigg(\frac{b-a}{2}t+\frac{a+b}{2}\bigg) f\bigg(\frac{b-a}{2}t+\frac{a+b}{2}\bigg) dt. \tag{6}$$

Details can be found in the textbook 
[Numerical Mathematics and Computing, Section 6.2](https://web.ma.utexas.edu/CNA/NMC6/).
"""

# ╔═╡ 656b359d-9504-4fa3-a327-d0977e10d642
mapnodes(x,a,b)=(b-a)*x/2 .+(a+b)/2

# ╔═╡ 427bfd71-5f3d-4b0b-b6f8-3753b5844a87
md"""
### Existing routines

Professional routines for numerical integration are rather complex, and majority of the programs have some of them built-in:

* Matlab's command `quad` uses adaptive Simpson's formula. 
* Julia's package [`QuadGK.jl`](https://github.com/JuliaMath/QuadGK.jl) contains the function `quadgk()` which approximates the integral using $O(n^2)$ operations with Gauss-Kronrod method.
* Julia also has the package [`FastGaussQuadrature.jl`](https://github.com/ajt60gaibb/FastGaussQuadrature.jl) which computes points and weights for a given number of points $n$ and the given weight function. Using these points and weights, the integral is easily approximated by dot product in $O(n)$ operations.
"""

# ╔═╡ cb4a8760-9ba7-4bb2-beea-9876a680829f
# ?quadgk

# ╔═╡ 4feab23e-b51b-4818-8343-a73deb67c942
# 1/8 of the circumference of elipse
quadgk(f₁,0,π/2)

# ╔═╡ fa5a8850-c845-4328-94ca-cef66b4b52dd
# Number π
quadgk(f₂,0,1)

# ╔═╡ ac14acbb-5e64-44ca-9922-68bd3a45b0a2
# Interval can be infinite
quadgk(x->exp(-x),0,Inf)

# ╔═╡ 4c7f59a4-0e11-4dd1-8d08-3c12b25223f5
varinfo(FastGaussQuadrature)

# ╔═╡ 5fbe0581-65b8-44b4-94f5-35950d579c92
# For example
methods(gausschebyshev)

# ╔═╡ b9a58890-f3a2-49bb-836c-95a2971c9da7
gausschebyshev(16)

# ╔═╡ 4ae529de-6df8-4d38-917f-7b972904204d
# We compute the integrals. In our case ω(x)=1,
# so we need Legendre polynomial.
ξ,ω=gausslegendre(32)

# ╔═╡ 988b6a6c-7732-47cd-915d-a26f44e1b3f0
# 1/8 of the circumference of the ellipse; a=0, b=π/2
(π/2-0)/2*dot(ω,map(f₁,mapnodes(ξ,0,π/2)))

# ╔═╡ 5fa9e09b-569d-4ad4-8dd6-ffc93dc6cca9
# Number π; a=0, b=1
(1-0)/2*dot(ω,map(f₂,mapnodes(ξ,0,1)))

# ╔═╡ 75fe9d8a-5e70-4d31-98f2-0ea82c511698
md"""
## Clenshaw-Curtis quadrature

With the substitution $x=\cos\theta$, it holds

$$
I\equiv \int\limits_{-1}^1f(x)\, dx =\int\limits_0^\pi f(\cos\theta)\sin\theta \, d\theta.$$

The integral on the right hand side is approximated by integrating Fourier series of the even (continuous) extension of the integrand:

$$
I\approx a_0+\sum_{k=1}^n \frac{2a_{2k}}{1-(2k)^2},$$

where the coefficients $a_k$ are computed as

$$
a_k=\frac{2}{\pi}\int\limits_0^\pi f(\cos\theta)\cos(k\theta)\,d\theta.$$

The coefficients can be computed using numerical integration or using Fast Fourier Transform (FFT), which is a lot faster. For details see [Clenshaw-Curtis Quadrature](https://en.wikipedia.org/wiki/Clenshaw%E2%80%93Curtis_quadrature).

The integration interval $[a,b]$ is transformed to the interval $[-1,1]$ using formula (6).
"""

# ╔═╡ 3d03b422-76c1-4f9b-ad69-12c07f2d4403
function ClenshawCurtis(f::Function,a::Number,b::Number,n::Int64)
    # Implementation using numerical integration
    z=Vector{Float64}(undef,n)
    g(x)=f(mapnodes(x,a,b))
    for i=1:n
        h(x)=g(cos(x))*cos(2*(i-1)*x)
        z[i]=2*quadgk(h,0,pi)[1]/pi
    end
    return (z[1]+2*sum([z[i]/(1-4*(i-1)^2) for i=2:n]))*(b-a)/2
end

# ╔═╡ cc0b62ab-821c-4b16-9c73-ec131176e296
# 1/8 of the circumference of the ellipse
ClenshawCurtis(f₁,0,pi/2,8)

# ╔═╡ 6ab103b9-e138-4c57-8730-83fe6313bbb0
# Number π, try 8 and 16
ClenshawCurtis(f₂,0,1,8)

# ╔═╡ 35dbcd72-9f1d-43a2-9134-4116eefacc97
 π

# ╔═╡ 307e81fa-5384-4aad-8fe2-7df9ab5bce0d
# "Improper" integral
ClenshawCurtis(x->exp(-x),0,1000,50)

# ╔═╡ 748ace02-d80f-4952-883b-de220083d58a
function ClenshawCurtisFFT(f::Function,a::Number,b::Number,n::Int64)
    # Fast implementation using fft(), 2^n is the number of points
    g(x)=f(mapnodes(x,a,b))
    w=map(x->g(cos(x)),range(0,stop=2*pi,length=2^n))
    w[1]=(w[1]+w[end])/2
    z=real(fft(w))
    z/=2.0^(n-1)
    return (z[1]+2*sum([z[i]/(1-(i-1)^2) for i=3:2:2^(n-1)]))*(b-a)/2
end

# ╔═╡ 44406f97-1a7e-490b-9e75-3ab42e7f353f
ClenshawCurtisFFT(f₁,0,pi/2,4)

# ╔═╡ 41654b34-29c9-44d9-a0a2-d30a77d1b29c
ClenshawCurtisFFT(f₂,0,1,16),pi

# ╔═╡ 15905e38-d8f3-450a-9bcc-1d6b7a4fcae1
ClenshawCurtisFFT(x->exp(-x),0,1000,18)

# ╔═╡ a46784f5-e834-4671-9689-fd7ec3fada6b


# ╔═╡ Cell order:
# ╟─c42b0321-1003-4935-857e-a8422fabd485
# ╟─5071e7f3-d353-4429-8885-ca9431edf5a9
# ╟─04f1b0f8-10bc-4aaf-bbf2-c411d740db0e
# ╟─7a11ccd0-68c1-46d3-9249-3cc9f867ef65
# ╠═e9054e4d-0e5a-4b26-a323-895ec36ddf73
# ╠═0283137b-539e-43c2-b2e4-34f1890d78b9
# ╟─4baefc88-ba8b-43f7-8cec-48bf4d80c314
# ╠═d0135ecd-9e39-4233-ab0d-74ac85ad4e2f
# ╠═c0ba2cf8-7751-4bf2-8fd0-c4bb69e95504
# ╠═e370021d-2f28-4ec2-85c6-6bb094d242d1
# ╠═faab641c-5a27-4e36-a4dd-c8e82272eae9
# ╠═c58f0660-6dd6-4045-bb33-a74351c48f73
# ╠═a7cb847f-bc7d-4299-9b13-37a698a1f28b
# ╟─aa07f01e-3bb9-4b47-b197-aba3c1310a2c
# ╠═311b6a2a-021b-4fc2-b9b0-5bafa0b1d4a6
# ╠═1c13d010-afbf-45b4-8210-a7cfea8641a4
# ╠═94267bb3-624e-49b4-bbc8-2c24aacbd095
# ╠═cb4df66a-1b00-403f-b49b-7a275dc2929f
# ╠═5cb0017b-1722-4f38-9f53-76ff788fce9e
# ╠═d1b0fca4-b727-4670-ad25-3545dfb62629
# ╠═29e800e1-b38e-4dbf-a1f1-6084b35ef762
# ╠═ad0d05f0-8e66-4410-956b-bf6802ada45a
# ╠═11c3d407-0923-4d65-85a2-b30911fd75d2
# ╟─e7b2a060-7ae9-4418-a434-1693f73d04ac
# ╠═656b359d-9504-4fa3-a327-d0977e10d642
# ╟─427bfd71-5f3d-4b0b-b6f8-3753b5844a87
# ╠═7d2eceff-5dc8-4088-8504-8fa469948983
# ╠═cb4a8760-9ba7-4bb2-beea-9876a680829f
# ╠═4feab23e-b51b-4818-8343-a73deb67c942
# ╠═fa5a8850-c845-4328-94ca-cef66b4b52dd
# ╠═ac14acbb-5e64-44ca-9922-68bd3a45b0a2
# ╠═4c7f59a4-0e11-4dd1-8d08-3c12b25223f5
# ╠═5fbe0581-65b8-44b4-94f5-35950d579c92
# ╠═b9a58890-f3a2-49bb-836c-95a2971c9da7
# ╠═4ae529de-6df8-4d38-917f-7b972904204d
# ╠═988b6a6c-7732-47cd-915d-a26f44e1b3f0
# ╠═5fa9e09b-569d-4ad4-8dd6-ffc93dc6cca9
# ╟─75fe9d8a-5e70-4d31-98f2-0ea82c511698
# ╠═3d03b422-76c1-4f9b-ad69-12c07f2d4403
# ╠═cc0b62ab-821c-4b16-9c73-ec131176e296
# ╠═6ab103b9-e138-4c57-8730-83fe6313bbb0
# ╠═35dbcd72-9f1d-43a2-9134-4116eefacc97
# ╠═307e81fa-5384-4aad-8fe2-7df9ab5bce0d
# ╠═456ae14f-d92c-4561-bcfd-d1b9849b8244
# ╠═748ace02-d80f-4952-883b-de220083d58a
# ╠═44406f97-1a7e-490b-9e75-3ab42e7f353f
# ╠═41654b34-29c9-44d9-a0a2-d30a77d1b29c
# ╠═15905e38-d8f3-450a-9bcc-1d6b7a4fcae1
# ╠═a46784f5-e834-4671-9689-fd7ec3fada6b
