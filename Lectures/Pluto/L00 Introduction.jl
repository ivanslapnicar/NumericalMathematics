### A Pluto.jl notebook ###
# v0.12.4

using Markdown
using InteractiveUtils

# ╔═╡ a1bdf5a0-0ee6-11eb-16dc-c55bc6758a85
# Yet another package, follow instructions if not installed.
using Polynomials

# ╔═╡ a481f540-027d-11eb-10b5-5bd4e827e274
md"""
# Introduction

__Numerical analysis is the science of computing the solutions of problems that are posed
mathematically in the field of real or complex numbers.__

Let us just give a couple of examples.
"""

# ╔═╡ c2a8c4e0-027d-11eb-26fe-1f7bf12091f4
md"""
## Computation of $\pi \approx 3.14159265358979$

$\pi$ is known as the ratio of the circumference of a circle and its diameter.
"""

# ╔═╡ eab6ace0-027d-11eb-0f82-db588f0a253b
md"""
__Archimedes Method__

Inscribe a regular $n$-gon in a circle of radius $1$. Compute the perimeter of its
upper half. This is easiest to do if $n=2^k$. Using geometric reasoning,

$$
p(n) = 2n \sin \: \frac{\pi}{2n}.$$

Even though this refers to $\pi$ we can compute this value for $n=1,2,\cdots, $ without
computing $\pi$ and without punching the $\sin$ button on our calculator. (Computing
the $\sin$ function requires you to know $\pi$!)
"""

# ╔═╡ 3d725c70-027b-11eb-10cb-5b6e56311aa1
md"""
From geometric reasoning, we know that

$$
\sin \: \frac{\pi}{2} =1, \quad \sin \: \frac{\pi}{4} = \frac{\sqrt{2}}{2}.$$

So

$$
p(1)= 2 \cdot \sin \: \frac{\pi}{2} =2, \quad p(2) = 4 \cdot \sin \: \frac{\pi}{4} = 2\sqrt{2}= 2.828427125$$

but what is $p(4) = 8 \cdot \sin \: \displaystyle\frac{\pi}{8}$. Use half angle formulas!

\begin{eqnarray}
\sin \: \frac{\theta}{2} = \sqrt{\frac{1-\cos \: \theta }{2}}, \quad \cos \: \theta = \sqrt{1-\sin^2 \: \theta } .
\end{eqnarray}

That yields

$$
\sin \: \frac{\pi}{8} = 0.382683432.$$

Thus

$$
p(4) = 8 \cdot \sin \: \frac{\pi}{8}=3.061467459.$$

Continuing this process we get

n | $\sin\displaystyle\frac{\pi}{2n}$ | $p(n)$
:---|:---|:---
1 | 1 | 2
2 | $\sqrt{0.5}$ | 2.82842712
4 |  0.382683432 | 3.061467459
8 | 0.195090322 | 3.121445152
16 | 0.09801714 | 3.136548491
32 | 0.049067674 | 3.140331157

This method is slow, but "sure". 
"""

# ╔═╡ 3dea5041-85b9-4d8b-9dcf-2d39c8481314
steps=10

# ╔═╡ a47b53a0-027f-11eb-256d-f7d5ce0f97e3
let
	c=0.0
	for n=1:steps
    	s=sqrt((1-c)/2)
    	c=sqrt(1-s^2)
    	println("2n = ", 2^(n+1), ", približna vrijednost od π = ",2^(n+1)*s)
	end
end

# ╔═╡ 4d68f850-027b-11eb-14a1-1741ea4ce7c2
# Točan π
π

# ╔═╡ 6ac49527-05f0-4a0a-ad99-aba5d4b07c17
md"""
__Problem.__ Increase the number of steps and explain what happens.
"""

# ╔═╡ 8443676d-70db-4f63-8b8a-a7e96175f813
md"""
Later, we will give a modern enhancement that
is faster. Letting $h=1/(2n)$ we have that

$$
p(n) = \frac{\sin \: \pi h}{h} = \pi -a_2 h^2 + a_4 h^4 - \cdots$$

where $a_k = \pi^{k+1}/(k+1)!$. Thus this converges to $\pi$ at roughly the rate of
$O(h^2)$. This is an _approximation problem_. There is no finite algorithm
to compute $\pi$, since it is a transcendental number (irrational and not the root of
any polynomial with integer coefficients). However, we can approximate it arbitrarily
well.
"""

# ╔═╡ bd7a0ad8-8b31-42cf-94e9-881711a3c40e
md"""
## Quadratic equation

The following problem poses completely different issues.

Let us compute roots of $p(x) = ax^2 + bx +c =0$ for constants $a$, $b$ and $c$.

The solution is just

$$
x_{1,2} = \frac{ -b \pm \sqrt{b^2 - 4ac}}{2a}. \tag{1}$$

How do we go about writing a code to compute this.

Take care of special cases.

__Case I.__ $a=0, b \neq 0$.

It is no longer a quadratic, it is linear. Only solution is

$$
x_1 = -c/b.$$

__Case II.__ $a=b=0$

If $c \neq 0$, no solution. If $c=0$, all $x$ are solutions.

The cases you spent time on in high school had to do with the discriminant
$b^2 - 4ac$.

__Case III.__ $b^2 -4ac < 0$. Two Complex Roots (not real).

$$
x_{1,2} = -\frac{b}{2a} \pm \mathbf{i} \frac{\sqrt{4ac-b^2}}{2a}, \quad \mathbf{i}^2 = -1.$$

__Case IV.__ $b^2 - 4ac =0$. One Double Root (real)

$$
x_1 = x_2 = -\frac{b}{2a}.$$

__Case V.__ $b^2 -4ac > 0$. Two Distinct Real Roots.
Use Formula (1)?
"""

# ╔═╡ b5fe78a7-5258-41ba-9544-c0cd64496ef5
begin
	a=1
	b=2
	c=10.0^(-17)
	x₁=(-b-sqrt(b*b-4*a*c))/(2*a)
	x₂=(-b+sqrt(b*b-4*a*c))/(2*a)
	
	x₁,x₂
end

# ╔═╡ bca1eb35-dd3d-41af-8130-a8035551ec49
md"""
The two real roots are (to about 17 digits)

$$
x_1 = -2, \quad x_2 = -5 \cdot 10^{-18}.$$

The above algorithm gets $x_1$ right, but $x_2 = 0$. The standard double precision floating-point number format, `Float64`, stores about $15$ decimal digits (54 binary digits) and in those 15 digits $\sqrt{b^2 -4ac}-b =0$. 
That is because we are subtracting two close numbers and one of these is approximate, so this difference is "all rounding error". A simple observation gets around this problem.

The "large" root of the quadratic in __Case V__ is

$$
x_1 = \frac{ -b - \mathrm{sign}(b) \sqrt{b^2 - 4ac}}{2a}$$

and the two roots satisfy

$$
x_1 x_2 = \frac{c}{a}.$$

Notice that except inside the square root, we are adding numbers of the same sign!
After some algebra, we get a formula for the small root

$$
x_2 = \frac{c}{a x_1} = \frac{-2c}{ b + \mathrm{sign}(b) \sqrt{b^2 - 4ac}}.$$

Using this formula we compute both roots to near machine precision.

In this example, we have an exact formula, but in floating point arithmetic the standard quadratic formula
yields results that are significantly different from what it yields in real arithmetic.
"""

# ╔═╡ 2f084ae2-1e38-4c0b-83b6-da9f025565b1
begin
	x₃=(-b-sign(b)*sqrt(b*b-4*a*c))/(2*a)
	x₄=-2*c/(b+sign(b)*sqrt(b*b-4*a*c))
	x₃,x₄
end

# ╔═╡ f8cf034f-b149-4ff9-8a7e-bc8143581d51
md"""
The following function implements all five cases. Try out the function with different inputs and cover all five cases.
"""

# ╔═╡ 0738969f-7c4f-4f8e-92de-62cb0cd6510c
function quadroots(a,b,c)
    # Function to find the roots of the quadratic equation
    # given coefficients a,b, and c.
    # This function takes no account of scaling. 
    if a==0
        #  Check odd cases when a=0 
        if b==0
            if c==0
                return "all numbers are roots a=b=c=0"
            else
                return "no roots a=b=0, c ne 0"
            end
        else
            x₁=-c/b
            x₂=x₁
            ier="one root a=0"
        end
    else
        Δ = b*b-4*a*c
        if Δ < 0
            # Two complex roots computed with real arithmetic
            ximaginary=sqrt(-Δ)/(2*a)
            xreal=-b/(2*a)
            x₁=xreal+im*ximaginary
            # x₂ is the complex conjugate of x₁, 
            # x₂ = xreal - im*ximaginary
            x₂=conj(x₁)
        else
            if b==0
                # Since Julia handles complex arithmetic without a
                # blink, we can just use the formula in this case.
                x₁=sqrt(-c)/a
                x₂=-x₁
            else
                # Case where there are two real roots.
                x₁=(-b-sign(b)*sqrt(Δ))/(2*a)
                x₂=-2*c/(b+sign(b)*sqrt(Δ))
            end
        end
        ier="roots are good"
    end
    x₁, x₂, ier
end

# ╔═╡ cd29c4bc-9cfa-40bc-bc80-74ee5a635a94
quadroots(1,0,7)

# ╔═╡ af57717b-4dca-4b54-8272-e98fe76fd845
quadroots(a,b,c)

# ╔═╡ 8c4dfbc0-0ee6-11eb-2055-a7fdabff7984
md"""
## Julia is fast and open

`Julia` is __fast__ and is the first language to solve the "two language problem", for example using Matlab or Python for development and C, C++ or C# for speed. Every function is compiled after the first call using [JIT compiler](https://en.wikipedia.org/wiki/Just-in-time_compilation), in this case [LLVM](https://en.wikipedia.org/wiki/LLVM). 

`Julia` is  __open__ since
* The complete source code is always accessible on GitHub.
* MIT license.
* Macro `@which` makes orientation easy.
* LLVM and assembler code are easily displayed.
"""

# ╔═╡ b57ce420-0ee6-11eb-3909-9342f41a5255
# Contents of a package
varinfo(Polynomials)

# ╔═╡ c0a55620-0ee6-11eb-2297-d56f329d7bd1
?Polynomial

# ╔═╡ d3cbd530-0ee6-11eb-1a80-4ded54e128cd
p=Polynomial([c,b,a])

# ╔═╡ e3f0fb20-0ee6-11eb-0c86-cfbcb66da20f
roots(p)

# ╔═╡ ee73691e-0ee6-11eb-2953-23bd44f54b81
@which roots(p)

# ╔═╡ f9f65730-0ee6-11eb-0650-0b7c48826738
md"
You can also take a look at the corresponding GitHub file for the
most recent version.
"

# ╔═╡ 0c3d24a0-0ee7-11eb-1781-4b863e8378de
# Output is in Julia terminal
@code_llvm quadroots(1,0,7)

# ╔═╡ 1fe64130-0ee7-11eb-3db0-bbd29c571cda
# Output is in Julia terminal
@code_native quadroots(1,0,7)

# ╔═╡ Cell order:
# ╟─a481f540-027d-11eb-10b5-5bd4e827e274
# ╟─c2a8c4e0-027d-11eb-26fe-1f7bf12091f4
# ╟─eab6ace0-027d-11eb-0f82-db588f0a253b
# ╟─3d725c70-027b-11eb-10cb-5b6e56311aa1
# ╠═3dea5041-85b9-4d8b-9dcf-2d39c8481314
# ╠═a47b53a0-027f-11eb-256d-f7d5ce0f97e3
# ╠═4d68f850-027b-11eb-14a1-1741ea4ce7c2
# ╟─6ac49527-05f0-4a0a-ad99-aba5d4b07c17
# ╟─8443676d-70db-4f63-8b8a-a7e96175f813
# ╟─bd7a0ad8-8b31-42cf-94e9-881711a3c40e
# ╠═b5fe78a7-5258-41ba-9544-c0cd64496ef5
# ╟─bca1eb35-dd3d-41af-8130-a8035551ec49
# ╠═2f084ae2-1e38-4c0b-83b6-da9f025565b1
# ╟─f8cf034f-b149-4ff9-8a7e-bc8143581d51
# ╠═0738969f-7c4f-4f8e-92de-62cb0cd6510c
# ╠═cd29c4bc-9cfa-40bc-bc80-74ee5a635a94
# ╠═af57717b-4dca-4b54-8272-e98fe76fd845
# ╟─8c4dfbc0-0ee6-11eb-2055-a7fdabff7984
# ╠═a1bdf5a0-0ee6-11eb-16dc-c55bc6758a85
# ╠═b57ce420-0ee6-11eb-3909-9342f41a5255
# ╠═c0a55620-0ee6-11eb-2297-d56f329d7bd1
# ╠═d3cbd530-0ee6-11eb-1a80-4ded54e128cd
# ╠═e3f0fb20-0ee6-11eb-0c86-cfbcb66da20f
# ╠═ee73691e-0ee6-11eb-2953-23bd44f54b81
# ╟─f9f65730-0ee6-11eb-0650-0b7c48826738
# ╠═0c3d24a0-0ee7-11eb-1781-4b863e8378de
# ╠═1fe64130-0ee7-11eb-3db0-bbd29c571cda
