### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ de15c2d3-d314-47cc-9bd1-9c7f3109b92e
using PlutoUI, Polynomials

# ╔═╡ 7b2b10d7-46d6-406b-ab2b-951fbaedf8e9
TableOfContents(title="📚 Table of Contents", aside=true)

# ╔═╡ a481f540-027d-11eb-10b5-5bd4e827e274
md"""
# Introduction

__Numerical analysis is the science of computing the solutions of problems that are posed
mathematically in the field of real or complex numbers.__

Let us just give a couple of examples.
"""

# ╔═╡ c2a8c4e0-027d-11eb-26fe-1f7bf12091f4
md"""
# Computation of $\pi \approx 3.14159265358979$

 $\pi$ is known as the ratio of the circumference of a circle and its diameter.
"""

# ╔═╡ eab6ace0-027d-11eb-0f82-db588f0a253b
md"""
## Archimedes Method

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

$$
\sin \: \frac{\theta}{2} = \sqrt{\frac{1-\cos \: \theta }{2}}, \quad \cos \: \theta = \sqrt{1-\sin^2 \: \theta }.$$

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
steps=30

# ╔═╡ a47b53a0-027f-11eb-256d-f7d5ce0f97e3
let
	c=0.0
	for n=1:steps
    	s=sqrt((1-c)/2)
    	c=sqrt(1-s^2)
    	println("2n = ", 2^(n+1), ", approximation of π = ",2^(n+1)*s)
	end
end

# ╔═╡ 4d68f850-027b-11eb-14a1-1741ea4ce7c2
# Accurate π
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
$O(h^2)$. This is an __approximation problem__. There is no finite algorithm
to compute $\pi$, since it is a transcendental number (irrational and not the root of
any polynomial with integer coefficients). However, we can approximate it arbitrarily
well.
"""

# ╔═╡ bd7a0ad8-8b31-42cf-94e9-881711a3c40e
md"""
# Quadratic equation

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
# Julia is fast and open

`Julia` is __fast__ and is the first language to solve the "two language problem", for example using Matlab or Python for development and C, C++ or C# for speed. Every function is compiled after the first call using [JIT compiler](https://en.wikipedia.org/wiki/Just-in-time_compilation), in this case [LLVM](https://en.wikipedia.org/wiki/LLVM).

`Julia` is  __open__ since
* The complete source code is always accessible on GitHub.
* MIT license.
* Macro `@which` makes orientation easy.
* LLVM and assembler code are easily displayed.
"""

# ╔═╡ b57ce420-0ee6-11eb-3909-9342f41a5255
# Contents of package Polynomials
varinfo(Polynomials)

# ╔═╡ c0a55620-0ee6-11eb-2297-d56f329d7bd1
# ?Polynomial

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

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Polynomials = "f27b6e38-b328-58d1-80ce-0feddd5e7a45"

[compat]
PlutoUI = "~0.7.9"
Polynomials = "~2.0.14"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.1+0"

[[Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[ExprTools]]
git-tree-sha1 = "b7e3d17636b348f005f11040025ae8c6f645fe92"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.6"

[[FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[Intervals]]
deps = ["Dates", "Printf", "RecipesBase", "Serialization", "TimeZones"]
git-tree-sha1 = "323a38ed1952d30586d0fe03412cde9399d3618b"
uuid = "d8418881-c3e1-53bb-8760-2df7ec849ed5"
version = "1.5.0"

[[JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "8076680b162ada2a031f707ac7b4953e30667a37"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.2"

[[LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

[[Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.0+0"

[[Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[Mocking]]
deps = ["ExprTools"]
git-tree-sha1 = "748f6e1e4de814b101911e64cc12d83a6af66782"
uuid = "78c3b35d-d492-501b-9361-3d52fe80e533"
version = "0.7.2"

[[MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.2.1"

[[MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "3927848ccebcc165952dc0d9ac9aa274a87bfe01"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "0.2.20"

[[NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.20+0"

[[Parsers]]
deps = ["Dates"]
git-tree-sha1 = "438d35d2d95ae2c5e8780b330592b6de8494e779"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.0.3"

[[Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.8.0"

[[PlutoUI]]
deps = ["Base64", "Dates", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "Suppressor"]
git-tree-sha1 = "44e225d5837e2a2345e69a1d1e01ac2443ff9fcb"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.9"

[[Polynomials]]
deps = ["Intervals", "LinearAlgebra", "MutableArithmetics", "RecipesBase"]
git-tree-sha1 = "0bbfdcd8cda81b8144de4be8a67f5717e959a005"
uuid = "f27b6e38-b328-58d1-80ce-0feddd5e7a45"
version = "2.0.14"

[[Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[RecipesBase]]
git-tree-sha1 = "44a75aa7a527910ee3d1751d1f0e4148698add9e"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.1.2"

[[Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[Suppressor]]
git-tree-sha1 = "a819d77f31f83e5792a76081eee1ea6342ab8787"
uuid = "fd094767-a336-5f1f-9728-57cf17d0bbfb"
version = "0.2.0"

[[TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.0"

[[Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.1"

[[Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[TimeZones]]
deps = ["Dates", "Future", "LazyArtifacts", "Mocking", "Pkg", "Printf", "RecipesBase", "Serialization", "Unicode"]
git-tree-sha1 = "6c9040665b2da00d30143261aea22c7427aada1c"
uuid = "f269a46b-ccf7-5d73-abea-4c690281aa53"
version = "1.5.7"

[[UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.12+3"

[[libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.1.1+0"

[[nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"
"""

# ╔═╡ Cell order:
# ╠═de15c2d3-d314-47cc-9bd1-9c7f3109b92e
# ╠═7b2b10d7-46d6-406b-ab2b-951fbaedf8e9
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
# ╠═b57ce420-0ee6-11eb-3909-9342f41a5255
# ╠═c0a55620-0ee6-11eb-2297-d56f329d7bd1
# ╠═d3cbd530-0ee6-11eb-1a80-4ded54e128cd
# ╠═e3f0fb20-0ee6-11eb-0c86-cfbcb66da20f
# ╠═ee73691e-0ee6-11eb-2953-23bd44f54b81
# ╟─f9f65730-0ee6-11eb-0650-0b7c48826738
# ╠═0c3d24a0-0ee7-11eb-1781-4b863e8378de
# ╠═1fe64130-0ee7-11eb-3db0-bbd29c571cda
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
