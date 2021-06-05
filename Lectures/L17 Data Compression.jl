### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ b9c7ca7b-a5ca-4b9e-bb1c-5f237aa98894
begin
	import Pkg
	Pkg.activate(mktempdir())
    Pkg.add([
        Pkg.PackageSpec(name="PlutoUI"),
		Pkg.PackageSpec(name="Images"),
		Pkg.PackageSpec(name="Plots")
    ])
end

# ╔═╡ 1dd2098b-8d4d-4fa4-a794-49badf9790c2
begin
	using PlutoUI, Images, LinearAlgebra, Plots
	plotly()
end

# ╔═╡ c94e89d0-1049-4b0d-aa26-d5a1eda436b9
md"""
# Data Compression

QR factorization with column pivoting can be used for __compression of data__.

Diagonal elements of the matrix $R$ decrease in absolute value, so we can cut off parts of matrices $Q$ and $R$ which we deem unimportant.

Here is an example of image compression.
"""

# ╔═╡ 3d91d5e7-e134-464a-a0a4-032a28254715
# Create directory
if !isdir("files")
	mkdir("files")
end

# ╔═╡ ccbc114d-b5ff-4874-be91-008bc1b8e794
# Download the image
download("https://ivanslapnicar.github.io/NumericalMathematics/files/P8040001a.jpg",
	"files/P8040001a.jpg")

# ╔═╡ 0fa72a51-e9d4-4d23-a0a1-e477d2fb1611
# Load the image
img=load("files/P8040001a.jpg")

# ╔═╡ 464d999f-f530-4cc0-b9c8-fa094d8a138c
# Description of data
typeof(img)

# ╔═╡ 86818c56-6d28-463a-aa87-3c0a715ef8ce
# Very first pixel (upper left corner)
img[1,1]

# ╔═╡ ce428ce6-9bc5-4ce2-bb17-3bb8173fc937
# In numbers
show(img[1,1])

# ╔═╡ 1dfeac3f-51df-4fcf-919c-1675ae331189
# Split the image into  R, G and B component
channels=channelview(img)

# ╔═╡ c3ae4d9d-aa93-4ad8-af1d-fad775d9c9bd
begin
	Red=channels[1,:,:]
	Green=channels[2,:,:]
	Blue=channels[3,:,:]
end

# ╔═╡ c1cbc1a9-0a1e-4ac2-90fb-e9ce3c20be62
# View the blue channel (in gray)
colorview(Gray,Blue)

# ╔═╡ c0cd1d7e-1fb9-4ae8-aa8d-794c36b85837
begin
	# Compute QR factorization WITH pivoting of each channel
	R=qr(Red,Val(true))
	G=qr(Green,Val(true))
	B=qr(Blue,Val(true));
end

# ╔═╡ 32420d34-f066-4fdc-a488-15c4482640d8
# Residual
norm(R.Q*R.R[:,invperm(R.p)]-float(Red))

# ╔═╡ c16bd7ca-a403-4f1b-adba-f0e9a028b14b
# Plot diagonal elements of R
scatter(abs.(diag(R.R)),
    title="Diagonal elements of R",legend=false)

# ╔═╡ d7fd21c0-25ae-11eb-217c-51b165642c9b
md"
k = $(@bind k Slider(10:10:300, show_value=true, default=100))
"

# ╔═╡ feebc2d6-e495-4145-8e90-3b0d7d073d04
begin
	# Compute approximations of rank k for each channel, 
	# RedC, GreenC, and BlueC, respectively.
	# Function Matrix() is needed for faster manipulation with Q
	RedC=Matrix(R.Q)[:,1:k]*R.R[1:k,invperm(R.p)]
	GreenC=Matrix(G.Q)[:,1:k]*G.R[1:k,invperm(G.p)]
	BlueC=Matrix(B.Q)[:,1:k]*B.R[1:k,invperm(B.p)]
end

# ╔═╡ 955201fd-cb90-4c0f-bf50-ac1ca1850d26
# The compressed image
colorview(RGB, RedC, GreenC, BlueC)

# ╔═╡ 5bb5fe50-2466-11eb-175c-e1ff94b7840c
k

# ╔═╡ b1c18ee7-4b8f-4b44-a960-d759040a29cd
norm(Red-RedC)/norm(Red)

# ╔═╡ Cell order:
# ╠═b9c7ca7b-a5ca-4b9e-bb1c-5f237aa98894
# ╠═1dd2098b-8d4d-4fa4-a794-49badf9790c2
# ╟─c94e89d0-1049-4b0d-aa26-d5a1eda436b9
# ╠═3d91d5e7-e134-464a-a0a4-032a28254715
# ╠═ccbc114d-b5ff-4874-be91-008bc1b8e794
# ╠═0fa72a51-e9d4-4d23-a0a1-e477d2fb1611
# ╠═464d999f-f530-4cc0-b9c8-fa094d8a138c
# ╠═86818c56-6d28-463a-aa87-3c0a715ef8ce
# ╠═ce428ce6-9bc5-4ce2-bb17-3bb8173fc937
# ╠═1dfeac3f-51df-4fcf-919c-1675ae331189
# ╠═c3ae4d9d-aa93-4ad8-af1d-fad775d9c9bd
# ╠═c1cbc1a9-0a1e-4ac2-90fb-e9ce3c20be62
# ╠═c0cd1d7e-1fb9-4ae8-aa8d-794c36b85837
# ╠═32420d34-f066-4fdc-a488-15c4482640d8
# ╠═c16bd7ca-a403-4f1b-adba-f0e9a028b14b
# ╠═feebc2d6-e495-4145-8e90-3b0d7d073d04
# ╠═5bb5fe50-2466-11eb-175c-e1ff94b7840c
# ╠═955201fd-cb90-4c0f-bf50-ac1ca1850d26
# ╟─d7fd21c0-25ae-11eb-217c-51b165642c9b
# ╠═b1c18ee7-4b8f-4b44-a960-d759040a29cd
