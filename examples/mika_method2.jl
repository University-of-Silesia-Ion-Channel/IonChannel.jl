### A Pluto.jl notebook ###
# v0.20.16

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 69b97b3d-80dd-41a6-b434-fd26fc2a6e39
begin
	include("../src/IonChannel.jl")
	using .IonChannel, PlutoUI
	import .IonChannel: plot, plot!, title!
end

# ╔═╡ 166394a2-736c-11f0-3403-4397f75a1ff3
md"""
## Loading necessary packages
"""

# ╔═╡ 59a7a871-c1ca-4ae3-be64-20ebe078eb53
md"""
## Step 1:  Calculation of the ion current’s probability distribution $\mathbb{P}$(I)
"""

# ╔═╡ ff9f76f5-3363-4c12-bb8d-74fa3efb422d
md"""
### Reading data from files
"""

# ╔═╡ 4963a611-376f-4411-a7a5-6a0ddaf1b01d
project_directory_files = cd(readdir, pwd());

# ╔═╡ b33321a8-40b2-4490-bad1-d0d484f72128
md"""
Pick data folder (has to be within the notebooks directory)

$(@bind data_folder Select(project_directory_files))
"""

# ╔═╡ d88ffeb9-5c1d-48d0-a121-b0bec21096c2
voltage_names= cd(readdir, pwd() * "/$(data_folder)/sampling/")

# ╔═╡ eb5160f5-90aa-402f-94b2-7e6fd7687b68
md"""
Pick membrane voltage: $(@bind voltage Select(voltage_names))
"""

# ╔═╡ ea94079f-d2a7-4324-8407-779e71ae9b3f
begin
	path_data = pwd() * "/$(data_folder)/sampling/$(voltage)/";
	path_dwell_times = pwd() * "/$(data_folder)/dwell_times/$(voltage)/";
end;

# ╔═╡ 8715a226-5efb-4883-845e-322b3f3bbf65
begin
	data_filenames = cd(readdir, path_data)
	clean_filenames = filter(
		fname -> occursin(r"^ce\d+\.txt$", fname),
		data_filenames
	)
end

# ╔═╡ 07ac8432-6ad0-4510-8b21-1100c75e84bc
md"""
Data file: $(@bind data_file Select(clean_filenames))
"""

# ╔═╡ c1366c40-186f-4b49-885d-7aeef515f9cf
begin
	local dt = split(data_file, '.')
	dt[1] = dt[1]*"dwell_timesy"
	md"""
	Dwell times file: $(dwell_times_file = join(dt, '.'))
	"""
end

# ╔═╡ 39edeabe-b2a7-47bc-aef5-66ded8be98e4
begin
	data_file_path = path_data * data_file
	dwell_times_path = path_dwell_times * dwell_times_file
end;

# ╔═╡ feb8b2a0-6845-43bc-8ddc-0d94b4d9feb1
x, y = read_data(data_file_path, dwell_times_path)

# ╔═╡ 85ba53f6-8432-4141-8f87-7ad9fa7f3070
begin
	md"""
	Pick how many points to idealize (1000:$(length(x)))
	
	$(@bind data_size NumberField(1000:length(x);default=50000))
	"""
	data_size = UInt32(data_size)
end

# ╔═╡ b38496ce-8e02-4d93-8b9a-ea31c2f5727e
begin
	Δt = Float32(1e-4)
	max_time = data_size*Δt
	data = get_specified_datapoints(x, y, Δt, data_size)
end

# ╔═╡ 0bd8638f-7d99-4fe7-b33c-7aae661704c2
md"""
### Standardizing data
"""

# ╔═╡ 2d20073a-5841-44dd-906d-d46c2aaacfa8
begin
	normalized_data = normalize_data(data)
	data["x"] = normalized_data
end

# ╔═╡ 65f9c674-d27f-4c5f-b819-a6ce49a531a4
md"""
### Creating histogram of the standardized data
"""

# ╔═╡ d46741a5-f0ac-4420-9d6a-c779e807916f
md"""
Number of bins $(@bind bins Slider(50:300; default=100, show_value=true))
"""

# ╔═╡ 7f5842e8-f126-494f-8623-979a00d67cd4
begin
	# histogram of the standardized data
	histogram_of_data = histogram_calculator(data["x"], bins)
end

# ╔═╡ b999775a-98a5-4381-a94e-a562c8198266
plot(histogram_of_data)

# ╔═╡ 41c55ea7-fe49-4b6b-a9e7-2ec7bddde8d2
md"""
### Computing probabilities to form the probability distribution $\mathbb{P}$(I)
"""

# ╔═╡ 6560cefa-9f7a-4248-8d27-337436db4f77
md"""
Choose the type of histogram (different visibility - nothing changes):
$(@bind histogram_type Select(["regular", "log"]))
"""

# ╔═╡ 77189139-305a-4815-bd92-ea2c360f57d4
probability_histogram = calculate_probability_histogram(histogram_of_data)

# ╔═╡ fefa50e6-5d1f-4a06-84cc-3a4ee7ba38ae
plot(probability_histogram)

# ╔═╡ 25f0f4a8-c40c-4113-913c-df5d13768cc5
md"""
## Step 2: Finding minimum value between peaks of a probability distribution histogram $\mathbb{P}$(I)
"""

# ╔═╡ 69bbd9b5-d5c4-4318-a2bc-912efc7b3bdc
histogram_analysis = analyze_histogram_peaks(probability_histogram)

# ╔═╡ e1f85d37-c06b-4ffc-94de-fdd0726d6529
md"""
## Step 3: Calculating dwell times using the minimum threshold
"""

# ╔═╡ 8928f413-52cc-4c1a-8b34-f3e75e0e030d
# ╠═╡ disabled = true
#=╠═╡
ϵ = 0.0
  ╠═╡ =#

# ╔═╡ c7b53fe7-4e49-4822-99a3-9bc3f04618dd
md"""
ϵ (The larger the noise the larger the ϵ and vice versa) 

$(@bind ϵ Slider(0.0:0.01:1; default=0.0, show_value=true))
"""

# ╔═╡ 2f98854f-902f-42d7-8d7a-4e045cf4b6ff
begin
	point_max1 = point(histogram_analysis, :pmax1_index, :pmax1)
	point_min = point(histogram_analysis, :pmin_index, :pmin)
	point_max2 = point(histogram_analysis, :pmax2_index, :pmax2)
end

# ╔═╡ 3dee373c-abfc-4a9f-9944-62f37fd7a7dc
begin
	line1 = line(point_max1, point_min)
	line2 = line(point_max2, point_min)
end

# ╔═╡ 5f38a022-197d-4919-858e-e4afddd04c74
md"""
## Step 4: Dividing the time series of ion channel into idealization + noise, using the minimum between peaks
"""

# ╔═╡ d5fc85f2-98cf-468b-8b72-e234a8d942ac
md"""
## Step 5: Optimizing the theshold current
"""

# ╔═╡ 3329d47a-f758-4d6e-84bc-dab7ee93b786
begin
	method = MikaMethod(ϵ, bins)
	optimized_data = calculate_method(normalized_data, method, Δt)
	# method_function(method)
end

# ╔═╡ 6ff04a66-a473-4e65-9f12-45485e4764a4
noise_test(optimized_data.noise)

# ╔═╡ 900280c6-f872-475f-90c2-3f73c812241b
md"""
## Checking the accuracy of the method
"""

# ╔═╡ 92eeaa84-9db4-41ba-88f0-899ff86ed8fb
begin
	state1_dwell_times_approx = dwell_times_approx(optimized_data)[1:2:end]
	state2_dwell_times_approx = dwell_times_approx(optimized_data)[2:2:end]
	state1_dwell_times = data["dwell times"][1:2:end]
	state2_dwell_times = data["dwell times"][2:2:end]
end

# ╔═╡ b60933b9-7eaa-4402-a74a-6520fde3f17a
begin
	md"""
	Number of bins $(@bind dt_bins Slider(10:150; default=100, show_value=true))
	"""
	dt_bins = UInt16(dt_bins)
end

# ╔═╡ f20d660a-5bab-4f60-aa60-35d084357c92
begin
	mean²error, h_dwell_times, h_dwell_times_approx = calculate_mean_square_error(data, dwell_times_approx(optimized_data), dt_bins)
end

# ╔═╡ 597c023d-e915-4846-bfd8-407e309c3852
md"""
Mean squared error $(mean²error)
"""

# ╔═╡ 89c354f9-341d-41d0-8e3c-c254003c37b7
md"""
Mean squared error $(mean²error)
"""

# ╔═╡ bdf93682-7c33-47f1-b93e-92eb9dacd4d4
begin
	plot(h_dwell_times; alpha=0.5, label="Exact")
	plot!(h_dwell_times_approx; alpha=0.5, label="Approximated")
	title!("Histogram of dwell times")
end

# ╔═╡ 5df65501-ecfa-4165-a317-3c5dcc66a0e8
md"""
Mean squared error $(mean²error)
"""

# ╔═╡ 3d798649-3698-496f-834e-fded75c2ea01
show_threshold_on_plot(probability_histogram, histogram_analysis, optimized_data)

# ╔═╡ e8c75a4d-a8ac-4ceb-a34e-92a09df58e85
sum(histogram_analysis.weights[histogram_analysis.pmax1_index:histogram_analysis.pmin_index])

# ╔═╡ ebad35d5-08fb-474a-9b82-8527691b8b17
sum(histogram_analysis.weights[histogram_analysis.pmin_index:histogram_analysis.pmax2_index])

# ╔═╡ a9955f5c-7672-411a-a341-c9484929a22a
md"""
## Show what's going on
"""

# ╔═╡ 79bf72be-62d7-43d7-a4e4-ad404ae88e75
begin
	md"""
	Left range index $(@bind N_left Slider(0:data_size; default=0.0, show_value=true))
	"""
end

# ╔═╡ 3fbc7bd1-9d53-432d-9489-13c25c862387
begin
	md"""
	Right range index $(@bind N_right Slider(N_left+500:data_size; default=100, show_value=true))
	"""
end

# ╔═╡ 27f17631-be59-4a81-a39d-4e510c214f98
T_left = trunc(N_left * Δt; digits=4)

# ╔═╡ 025a18f7-6e37-4724-851c-7e9d2c3755be
T_right = trunc(N_right *Δt; digits=4)

# ╔═╡ 39bfa301-3259-496f-89b3-d0b38c793e3b
plot_idealization_representation(data, optimized_data, T_left, T_right, Δt)

# ╔═╡ dd0870d6-f5a4-4127-998e-5e150e020535
show_approx_on_plot(data, optimized_data, T_left, T_right, Δt)

# ╔═╡ 32d0372c-772c-493d-9eb2-d57869ffb745
n = noise(optimized_data)

# ╔═╡ 8d96a05f-f1a1-4f80-a2ef-b81ad6c504ba
# ╠═╡ disabled = true
#=╠═╡
mean_error(method, Δt, UInt32(225000), true)
  ╠═╡ =#

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
PlutoUI = "~0.7.71"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.6"
manifest_format = "2.0"
project_hash = "0c76a76c3ac8f04e01e91e0dc955aee1f9d81e4a"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "67e11ee83a43eb71ddc950302c53bf33f0690dfe"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.12.1"

    [deps.ColorTypes.extensions]
    StyledStringsExt = "StyledStrings"

    [deps.ColorTypes.weakdeps]
    StyledStrings = "f489334b-da3d-4c2e-b8f0-e476e12c162b"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "b6d6bfdd7ce25b0f9b2f6b3dd56b2673a66c8770"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.5"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.6.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.7.2+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.11.0"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.MIMEs]]
git-tree-sha1 = "c64d943587f7187e751162b3b84445bbbd79f691"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "1.1.0"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.6+0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.12.12"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.27+1"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "7d2f8f21da5db6a806faf7b9b292296da42b2810"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.3"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.11.0"

    [deps.Pkg.extensions]
    REPLExt = "REPL"

    [deps.Pkg.weakdeps]
    REPL = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Downloads", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "8329a3a4f75e178c11c1ce2342778bcbbbfa7e3c"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.71"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "0f27480397253da18fe2c12a4ba4eb9eb208bf3d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.5.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

    [deps.Statistics.weakdeps]
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.Tricks]]
git-tree-sha1 = "372b90fe551c019541fafc6ff034199dc19c8436"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.12"

[[deps.URIs]]
git-tree-sha1 = "bef26fb046d031353ef97a82e3fdb6afe7f21b1a"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.6.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.59.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"
"""

# ╔═╡ Cell order:
# ╟─166394a2-736c-11f0-3403-4397f75a1ff3
# ╠═69b97b3d-80dd-41a6-b434-fd26fc2a6e39
# ╟─59a7a871-c1ca-4ae3-be64-20ebe078eb53
# ╟─ff9f76f5-3363-4c12-bb8d-74fa3efb422d
# ╠═4963a611-376f-4411-a7a5-6a0ddaf1b01d
# ╠═b33321a8-40b2-4490-bad1-d0d484f72128
# ╟─ea94079f-d2a7-4324-8407-779e71ae9b3f
# ╠═d88ffeb9-5c1d-48d0-a121-b0bec21096c2
# ╟─eb5160f5-90aa-402f-94b2-7e6fd7687b68
# ╠═8715a226-5efb-4883-845e-322b3f3bbf65
# ╠═07ac8432-6ad0-4510-8b21-1100c75e84bc
# ╟─c1366c40-186f-4b49-885d-7aeef515f9cf
# ╠═39edeabe-b2a7-47bc-aef5-66ded8be98e4
# ╟─597c023d-e915-4846-bfd8-407e309c3852
# ╠═feb8b2a0-6845-43bc-8ddc-0d94b4d9feb1
# ╠═85ba53f6-8432-4141-8f87-7ad9fa7f3070
# ╠═b38496ce-8e02-4d93-8b9a-ea31c2f5727e
# ╟─0bd8638f-7d99-4fe7-b33c-7aae661704c2
# ╠═2d20073a-5841-44dd-906d-d46c2aaacfa8
# ╟─65f9c674-d27f-4c5f-b819-a6ce49a531a4
# ╟─d46741a5-f0ac-4420-9d6a-c779e807916f
# ╠═7f5842e8-f126-494f-8623-979a00d67cd4
# ╠═b999775a-98a5-4381-a94e-a562c8198266
# ╟─41c55ea7-fe49-4b6b-a9e7-2ec7bddde8d2
# ╟─6560cefa-9f7a-4248-8d27-337436db4f77
# ╠═77189139-305a-4815-bd92-ea2c360f57d4
# ╠═fefa50e6-5d1f-4a06-84cc-3a4ee7ba38ae
# ╟─25f0f4a8-c40c-4113-913c-df5d13768cc5
# ╠═69bbd9b5-d5c4-4318-a2bc-912efc7b3bdc
# ╟─e1f85d37-c06b-4ffc-94de-fdd0726d6529
# ╠═8928f413-52cc-4c1a-8b34-f3e75e0e030d
# ╟─c7b53fe7-4e49-4822-99a3-9bc3f04618dd
# ╟─89c354f9-341d-41d0-8e3c-c254003c37b7
# ╠═6ff04a66-a473-4e65-9f12-45485e4764a4
# ╠═2f98854f-902f-42d7-8d7a-4e045cf4b6ff
# ╠═3dee373c-abfc-4a9f-9944-62f37fd7a7dc
# ╟─5f38a022-197d-4919-858e-e4afddd04c74
# ╟─d5fc85f2-98cf-468b-8b72-e234a8d942ac
# ╠═3329d47a-f758-4d6e-84bc-dab7ee93b786
# ╟─900280c6-f872-475f-90c2-3f73c812241b
# ╠═92eeaa84-9db4-41ba-88f0-899ff86ed8fb
# ╠═b60933b9-7eaa-4402-a74a-6520fde3f17a
# ╠═f20d660a-5bab-4f60-aa60-35d084357c92
# ╠═bdf93682-7c33-47f1-b93e-92eb9dacd4d4
# ╟─5df65501-ecfa-4165-a317-3c5dcc66a0e8
# ╠═3d798649-3698-496f-834e-fded75c2ea01
# ╠═e8c75a4d-a8ac-4ceb-a34e-92a09df58e85
# ╠═ebad35d5-08fb-474a-9b82-8527691b8b17
# ╟─a9955f5c-7672-411a-a341-c9484929a22a
# ╟─79bf72be-62d7-43d7-a4e4-ad404ae88e75
# ╟─3fbc7bd1-9d53-432d-9489-13c25c862387
# ╟─27f17631-be59-4a81-a39d-4e510c214f98
# ╟─025a18f7-6e37-4724-851c-7e9d2c3755be
# ╠═39bfa301-3259-496f-89b3-d0b38c793e3b
# ╠═dd0870d6-f5a4-4127-998e-5e150e020535
# ╠═32d0372c-772c-493d-9eb2-d57869ffb745
# ╠═8d96a05f-f1a1-4f80-a2ef-b81ad6c504ba
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
