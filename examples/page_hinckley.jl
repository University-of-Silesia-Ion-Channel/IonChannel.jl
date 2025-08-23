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

# ╔═╡ dbd814ae-a166-4096-a3bc-69a169aa1e5a
begin
	include("../src/IonChannel.jl")
	using .IonChannel, PlutoUI
	import .IonChannel: plot, plot!, title!
end

# ╔═╡ 91ef147a-729a-11f0-1157-03caaf19ff7b
md"""
## Loading necessary packages
"""

# ╔═╡ 4cf37e3f-6ac4-4a33-9b45-071c8d4e2347
md"""
# Page Hinkley Method
"""

# ╔═╡ e4f89e4e-6d7e-40c2-9ce0-25924dcb6fee
md"""
### Reading data from files
"""

# ╔═╡ 7044c2e9-cbd0-4b34-a476-3e03624a730c
project_directory_files = cd(readdir, pwd());

# ╔═╡ e5c455d4-6313-48ad-a2ba-e72b119429fe
md"""
Pick data folder (has to be within the notebooks directory)

$(@bind data_folder Select(project_directory_files))
"""

# ╔═╡ 1dbb9c49-a92c-4f7a-87cd-72817cbb99cb
voltage_names= cd(readdir, pwd() * "/$(data_folder)/sampling/");

# ╔═╡ 3ce29b12-01f6-48c6-95e6-233b2109d55f
voltage_names

# ╔═╡ 1572f0e7-f6cc-424f-8292-0ee3d9b2f77d
md"""
Pick membrane voltage: $(@bind voltage Select(voltage_names))
"""

# ╔═╡ b277b495-5381-4591-af60-89ccc8aa80f6
begin
	path_data = pwd() * "/$(data_folder)/sampling/$(voltage)/";
	path_dwell_times = pwd() * "/$(data_folder)/dwell_times/$(voltage)/";
end;

# ╔═╡ b22ddb25-23a9-4a3c-80a7-5b8984a95c1d
begin
	data_filenames = cd(readdir, path_data)[2:2:end];
	dwelltimes_filenames = cd(readdir, path_dwell_times)[1:2:end];
end;

# ╔═╡ 43f2d875-7e65-4e95-ae13-01c606180bcd
md"""
Data file: $(@bind data_file Select(data_filenames))
"""

# ╔═╡ 1a751b4a-16a8-4270-9e4f-54c200b0a844
begin
	local dt = split(data_file, '.')
	dt[1] = dt[1]*"dwell_timesy"
	md"""
	Dwell times file: $(dwell_times_file = join(dt, '.'))
	"""
end

# ╔═╡ b9479dd4-a894-4c02-857f-facbfae2ab0e
begin
	data_file_path = path_data * data_file
	dwell_times_path = path_dwell_times * dwell_times_file
end;

# ╔═╡ 8389df5c-2efe-47f6-8a9c-6e2a230b2b4f
x, y = read_data(data_file_path, dwell_times_path)

# ╔═╡ 471b5827-a5cc-4cf0-9134-62b77750d473
# ╠═╡ disabled = true
#=╠═╡
begin
	idealizations = create_idealizations()
	open("idealizations.txt", "w") do f
	    for (key, value) in idealizations
	        println(f, "$key")
			values = join(value, ", ")
			println(f, "$values")
	    end
	end
end
  ╠═╡ =#

# ╔═╡ f66d33ae-2c3e-4661-84fe-8de5a29ae533
begin
	md"""
	Pick how many points to idealize (1000:$(length(x)))
	
	$(@bind data_size NumberField(1000:1000:length(x);default=50000))
	"""
	data_size = UInt32(data_size)
end

# ╔═╡ a2bb9e0a-0566-468a-a4c9-81e6cbe00d86


# ╔═╡ bd9f4e7c-0069-4b61-bd71-4e149c1f6aff
md"""
### Standardizing data
"""

# ╔═╡ 47620f98-a1b9-4d39-9911-26c27edfe4bf
Δt::Float32 = 1e-4

# ╔═╡ 02435706-2495-49ff-b313-0021cfab445a
begin
	data = get_specified_datapoints(x, y, Δt, data_size)
	normalized_data = normalize_data(data)
	data["x"] = normalized_data
end

# ╔═╡ fef2cc53-d97a-4695-ac71-2caae03922d5
md"""
### Plotting data
"""

# ╔═╡ 18a20f62-284b-42ca-bad3-ebef333cfda8
begin
	md"""
	Left range index $(@bind N_left Slider(0:data_size; default=0.0, show_value=true))
	"""
end

# ╔═╡ a61b3892-d50e-46ed-8a04-fcbe1f11e43e
begin
	md"""
	Right range index $(@bind N_right Slider(N_left:data_size; default=N_left+500, show_value=true))
	"""
end

# ╔═╡ 2cc62dc0-75d2-46a3-90c2-c3dd08c3a3e2
T_left = trunc(N_left * Δt; digits=4)

# ╔═╡ f558c3bb-a336-455b-8397-8d8e8c6707c7
T_right = trunc(N_right * Δt ;digits=4)

# ╔═╡ 5123c61b-7049-4af9-b2db-2800414e3ad2
md"""
## Implementation of Page Hinckley Method
"""

# ╔═╡ e1010167-91af-4670-bef4-3b2824789aec
md"""
Tolerance: $(@bind δ_ Slider(0.0:0.01:2.0; default=0.0, show_value=true))

Detection parameter: $(@bind λ_ Slider(0.0:0.01:10.0; default=1.0, show_value=true))
"""

# ╔═╡ 90e2c2e2-652e-492f-beb1-58462466c98c
md"""
## Checking the accuracy of the method
"""

# ╔═╡ 4e60a773-afe6-4e84-8f31-d4e4c6342b02
begin
	md"""
	Number of bins $(@bind dt_bins Slider(10:150; default=100, show_value=true))
	"""
	dt_bins = UInt16(dt_bins)
end

# ╔═╡ f2ac42ea-ebb0-42c4-9269-dccdcd575e74
m = MeanDeviationMethod(δ_, λ_)

# ╔═╡ 99f08837-ddcd-4eb3-9550-3aff9f257e8c
begin
	method_output = calculate_method(normalized_data, m, Δt)
	mean²error, h_dwell_times, h_dwell_times_approx = calculate_mean_square_error(data, method_output.dwell_times_approx, dt_bins)
end

# ╔═╡ 3f81e4ec-8f4c-41de-be8e-ceef9ba69125
md"""
Mean squared error $(mean²error)
"""

# ╔═╡ a93aabed-114e-48a9-a322-f16579fd19f5
begin
	show_approx_on_plot(data, method_output, T_left, T_right, Δt)
end

# ╔═╡ 0e60f2af-655f-4079-915f-656179159d08
plot_idealization_representation(data, method_output, T_left, T_right, Δt)

# ╔═╡ 090447a8-0469-4199-a2b5-f898e44e40c3
method_output

# ╔═╡ 94dc7848-5852-4f9f-b7d7-1907aba39305
md"""
Mean squared error $(mean²error)
"""

# ╔═╡ 3b891242-97c6-4bfe-b4e2-bdb72639a688
begin
	plot(h_dwell_times; alpha=0.5, label="Exact")
	plot!(h_dwell_times_approx; alpha=0.5, label="Approximated")
	title!("Histogram of dwell times")
end

# ╔═╡ 405c80b1-d86f-4139-acb8-57bb04513573
md"""
Mean squared error $(mean²error)
"""

# ╔═╡ 512bfa00-d3e9-4eae-8d40-403dd96d17c8
# ╠═╡ disabled = true
#=╠═╡
mean_error(m, Δt, UInt32(225000), true)
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
# ╟─91ef147a-729a-11f0-1157-03caaf19ff7b
# ╠═dbd814ae-a166-4096-a3bc-69a169aa1e5a
# ╟─4cf37e3f-6ac4-4a33-9b45-071c8d4e2347
# ╠═e4f89e4e-6d7e-40c2-9ce0-25924dcb6fee
# ╠═7044c2e9-cbd0-4b34-a476-3e03624a730c
# ╠═e5c455d4-6313-48ad-a2ba-e72b119429fe
# ╠═1dbb9c49-a92c-4f7a-87cd-72817cbb99cb
# ╟─3ce29b12-01f6-48c6-95e6-233b2109d55f
# ╠═1572f0e7-f6cc-424f-8292-0ee3d9b2f77d
# ╠═b277b495-5381-4591-af60-89ccc8aa80f6
# ╠═b22ddb25-23a9-4a3c-80a7-5b8984a95c1d
# ╠═43f2d875-7e65-4e95-ae13-01c606180bcd
# ╠═1a751b4a-16a8-4270-9e4f-54c200b0a844
# ╟─3f81e4ec-8f4c-41de-be8e-ceef9ba69125
# ╠═b9479dd4-a894-4c02-857f-facbfae2ab0e
# ╠═8389df5c-2efe-47f6-8a9c-6e2a230b2b4f
# ╟─471b5827-a5cc-4cf0-9134-62b77750d473
# ╠═f66d33ae-2c3e-4661-84fe-8de5a29ae533
# ╠═a2bb9e0a-0566-468a-a4c9-81e6cbe00d86
# ╠═bd9f4e7c-0069-4b61-bd71-4e149c1f6aff
# ╠═47620f98-a1b9-4d39-9911-26c27edfe4bf
# ╠═02435706-2495-49ff-b313-0021cfab445a
# ╟─fef2cc53-d97a-4695-ac71-2caae03922d5
# ╠═18a20f62-284b-42ca-bad3-ebef333cfda8
# ╠═a61b3892-d50e-46ed-8a04-fcbe1f11e43e
# ╠═2cc62dc0-75d2-46a3-90c2-c3dd08c3a3e2
# ╠═f558c3bb-a336-455b-8397-8d8e8c6707c7
# ╠═a93aabed-114e-48a9-a322-f16579fd19f5
# ╠═0e60f2af-655f-4079-915f-656179159d08
# ╠═090447a8-0469-4199-a2b5-f898e44e40c3
# ╟─5123c61b-7049-4af9-b2db-2800414e3ad2
# ╟─e1010167-91af-4670-bef4-3b2824789aec
# ╟─94dc7848-5852-4f9f-b7d7-1907aba39305
# ╟─90e2c2e2-652e-492f-beb1-58462466c98c
# ╠═4e60a773-afe6-4e84-8f31-d4e4c6342b02
# ╠═f2ac42ea-ebb0-42c4-9269-dccdcd575e74
# ╠═99f08837-ddcd-4eb3-9550-3aff9f257e8c
# ╠═3b891242-97c6-4bfe-b4e2-bdb72639a688
# ╟─405c80b1-d86f-4139-acb8-57bb04513573
# ╠═512bfa00-d3e9-4eae-8d40-403dd96d17c8
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
