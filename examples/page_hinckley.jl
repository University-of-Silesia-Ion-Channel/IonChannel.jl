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
	using .IonChannel
	import .IonChannel: plot, plot!, title!
	import Pkg
	Pkg.activate(".")
	using PlutoUI
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

# ╔═╡ d7f7458a-a23e-4fe9-872f-1c88efe47c5c
md"""
Pick data type:

$(@bind data_type Select(["txt", "pickle"]))
"""

# ╔═╡ 1dbb9c49-a92c-4f7a-87cd-72817cbb99cb
if data_type == "txt"
	data_names= cd(readdir, pwd() * "/$(data_folder)/sampling/");
elseif data_type == "pickle"
	data_names = cd(readdir, pwd() * ("/$(data_folder)/pickles/"));
end;

# ╔═╡ 1572f0e7-f6cc-424f-8292-0ee3d9b2f77d
if data_type == "txt"
	md"""
	Pick membrane voltage: $(@bind voltage Select(data_names))
	"""
else
	md"""
	Pick pickle type: $(@bind pickle_sub Select(data_names))
	"""
end

# ╔═╡ b277b495-5381-4591-af60-89ccc8aa80f6
begin
	if data_type == "txt"
		path_data = pwd() * "/$(data_folder)/sampling/$(voltage)/";
		path_dwell_times = pwd() * "/$(data_folder)/dwell_times/$(voltage)/";
		data_filenames = cd(readdir, path_data)[2:2:end];
		dwelltimes_filenames = cd(readdir, path_dwell_times)[1:2:end];
	elseif data_type == "pickle"
		path_data = pwd() * "/$(data_folder)/pickles/$(pickle_sub)/";
		path_dwell_times = ""
		data_filenames = cd(readdir, path_data);
	end
end;

# ╔═╡ 04dbc378-3548-4d7b-bb65-f9317c093195
data_filenames

# ╔═╡ 43f2d875-7e65-4e95-ae13-01c606180bcd
md"""
Data file: $(@bind data_file Select(data_filenames))
"""

# ╔═╡ 1f91879c-73cf-4575-b83d-fea8252e07bd
data_names

# ╔═╡ 52e04747-9111-427f-9dc1-4f33640c8eef
data_filenames

# ╔═╡ a5a2e939-8017-416e-9efe-5955ec467398
function get_starting_points_pickle(data_names, data_folder)
	
	
	what_first_dict = Dict()
	for d_name in data_names
		path_data = pwd() * "/$(data_folder)/pickles/$(d_name)/";
		path_dwell_times = ""
		data_filenames = cd(readdir, path_data)
		starting_points = Dict()
		for d_file in data_filenames
			@info "$d_file"
			# split_d_file = split(d_file, "_")
			data_file_path = path_data * d_file
			dwell_times_path = ""
			x, _ = read_data(data_file_path, dwell_times_path)
			starting_points[d_file] = x[1] == -1.0 ? 0 : 1
			@info "$(starting_points[d_file])"
		end
		what_first_dict[d_name] = starting_points
	end
	what_first_dict
end

# ╔═╡ c5c4380d-90b8-4ae4-826f-b687b08e30df
# ╠═╡ disabled = true
#=╠═╡
get_starting_points_pickle(data_names, data_folder)
  ╠═╡ =#

# ╔═╡ 1a751b4a-16a8-4270-9e4f-54c200b0a844
begin
	if data_type == "txt"
	local dt = split(data_file, '.')
	dt[1] = dt[1]*"dwell_timesy"
	md"""
	Dwell times file: $(dwell_times_file = join(dt, '.'))
	"""
	end
end

# ╔═╡ b9479dd4-a894-4c02-857f-facbfae2ab0e
begin
	if data_type == "txt"
		data_file_path = path_data * data_file
		dwell_times_path = path_dwell_times * dwell_times_file
	else
		data_file_path = path_data * data_file
		dwell_times_path = ""
	end
end

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

# ╔═╡ bd9f4e7c-0069-4b61-bd71-4e149c1f6aff
md"""
### Standardizing data
"""

# ╔═╡ 47620f98-a1b9-4d39-9911-26c27edfe4bf
Δt::Float32 = 1e-4

# ╔═╡ 8389df5c-2efe-47f6-8a9c-6e2a230b2b4f
begin
	x, y = read_data(data_file_path, dwell_times_path)
	if data_type == "pickle"
		y = Δt .* y .* 1000
	end
end

# ╔═╡ f66d33ae-2c3e-4661-84fe-8de5a29ae533
begin
	md"""
	Pick how many points to idealize (1000:$(length(x)))
	
	$(@bind data_size NumberField(1000:1000:length(x);default=50000))
	"""
	
end

# ╔═╡ 02435706-2495-49ff-b313-0021cfab445a
begin
	data = get_specified_datapoints(x, y, Δt, UInt32(data_size))
	normalized_data = normalize_data(data)
	data["x"] = normalized_data
end

# ╔═╡ 69406039-4a56-4163-b500-c13ad6618520
data

# ╔═╡ 9bae896a-a583-4fd9-8cc7-2b32ff748820
cumsum(data["dwell times"])

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

# ╔═╡ 10720e98-1f90-4cbd-b373-5d21b8729c58
data["dwell times"]

# ╔═╡ 5123c61b-7049-4af9-b2db-2800414e3ad2
md"""
## Implementation of Page Hinckley Method
"""

# ╔═╡ e1010167-91af-4670-bef4-3b2824789aec
md"""
Tolerance: $(@bind δ_ Slider(0.0:0.01:2.0; default=0.0, show_value=true))
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
m = MeanDeviationMethod(δ_)

# ╔═╡ 99f08837-ddcd-4eb3-9550-3aff9f257e8c
begin
	method_output = calculate_method(normalized_data, m, Δt)
	mean²error, h_dwell_times, h_dwell_times_approx = calculate_mean_square_error(data, method_output.dwell_times_approx)
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

# ╔═╡ e5f04386-4d5b-409e-80cd-11c010a32ffe
data_folder_path = pwd() * "/$(data_folder)"

# ╔═╡ c56f0f6e-5264-4220-85a4-8f8fae6f5ffb
md"""
Pick file to idealize data $(@bind what_first_path Select(cd(readdir, data_folder_path)))
"""

# ╔═╡ 511a0cd0-3e04-46cf-9a36-52dc7082dab2
begin
	what_fitst_file_path = pwd() * "/$(data_folder)/$(what_first_path)"
	what_first_dict = Dict(
	    String(split(line,',')[1]) => parse(UInt8, split(line,',')[2])
	    for line in eachline(what_fitst_file_path)
	)
end

# ╔═╡ 218d6150-eb02-4edb-84fa-0579bd539e5e
begin
	actual_idealized_data = actual_idealize_data(data, what_first_dict, data_file, Δt)
	accuracy_of_idealization(actual_idealized_data, method_output.idealized_data)
end

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

# ╔═╡ 1d8c4a52-ed04-4035-a2be-887e457368bc
begin
	hist = histogram_calculator(data["x"], Int16(100))
	prob_hist = calculate_probability_histogram(hist)
	analysis = analyze_histogram_peaks(prob_hist)
end

# ╔═╡ a16c263c-7439-49f0-86ad-dc772de5814c
plot(prob_hist)

# ╔═╡ a2378583-c520-4b38-a662-372a7594fc72


# ╔═╡ b92b0a30-db79-4361-865c-8a55e28852e9
begin
	dwell_times_distrib = IonChannel.fit(IonChannel.Exponential, data["dwell times"])
	plotrange = 0.0:0.001:maximum(vcat(data["dwell times"], method_output.dwell_times_approx))
	plot(plotrange, IonChannel.pdf(dwell_times_distrib, plotrange))
	dwell_times_approx_distrib = IonChannel.fit(IonChannel.Exponential, method_output.dwell_times_approx)
	plot!(plotrange, IonChannel.pdf(dwell_times_approx_distrib, plotrange))
end

# ╔═╡ ce812831-0e6e-4663-9bdd-d1fd26b5b81e
m2e = sum((IonChannel.pdf(dwell_times_distrib, plotrange) - IonChannel.pdf(dwell_times_approx_distrib, plotrange)).^2) / length(plotrange)

# ╔═╡ Cell order:
# ╟─91ef147a-729a-11f0-1157-03caaf19ff7b
# ╠═dbd814ae-a166-4096-a3bc-69a169aa1e5a
# ╟─4cf37e3f-6ac4-4a33-9b45-071c8d4e2347
# ╠═e4f89e4e-6d7e-40c2-9ce0-25924dcb6fee
# ╠═7044c2e9-cbd0-4b34-a476-3e03624a730c
# ╠═e5c455d4-6313-48ad-a2ba-e72b119429fe
# ╠═d7f7458a-a23e-4fe9-872f-1c88efe47c5c
# ╠═1dbb9c49-a92c-4f7a-87cd-72817cbb99cb
# ╠═1572f0e7-f6cc-424f-8292-0ee3d9b2f77d
# ╠═b277b495-5381-4591-af60-89ccc8aa80f6
# ╠═04dbc378-3548-4d7b-bb65-f9317c093195
# ╠═43f2d875-7e65-4e95-ae13-01c606180bcd
# ╠═1f91879c-73cf-4575-b83d-fea8252e07bd
# ╠═52e04747-9111-427f-9dc1-4f33640c8eef
# ╟─a5a2e939-8017-416e-9efe-5955ec467398
# ╠═c5c4380d-90b8-4ae4-826f-b687b08e30df
# ╠═1a751b4a-16a8-4270-9e4f-54c200b0a844
# ╟─3f81e4ec-8f4c-41de-be8e-ceef9ba69125
# ╠═b9479dd4-a894-4c02-857f-facbfae2ab0e
# ╠═8389df5c-2efe-47f6-8a9c-6e2a230b2b4f
# ╟─471b5827-a5cc-4cf0-9134-62b77750d473
# ╠═f66d33ae-2c3e-4661-84fe-8de5a29ae533
# ╠═bd9f4e7c-0069-4b61-bd71-4e149c1f6aff
# ╠═47620f98-a1b9-4d39-9911-26c27edfe4bf
# ╠═02435706-2495-49ff-b313-0021cfab445a
# ╠═69406039-4a56-4163-b500-c13ad6618520
# ╠═9bae896a-a583-4fd9-8cc7-2b32ff748820
# ╟─fef2cc53-d97a-4695-ac71-2caae03922d5
# ╠═18a20f62-284b-42ca-bad3-ebef333cfda8
# ╠═a61b3892-d50e-46ed-8a04-fcbe1f11e43e
# ╠═2cc62dc0-75d2-46a3-90c2-c3dd08c3a3e2
# ╠═f558c3bb-a336-455b-8397-8d8e8c6707c7
# ╠═a93aabed-114e-48a9-a322-f16579fd19f5
# ╠═0e60f2af-655f-4079-915f-656179159d08
# ╠═090447a8-0469-4199-a2b5-f898e44e40c3
# ╠═10720e98-1f90-4cbd-b373-5d21b8729c58
# ╟─5123c61b-7049-4af9-b2db-2800414e3ad2
# ╟─e1010167-91af-4670-bef4-3b2824789aec
# ╟─94dc7848-5852-4f9f-b7d7-1907aba39305
# ╟─90e2c2e2-652e-492f-beb1-58462466c98c
# ╠═4e60a773-afe6-4e84-8f31-d4e4c6342b02
# ╠═f2ac42ea-ebb0-42c4-9269-dccdcd575e74
# ╠═99f08837-ddcd-4eb3-9550-3aff9f257e8c
# ╠═e5f04386-4d5b-409e-80cd-11c010a32ffe
# ╠═c56f0f6e-5264-4220-85a4-8f8fae6f5ffb
# ╠═511a0cd0-3e04-46cf-9a36-52dc7082dab2
# ╠═218d6150-eb02-4edb-84fa-0579bd539e5e
# ╠═3b891242-97c6-4bfe-b4e2-bdb72639a688
# ╟─405c80b1-d86f-4139-acb8-57bb04513573
# ╠═1d8c4a52-ed04-4035-a2be-887e457368bc
# ╠═a16c263c-7439-49f0-86ad-dc772de5814c
# ╠═a2378583-c520-4b38-a662-372a7594fc72
# ╠═b92b0a30-db79-4361-865c-8a55e28852e9
# ╠═ce812831-0e6e-4663-9bdd-d1fd26b5b81e
