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
	import .IonChannel: plot, plot!, title!, hline!, vline!
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
# MDL Method
"""

# ╔═╡ e4f89e4e-6d7e-40c2-9ce0-25924dcb6fee
md"""
### Reading data from files
"""

# ╔═╡ 789c489f-51a4-4388-b294-3a7a289702e9
project_directory_files = cd(readdir, pwd());

# ╔═╡ 848fb88d-976c-49b9-9bff-e36f8ca8a03a
md"""
Pick data folder (has to be within the notebooks directory)

$(@bind data_folder Select(project_directory_files))
"""

# ╔═╡ 869346f8-9afd-4c53-8fb8-3c23e181a0af
md"""
Pick data type:

$(@bind data_type Select(["txt", "pickle"]))
"""

# ╔═╡ c419341d-631d-43ab-8fe0-87bf5684a580
if data_type == "txt"
	data_names= cd(readdir, pwd() * "/$(data_folder)/sampling/");
elseif data_type == "pickle"
	data_names = cd(readdir, pwd() * ("/$(data_folder)/pickles/"));
end;

# ╔═╡ 48d73436-d684-4715-b244-14ef8138e6d3
if data_type == "txt"
	md"""
	Pick membrane voltage: $(@bind voltage Select(data_names))
	"""
else
	md"""
	Pick pickle type: $(@bind pickle_sub Select(data_names))
	"""
end

# ╔═╡ f8985d91-0e26-4e73-adb6-a244e8d2264b
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

# ╔═╡ 8e621678-af94-416e-a0df-4c0b4b44ceb5
md"""
Data file: $(@bind data_file Select(data_filenames))
"""

# ╔═╡ 77220080-5ede-4641-9e8b-67066a37af2f
begin
	if data_type == "txt"
	local dt = split(data_file, '.')
	dt[1] = dt[1]*"dwell_timesy"
	md"""
	Dwell times file: $(dwell_times_file = join(dt, '.'))
	"""
	end
end

# ╔═╡ 874cda22-1080-4271-9c69-5335aa8ed041
begin
	if data_type == "txt"
		data_file_path = path_data * data_file
		dwell_times_path = path_dwell_times * dwell_times_file
	else
		data_file_path = path_data * data_file
		dwell_times_path = ""
	end
end

# ╔═╡ bd9f4e7c-0069-4b61-bd71-4e149c1f6aff
md"""
### Standardizing data
"""

# ╔═╡ 47620f98-a1b9-4d39-9911-26c27edfe4bf
Δt::Float32 = 1e-4

# ╔═╡ 2d89d12d-3b1b-4851-8e46-e83c3431f6a4
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
	
	$(@bind d_size NumberField(1000:1000:length(x);default=10000))
	"""
	
end

# ╔═╡ a2bb9e0a-0566-468a-a4c9-81e6cbe00d86
data_size = UInt32(d_size)

# ╔═╡ 02435706-2495-49ff-b313-0021cfab445a
begin
	data = get_specified_datapoints(x, y, Δt, data_size)
	normalized_data = normalize_data(data)
	data["x"] = normalized_data
end

# ╔═╡ 356e6684-8063-49e3-b635-df3d81f9e05e
cumsum(data["dwell times"])

# ╔═╡ 5123c61b-7049-4af9-b2db-2800414e3ad2
md"""
## Implementation of MDL Method
"""

# ╔═╡ fef2cc53-d97a-4695-ac71-2caae03922d5
md"""
### Plotting data
"""

# ╔═╡ 90e2c2e2-652e-492f-beb1-58462466c98c
md"""
## Checking the accuracy of the method
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

# ╔═╡ 512bfa00-d3e9-4eae-8d40-403dd96d17c8
# ╠═╡ disabled = true
#=╠═╡
mean_error(m, Δt, UInt32(225000), true)
  ╠═╡ =#

# ╔═╡ 18a20f62-284b-42ca-bad3-ebef333cfda8
begin
	md"""
	Left range index $(@bind N_left Slider(0:data_size; default=0.0, show_value=true))
	"""
end

# ╔═╡ 2cc62dc0-75d2-46a3-90c2-c3dd08c3a3e2
T_left = trunc(N_left * Δt; digits=4)

# ╔═╡ a61b3892-d50e-46ed-8a04-fcbe1f11e43e
begin
	md"""
	Right range index $(@bind N_right Slider(N_left:data_size-1; default=N_left+500, show_value=true))
	"""
end

# ╔═╡ f558c3bb-a336-455b-8397-8d8e8c6707c7
T_right = trunc(N_right * Δt ;digits=4)

# ╔═╡ e1010167-91af-4670-bef4-3b2824789aec
begin
	md"""
	Minimum segments to work on: $(@bind min_seg Slider(2:300; default=2, show_value=true))
	
	Threshold for `stepstat_mdl` $(@bind threshold Slider(0.00:0.01:1.0; default=1.0, show_value=true))
	
	Number of bins for a histogram $(@bind bins Slider(40:300; default=100, show_value=true))
	"""
end

# ╔═╡ 76268105-9dd6-4e2d-a1f5-f87e6d927a61
n_bins = UInt16(bins)

# ╔═╡ f2ac42ea-ebb0-42c4-9269-dccdcd575e74
m = MDLMethod(min_seg, threshold)

# ╔═╡ 99f08837-ddcd-4eb3-9550-3aff9f257e8c
begin
	method_output = calculate_method(normalized_data, m, Δt)
	mean²error, h_dwell_times, h_dwell_times_approx = calculate_mean_square_error(data, method_output.dwell_times_approx, n_bins)
end

# ╔═╡ c05856d0-e0e0-4ee0-87bb-57f2d566fb63
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

# ╔═╡ 218d6150-eb02-4edb-84fa-0579bd539e5e
begin
	actual_idealized_data = actual_idealize_data(data, what_first_dict, data_file, Δt)
	accuracy = accuracy_of_idealization(actual_idealized_data, method_output.idealized_data)
end

# ╔═╡ f7143895-5b25-43a2-aef5-21807505c128
md"""
Mean squared error $(mean²error)

Accuracy $(accuracy) or $(1 - accuracy)
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

Accuracy $(accuracy) or $(1 - accuracy)
"""

# ╔═╡ 980869db-e712-4344-978e-0205e24cc1d2
method_output

# ╔═╡ 6d9cf785-e38e-46e1-990d-7df46d794e68
begin
	plot_mdl_timestep(data, method_output, T_left, T_right, Δt)
	vline!(method_output.unfiltered_breaks, alpha=0.05)
end

# ╔═╡ d1355f80-ecf1-451b-9225-131b863e2f82
method_output.unfiltered_breaks

# ╔═╡ 1031bc9a-da74-4832-9be1-d97998f09821
md"""
Remarks:
Pickle data has way less activity than experimental data (i.e the channel changes states less frequently) making the method better, with use of longer segments. Experimental data because of its' activity requires segments as small as 2.
"""

# ╔═╡ d0a1be63-d41e-4f3b-9a79-e1c8fc1b2668
plot(IonChannel.histogram_calculator(data["x"]))

# ╔═╡ Cell order:
# ╟─91ef147a-729a-11f0-1157-03caaf19ff7b
# ╠═dbd814ae-a166-4096-a3bc-69a169aa1e5a
# ╟─4cf37e3f-6ac4-4a33-9b45-071c8d4e2347
# ╟─e4f89e4e-6d7e-40c2-9ce0-25924dcb6fee
# ╠═789c489f-51a4-4388-b294-3a7a289702e9
# ╠═848fb88d-976c-49b9-9bff-e36f8ca8a03a
# ╠═869346f8-9afd-4c53-8fb8-3c23e181a0af
# ╠═c419341d-631d-43ab-8fe0-87bf5684a580
# ╠═48d73436-d684-4715-b244-14ef8138e6d3
# ╠═f8985d91-0e26-4e73-adb6-a244e8d2264b
# ╠═8e621678-af94-416e-a0df-4c0b4b44ceb5
# ╠═77220080-5ede-4641-9e8b-67066a37af2f
# ╠═c05856d0-e0e0-4ee0-87bb-57f2d566fb63
# ╠═874cda22-1080-4271-9c69-5335aa8ed041
# ╠═2d89d12d-3b1b-4851-8e46-e83c3431f6a4
# ╠═f66d33ae-2c3e-4661-84fe-8de5a29ae533
# ╟─a2bb9e0a-0566-468a-a4c9-81e6cbe00d86
# ╟─bd9f4e7c-0069-4b61-bd71-4e149c1f6aff
# ╠═47620f98-a1b9-4d39-9911-26c27edfe4bf
# ╠═02435706-2495-49ff-b313-0021cfab445a
# ╠═356e6684-8063-49e3-b635-df3d81f9e05e
# ╟─5123c61b-7049-4af9-b2db-2800414e3ad2
# ╟─76268105-9dd6-4e2d-a1f5-f87e6d927a61
# ╟─f7143895-5b25-43a2-aef5-21807505c128
# ╟─fef2cc53-d97a-4695-ac71-2caae03922d5
# ╠═2cc62dc0-75d2-46a3-90c2-c3dd08c3a3e2
# ╠═f558c3bb-a336-455b-8397-8d8e8c6707c7
# ╠═a93aabed-114e-48a9-a322-f16579fd19f5
# ╠═0e60f2af-655f-4079-915f-656179159d08
# ╠═090447a8-0469-4199-a2b5-f898e44e40c3
# ╟─90e2c2e2-652e-492f-beb1-58462466c98c
# ╠═f2ac42ea-ebb0-42c4-9269-dccdcd575e74
# ╠═99f08837-ddcd-4eb3-9550-3aff9f257e8c
# ╠═e5f04386-4d5b-409e-80cd-11c010a32ffe
# ╠═c56f0f6e-5264-4220-85a4-8f8fae6f5ffb
# ╠═511a0cd0-3e04-46cf-9a36-52dc7082dab2
# ╠═218d6150-eb02-4edb-84fa-0579bd539e5e
# ╠═3b891242-97c6-4bfe-b4e2-bdb72639a688
# ╟─405c80b1-d86f-4139-acb8-57bb04513573
# ╠═512bfa00-d3e9-4eae-8d40-403dd96d17c8
# ╠═980869db-e712-4344-978e-0205e24cc1d2
# ╟─18a20f62-284b-42ca-bad3-ebef333cfda8
# ╟─a61b3892-d50e-46ed-8a04-fcbe1f11e43e
# ╠═e1010167-91af-4670-bef4-3b2824789aec
# ╠═6d9cf785-e38e-46e1-990d-7df46d794e68
# ╠═d1355f80-ecf1-451b-9225-131b863e2f82
# ╟─1031bc9a-da74-4832-9be1-d97998f09821
# ╠═d0a1be63-d41e-4f3b-9a79-e1c8fc1b2668
