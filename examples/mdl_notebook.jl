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
	
	$(@bind d_size NumberField(1000:1000:length(x);default=50000))
	"""
	
end

# ╔═╡ a2bb9e0a-0566-468a-a4c9-81e6cbe00d86
data_size = UInt32(d_size)

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
	    String(split(line,',')[1]) => parse(Int, split(line,',')[2])
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
	Right range index $(@bind N_right Slider(N_left:data_size; default=N_left+500, show_value=true))
	"""
end

# ╔═╡ f558c3bb-a336-455b-8397-8d8e8c6707c7
T_right = trunc(N_right * Δt ;digits=4)

# ╔═╡ e1010167-91af-4670-bef4-3b2824789aec
begin
	md"""
	Minimum segments to work on: $(@bind min_seg Slider(2:300; default=10, show_value=true))
	
	Threshold for `stepstat_mdl` $(@bind threshold Slider(0.00:0.01:1.0; default=0.8, show_value=true))
	
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
plot_mdl_timestep(data, method_output, T_left, T_right, Δt)

# ╔═╡ Cell order:
# ╟─91ef147a-729a-11f0-1157-03caaf19ff7b
# ╠═dbd814ae-a166-4096-a3bc-69a169aa1e5a
# ╟─4cf37e3f-6ac4-4a33-9b45-071c8d4e2347
# ╟─e4f89e4e-6d7e-40c2-9ce0-25924dcb6fee
# ╠═7044c2e9-cbd0-4b34-a476-3e03624a730c
# ╟─e5c455d4-6313-48ad-a2ba-e72b119429fe
# ╠═1dbb9c49-a92c-4f7a-87cd-72817cbb99cb
# ╟─3ce29b12-01f6-48c6-95e6-233b2109d55f
# ╟─1572f0e7-f6cc-424f-8292-0ee3d9b2f77d
# ╟─b277b495-5381-4591-af60-89ccc8aa80f6
# ╟─b22ddb25-23a9-4a3c-80a7-5b8984a95c1d
# ╟─43f2d875-7e65-4e95-ae13-01c606180bcd
# ╟─1a751b4a-16a8-4270-9e4f-54c200b0a844
# ╟─b9479dd4-a894-4c02-857f-facbfae2ab0e
# ╟─8389df5c-2efe-47f6-8a9c-6e2a230b2b4f
# ╟─471b5827-a5cc-4cf0-9134-62b77750d473
# ╟─f66d33ae-2c3e-4661-84fe-8de5a29ae533
# ╟─a2bb9e0a-0566-468a-a4c9-81e6cbe00d86
# ╟─bd9f4e7c-0069-4b61-bd71-4e149c1f6aff
# ╠═47620f98-a1b9-4d39-9911-26c27edfe4bf
# ╟─02435706-2495-49ff-b313-0021cfab445a
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
# ╟─e1010167-91af-4670-bef4-3b2824789aec
# ╠═6d9cf785-e38e-46e1-990d-7df46d794e68
