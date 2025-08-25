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
	using .IonChannel
	import .IonChannel: plot, plot!, title!
	import Pkg
	Pkg.activate(".")
	using PlutoUI
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
begin
	md"""
	Number of bins $(@bind bins Slider(50:300; default=100, show_value=true))
	"""
	bins = UInt16(bins)
end

# ╔═╡ 694c6bad-d8a9-4ff6-b994-eb355daec4f9


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
# ╠═d46741a5-f0ac-4420-9d6a-c779e807916f
# ╠═694c6bad-d8a9-4ff6-b994-eb355daec4f9
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
