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
	using PyCall, Plots
	import Pkg
	Pkg.activate(".")
	using PlutoUI
end

# ╔═╡ 91ef147a-729a-11f0-1157-03caaf19ff7b
md"""
## Loading necessary packages
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

# ╔═╡ 47620f98-a1b9-4d39-9911-26c27edfe4bf
Δt::Float32 = 1e-4

# ╔═╡ e6718caf-64b7-46fe-9712-48603ffb7e74
models_path = pwd() * "/models/"

# ╔═╡ c80dcf87-2c70-47bb-9dd2-9ab2a598c90a
models = cd(readdir, models_path)

# ╔═╡ 69fcd844-de1a-4b3d-83a3-68ea5906b9ab
data_folder_path = pwd() * "/$(data_folder)"

# ╔═╡ b3adf761-5b56-4b63-a3d6-4a4d87c640ea
md"""
Pick file to idealize data $(@bind what_first_path Select(cd(readdir, data_folder_path)))
"""

# ╔═╡ 87a7f110-9e23-443b-be5e-3934fb602e6f
begin
	what_fitst_file_path = pwd() * "/$(data_folder)/$(what_first_path)"
	what_first_dict = Dict(
	    String(split(line,',')[1]) => parse(Int, split(line,',')[2])
	    for line in eachline(what_fitst_file_path)
	)
end

# ╔═╡ b992d12d-6515-422b-a4fa-0763b22626fa
md"""
Choose model: $(@bind model_file Select(models))
"""

# ╔═╡ 1a6b9ad6-6d73-4313-a291-c79d345bdbc5
begin
	keras = pyimport("tensorflow.keras")
	model = keras.models.load_model(models_path * model_file)
end

# ╔═╡ 1572f0e7-f6cc-424f-8292-0ee3d9b2f77d
md"""
Pick membrane voltage: $(@bind voltage Select(voltage_names))
"""

# ╔═╡ b277b495-5381-4591-af60-89ccc8aa80f6
begin
	path_data = pwd() * "/$(data_folder)/sampling/$(voltage)/";
	path_dwell_times = pwd() * "/$(data_folder)/dwell_times/$(voltage)/";
end;

# ╔═╡ 551a682c-d4f7-4337-a735-e103014d957b
begin
	data_filenames = cd(readdir, path_data);
	clean_filenames = filter(
        fname -> occursin(r"^ce\d+\.txt$", fname),
        data_filenames
    )
end

# ╔═╡ bcc7cc74-788a-4e0a-b059-2f00ecd0aba5
data_filenames

# ╔═╡ eee4b6f1-a9f8-45d8-94cd-8bc36c1b5e1f
voltage_names

# ╔═╡ 43f2d875-7e65-4e95-ae13-01c606180bcd
md"""
Data file: $(@bind data_file Select(clean_filenames))
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

# ╔═╡ 4349a657-0eac-40e7-8b84-8696aaa3e690
begin
	md"""
	Pick how many points to idealize (1000:$(length(x)))
	
	$(@bind data_size NumberField(1000:1000:length(x);default=50000))
	"""
	data_size = UInt32(data_size)
end

# ╔═╡ 3f894fd5-2256-45a5-8ef9-039223f32014
begin
	data = get_specified_datapoints(x, y, Δt, data_size)
	normalized_data = normalize_data(data)
	data["x"] = normalized_data
end

# ╔═╡ e6f2282a-fe7c-4467-ae54-6a439a875ac8
function plot_idealization_for_methods(data::Dict{String, Vector{Float32}}, method_outputs::Vector{MethodOutput}, T_left::Float32, T_right::Float32, Δt::Float32)
    @assert T_left <= T_right "N_left must be less or equal to N_right"
    @assert T_right <= (length(data["x"]) - 1) * Δt "T_right exceeds data duration"

    N_left = max(1, Int(round(T_left / Δt)) + 1)
    N_right = min(length(data["x"]), Int(round(T_right / Δt)) + 1)
    N = N_right - N_left + 1

    time = range(T_left, T_left + Δt*(N-1), length=N)

	color_blue = :blue  # a clear, less saturated blue

  # lighter, semi-transparent orange
	
	plots = []
	for method_output in method_outputs
	    if typeof(method_output) <: MikaMethodOutput
			vals = sort(unique(method_output.idealized_data))
			mapped = (method_output.idealized_data .== vals[2])
	        y2 = mapped[N_left:N_right]
			push!(plots, plot(time, y2, color=color_blue, legend=false, title="Mika Method", titlefont=font(6), dpi=500))
	    elseif typeof(method_output) <: MeanDeviationMethodOutput
			y2 = method_output.idealized_data[N_left:N_right]
			push!(plots, plot(time, y2, color=color_blue, legend=false, title="Mean Deviation Method", titlefont=font(6), dpi=500))
		elseif typeof(method_output) <: DeepChannelMethodOutput
			y2 = method_output.idealized_data[N_left:N_right]
			push!(plots, plot(time, y2, color=color_blue, legend=false, title="Deep Channel Method", titlefont=font(6), dpi=500))
		elseif typeof(method_output) <: NaiveMethodOutput
			y2 = method_output.idealized_data[N_left:N_right]
			push!(plots, plot(time, y2, color=color_blue, legend=false, title="Naive Method", titlefont=font(6), dpi=500))
		elseif typeof(method_output) <: MDLMethodOutput
			y2 = method_output.idealized_data[N_left:N_right]
			push!(plots, plot(time, y2, color=color_blue, legend=false, title="MDL Method", titlefont=font(6), dpi=500))
	    end
	end
    y1 = data["x"][N_left:N_right]
	plt1 = plot(time, y1, color=:green, legend=false, title="Ion channel current plot", titlefont=font(6), dpi=500)

	y3 = actual_idealize_data(data, what_first_dict, data_file, Δt)[N_left:N_right]
	plot(plt1, plots..., layout=grid(length(method_outputs) + 1, 1), heights=[0.5, 0.125, 0.125, 0.125, 0.125];  size=(1280, 1280))
	plot!(time, fill(y3, length(method_outputs) + 1), color=:red, alpha=0.7, legend=false, titlefont=font(6), dpi=500, linestyle=:dash)
	

end

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

# ╔═╡ bb51a8e8-0d34-437a-b541-10d721b9ebfe
function plot_mdl(data::Dict{String, Vector{Float64}}, breakpoints::Vector{Float64}, T_left::Float64, T_right::Float64, Δt::Float64)
	@assert T_left <= T_right "T_left must be less or equal to T_right"
	@assert haskey(data, "x") "Data must contain 'x' key with raw signal data"
	@assert haskey(data, "dwell times") "Data must contain 'dwell times' key with dwell segment durations"
	@assert T_left >= 0 "T_left must be non-negative"
	
	T_diff = T_right - T_left
	N = Int(round(T_diff / Δt)) + 1
	time = range(T_left, T_right, length=N)
	
	N_left = max(1, Int(floor(T_left / Δt)) + 1)
	N_right = min(length(data["x"]), Int(ceil(T_right / Δt)) + 1)
	
	plot(time, data["x"][N_left:N_right], dpi=200, label="datapoints", color="green", seriestype=:line)
	
	cumulative_times = cumsum(data["dwell times"])
	indices = findall(t -> t >= T_left && t <= T_right, cumulative_times)
	breakpoints_to_draw = cumulative_times[indices]
	
	approx_indices = findall(t -> t >= T_left && t <= T_right, breakpoints)
	approx_breakpoints_to_draw = breakpoints[approx_indices]
	
	alpha_val = min(0.5, 1 / T_diff)
	vline!(breakpoints_to_draw; alpha=alpha_val, label="Exact", color="red")
	
	vline!(approx_breakpoints_to_draw, label="Approximated", alpha=alpha_val, color="blue")

end


# ╔═╡ 0fa305bd-3398-4abc-a9e8-d347d1997c36
md"""
MDL threshold $(@bind threshold Slider(0.01:0.01:1.0, default=0.8, show_value=true))
"""

# ╔═╡ c88f2283-9725-4996-ab54-745cf73e68b9
md"""
MDL minimum segments $(@bind min_seg Slider(1:300, default=300, show_value=true))
"""

# ╔═╡ 4a6f6699-1b84-4329-846d-08ae252cf762
# methods = [DeepChannelMethod(model), MeanDeviationMethod(0.0, 1.0), MikaMethod(0.0, 100), NaiveMethod(100)]
methods = [DeepChannelMethod(model), MDLMethod(min_seg, threshold, 100), MeanDeviationMethod(0.0, 1.0), MikaMethod(0.0, 100), NaiveMethod(100)]

# ╔═╡ 18d40559-432e-43d9-9027-7efcbc681a7d
begin
	error_outputs = []
	for method in methods
		@info "using $(method)"
		push!(error_outputs, mean_error(method, Δt, UInt32(225000), true))
	end
	error_outputs
end

# ╔═╡ fada3c74-3f60-41fa-b011-5e24e112f1ee
begin
	method_outputs = []
	for method in methods
		push!(method_outputs, calculate_method(data["x"], method, Δt))
	end
	method_outputs = Vector{MethodOutput}(method_outputs)
end

# ╔═╡ dddebe29-e457-41f0-a548-6c31842b9953
plot_idealization_for_methods(data, method_outputs, T_left, T_right, Δt)

# ╔═╡ e2c14ce0-514d-4158-8cd8-f9c417790f00
begin
	accuracy_table = []
	for method in methods
		@info "using $(method)"
		method_output = calculate_method(data["x"], method, Δt)
		actual_idealization = actual_idealize_data(data, what_first_dict, data_file, Δt)
		if typeof(method_output) <: MikaMethodOutput
			vals = sort(unique(method_output.idealized_data))
			mapped = (method_output.idealized_data .== vals[2])
			approx_idealization = Vector{UInt8}(mapped)
		else
			approx_idealization = method_output.idealized_data
		end
		push!(accuracy_table, accuracy_of_idealization(actual_idealization, approx_idealization))
	end
end

# ╔═╡ 1d0925c7-032a-4c2c-a4d8-8e4a25547024
accuracy_table

# ╔═╡ 4bed656f-91f1-470c-b397-dc35d997a086
method = MDLMethod(min_seg, threshold, 100)

# ╔═╡ Cell order:
# ╟─91ef147a-729a-11f0-1157-03caaf19ff7b
# ╠═dbd814ae-a166-4096-a3bc-69a169aa1e5a
# ╠═e4f89e4e-6d7e-40c2-9ce0-25924dcb6fee
# ╠═7044c2e9-cbd0-4b34-a476-3e03624a730c
# ╠═e5c455d4-6313-48ad-a2ba-e72b119429fe
# ╠═1dbb9c49-a92c-4f7a-87cd-72817cbb99cb
# ╟─3ce29b12-01f6-48c6-95e6-233b2109d55f
# ╠═b277b495-5381-4591-af60-89ccc8aa80f6
# ╠═bcc7cc74-788a-4e0a-b059-2f00ecd0aba5
# ╠═1a751b4a-16a8-4270-9e4f-54c200b0a844
# ╠═b9479dd4-a894-4c02-857f-facbfae2ab0e
# ╠═8389df5c-2efe-47f6-8a9c-6e2a230b2b4f
# ╟─471b5827-a5cc-4cf0-9134-62b77750d473
# ╠═4349a657-0eac-40e7-8b84-8696aaa3e690
# ╠═47620f98-a1b9-4d39-9911-26c27edfe4bf
# ╠═3f894fd5-2256-45a5-8ef9-039223f32014
# ╟─e6718caf-64b7-46fe-9712-48603ffb7e74
# ╠═c80dcf87-2c70-47bb-9dd2-9ab2a598c90a
# ╠═69fcd844-de1a-4b3d-83a3-68ea5906b9ab
# ╠═b3adf761-5b56-4b63-a3d6-4a4d87c640ea
# ╠═87a7f110-9e23-443b-be5e-3934fb602e6f
# ╠═b992d12d-6515-422b-a4fa-0763b22626fa
# ╠═1a6b9ad6-6d73-4313-a291-c79d345bdbc5
# ╠═4a6f6699-1b84-4329-846d-08ae252cf762
# ╠═18d40559-432e-43d9-9027-7efcbc681a7d
# ╠═fada3c74-3f60-41fa-b011-5e24e112f1ee
# ╠═e6f2282a-fe7c-4467-ae54-6a439a875ac8
# ╠═dddebe29-e457-41f0-a548-6c31842b9953
# ╟─551a682c-d4f7-4337-a735-e103014d957b
# ╟─1572f0e7-f6cc-424f-8292-0ee3d9b2f77d
# ╟─eee4b6f1-a9f8-45d8-94cd-8bc36c1b5e1f
# ╠═43f2d875-7e65-4e95-ae13-01c606180bcd
# ╟─18a20f62-284b-42ca-bad3-ebef333cfda8
# ╟─a61b3892-d50e-46ed-8a04-fcbe1f11e43e
# ╟─2cc62dc0-75d2-46a3-90c2-c3dd08c3a3e2
# ╟─f558c3bb-a336-455b-8397-8d8e8c6707c7
# ╠═e2c14ce0-514d-4158-8cd8-f9c417790f00
# ╠═1d0925c7-032a-4c2c-a4d8-8e4a25547024
# ╠═4bed656f-91f1-470c-b397-dc35d997a086
# ╠═bb51a8e8-0d34-437a-b541-10d721b9ebfe
# ╠═0fa305bd-3398-4abc-a9e8-d347d1997c36
# ╠═c88f2283-9725-4996-ab54-745cf73e68b9
