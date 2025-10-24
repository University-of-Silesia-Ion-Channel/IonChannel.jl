### A Pluto.jl notebook ###
# v0.20.19

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
	# include("../src/IonChannel.jl")
	# using .IonChannel
	# using PyCall
	import Pkg
	Pkg.activate(".")
	# Pkg.add("Revise")
	using Revise
	# Pkg.add("DataFrames")
	# Pkg.add("CSV")
	# Pkg.add("Plots")
	using PlutoUI, DataFrames, CSV, Plots
	using IonChannel
	Revise.revise()
end

# ╔═╡ 91ef147a-729a-11f0-1157-03caaf19ff7b
md"""
## Loading necessary packages
"""

# ╔═╡ da4ca310-423a-4dec-89e4-c508f143b225
plotly()

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

# ╔═╡ da9e0281-ea18-4f4d-84fd-2022826cdb12
md"""
### Reading data from files
"""

# ╔═╡ 234d65fb-c411-45a4-8cc3-3a32b4df8ec2
project_directory_files = cd(readdir, pwd());

# ╔═╡ bbb40b9b-d7df-40e6-af32-57c03420583e
md"""
Pick data folder (has to be within the notebooks directory)

$(@bind data_folder Select(project_directory_files))
"""

# ╔═╡ 07583db4-1eba-4084-9a65-0b31173f37d9
md"""
Pick data type:

$(@bind data_type Select(["txt", "pickle"]))
"""

# ╔═╡ 960fb302-6252-4b90-bf1f-f5d744c4c2ac
if data_type == "txt"
	data_names= cd(readdir, pwd() * "/$(data_folder)/sampling/");
elseif data_type == "pickle"
	data_names = cd(readdir, pwd() * ("/$(data_folder)/pickles/"));
end;

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
	    String(split(line,',')[1]) => parse(UInt8, split(line,',')[2])
	    for line in eachline(what_fitst_file_path)
	)
end

# ╔═╡ b992d12d-6515-422b-a4fa-0763b22626fa
md"""
Choose model: $(@bind model_file Select(models))
"""

# ╔═╡ 1a6b9ad6-6d73-4313-a291-c79d345bdbc5
begin
	keras = IonChannel.pyimport("tensorflow.keras")
	model = keras.models.load_model(models_path * model_file)
end

# ╔═╡ d3d9a860-61e4-48aa-8bb2-f7adde600f51
# ╠═╡ disabled = true
#=╠═╡
begin
	error_outputs = Dict([])
	for method in methods_pickle
		@info "using $(method)"
		m_table, m_acc, m_error = mean_error_pickle(method, Δt, UInt32(5000), true)
		error_outputs[string(split(string(typeof(method)), '.')[end])] = (m_table, m_acc, m_error)
	end
	error_outputs
end
  ╠═╡ =#

# ╔═╡ 946919c5-db55-4357-a57a-ae205e7f05e4
function vector_dict_to_df(d::Dict{String,Vector{Float32}})
	# Determine maximum length among all vectors
	maxlen = isempty(d) ? 0 : maximum(length.(values(d)))
	# Build a NamedTuple of columns with element type Union{Missing,Float32}
	cols = (; (Symbol(k) => Union{Missing,Float32}[ i <= length(v) ? v[i] : missing
											for i in 1:maxlen ]
				for (k, v) in d)...)
	DataFrame(cols)
end

# ╔═╡ 642fd764-2e35-435e-b553-d6f2afc6265a
# ╠═╡ disabled = true
#=╠═╡
begin
	
	df_errors = Dict()
	df_accuracies = Dict()
	df_final_summary = Dict()
	for (method, method_error) in error_outputs
		@info "$method"
		full_table, mean_accuracies, mean_errors = method_error
		# @info "$(full_table)"
		@info "$(keys(mean_errors))"
		temp_df_error = Dict()
		temp_df_accuracy = Dict()
		temp_df_summary = Dict()
		for pickle_type in keys(mean_errors)
			# save 
			temp_df_error[pickle_type] = vector_dict_to_df(full_table["errors"][pickle_type])
			temp_df_accuracy[pickle_type] = vector_dict_to_df(full_table["accuracies"][pickle_type])
			noise_levels = collect(keys(mean_accuracies[pickle_type]))
			@info "$noise_levels"
			@info "$([mean_errors[pickle_type][n][1] for n in noise_levels])"
			df_summary = DataFrame(
				noise_level = noise_levels,
				mean_error = [mean_errors[pickle_type][n][1] for n in noise_levels],
				mean_error_std = [mean_errors[pickle_type][n][2] for n in noise_levels],
				mean_accuracy = [get(mean_accuracies[pickle_type], n, NaN32)[1] for n in noise_levels],
				mean_accuracy_std = [get(mean_accuracies[pickle_type], n, NaN32)[2] for n in noise_levels]
			)
			@info "$df_summary"
			temp_df_summary[pickle_type] = df_summary
		end
		df_errors[method] = temp_df_error
		df_accuracies[method] = temp_df_accuracy

		df_final_summary[method] = temp_df_summary
		# full table
		# errors/accuracies -> pickle_type(m20/p20) -> VL/L/M/H/VH -> Vector
	
		# mean_accuracies
		# pickle_type(m20/p20) -> VL/L/M/H/VH -> Vector
	
		# mean_errors
		# pickle_type(m20/p20) -> VL/L/M/H/VH -> Vector
		
	end
end
  ╠═╡ =#

# ╔═╡ efe8289b-4b3d-4a8d-9f82-cdd8e8a8dd34
# ╠═╡ disabled = true
#=╠═╡
for dir in keys(df_errors)
	# @info "$dir"
	for pickle_type in keys(df_errors[dir])
		# @info "$pickle_type"
		try
			mkdir("../outp/$(dir)")
		catch e
			println("Directory already exists")
		end
		try
			mkdir("../outp/$(dir)/$(pickle_type)")
		catch e
			println("Directory already exists")
		end
		CSV.write("../outp/$(dir)/$(pickle_type)/errors.csv", df_errors[dir][pickle_type])
		CSV.write("../outp/$(dir)/$(pickle_type)/accuracies.csv", df_accuracies[dir][pickle_type])
		CSV.write("../outp/$(dir)/$(pickle_type)/mean_errors_and_accuracies.csv", df_final_summary[dir][pickle_type])
	end
end
  ╠═╡ =#

# ╔═╡ fcdf4e9d-d2bc-45d9-989c-d0f38b8f5774
# ╠═╡ disabled = true
#=╠═╡
for filename in keys(error_outputs)
	# mkdir("../outp/$(filename)")
	CSV.write("../outp/$(filename)/mean_errors.csv", error_outputs[filename][1])
	CSV.write("../outp/$(filename)/mean_accuracies.csv", error_outputs[filename][2])
	CSV.write("../outp/$(filename)/mean_accuracies_voltages.csv", error_outputs[filename][3])
end
  ╠═╡ =#

# ╔═╡ 4a6f6699-1b84-4329-846d-08ae252cf762
begin
	# methods = [DeepChannelMethod(model), MeanDeviationMethod(0.0, 1.0), MikaMethod(0.0, 100), NaiveMethod(100)]
	methods_txt = [MDLMethod(UInt16(2), Float32(0.51)), MikaMethod(), DeepChannelMethod(model), MeanDeviationMethod(0.0), NaiveMethod()]
	methods_pickle = deepcopy(methods_txt)
	methods_pickle[1] = MDLMethod(UInt16(100), Float32(0.15))
end

# ╔═╡ 18d40559-432e-43d9-9027-7efcbc681a7d
begin
	error_outputs_txt = Dict([])
	for method in methods_txt
		@info "using $(method)"
		m_table, m_acc, m_error = mean_error_txt(method, Δt, UInt32(225000), true)
		error_outputs_txt[string(split(string(typeof(method)), '.')[end])] = (m_table, m_acc, m_error)
	end
	error_outputs_txt
end

# ╔═╡ 4e86a3c1-39ac-4075-98a5-2f23af247b0b
begin
	
	df_errors = Dict()
	df_accuracies = Dict()
	df_final_summary = Dict()
	for (method, method_error) in error_outputs_txt
		@info "$method"
		full_table, mean_accuracies, mean_errors = method_error		
		# save 
		temp_df_error = vector_dict_to_df(full_table["errors"])
		temp_df_accuracy = vector_dict_to_df(full_table["accuracies"])
		noise_levels = collect(keys(mean_accuracies))
		@info "$noise_levels"
		@info "$([mean_errors[n][1] for n in noise_levels])"
		df_summary = DataFrame(
			noise_level = noise_levels,
			mean_error = [mean_errors[n][1] for n in noise_levels],
			mean_error_std = [mean_errors[n][2] for n in noise_levels],
			mean_accuracy = [get(mean_accuracies, n, NaN32)[1] for n in noise_levels],
			mean_accuracy_std = [get(mean_accuracies, n, NaN32)[2] for n in noise_levels]
		)
		@info "$df_summary"
		temp_df_summary = df_summary
		
		df_errors[method] = temp_df_error
		df_accuracies[method] = temp_df_accuracy

		df_final_summary[method] = temp_df_summary
		# full table
		# errors/accuracies -> pickle_type(m20/p20) -> VL/L/M/H/VH -> Vector
	
		# mean_accuracies
		# pickle_type(m20/p20) -> VL/L/M/H/VH -> Vector
	
		# mean_errors
		# pickle_type(m20/p20) -> VL/L/M/H/VH -> Vector
		
	end
end

# ╔═╡ 45496c71-959e-451a-aa2f-b91263535cd9
df_accuracies

# ╔═╡ 178e295d-fb7a-4e09-bc81-eccde9e8c386
keys(df_errors)

# ╔═╡ 6159272b-fd43-4dc9-bb2f-bceba7329003
df_final_summary["MikaMethod"]

# ╔═╡ 8b5726e5-1de9-44e1-bac7-68af7b899b24
begin
	# dir = "experimental"
	# @info "$dir"
	direct = "../outp/experimental"
	try
		mkdir(direct)
	catch e
		println("Directory already exists")
	end
	for dir in keys(df_errors)
		# @info "$pickle_type"
		try
			mkdir(direct * "/$(dir)")
		catch e
			println("Directory already exists")
		end
		
		CSV.write(direct * "/$(dir)/errors.csv", df_errors[dir])
		CSV.write(direct * "/$(dir)/accuracies.csv", df_accuracies[dir])
		CSV.write(direct * "/$(dir)/mean_errors_and_accuracies.csv", df_final_summary[dir])
	end
end

# ╔═╡ 8f012c3e-d5f1-4267-a686-ee458a833a2a
md"""
### Checking signal to noise ratio of a given file
"""

# ╔═╡ 9f241a6f-27f2-44b5-887d-81d469294550
if data_type == "txt"
	md"""
	Pick membrane voltage: $(@bind voltage Select(data_names))
	"""
else
	md"""
	Pick pickle type: $(@bind pickle_sub Select(data_names))
	"""
end

# ╔═╡ 1599b194-dfb4-4f3e-bc4b-c62857f74a54
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

# ╔═╡ e336e5e0-94a9-4093-8969-276ae308988f
md"""
Data file: $(@bind data_file Select(data_filenames))
"""

# ╔═╡ 61740328-3bbe-4a40-ad22-5dc1b9e1eb76
begin
	if data_type == "txt"
		local dt = split(data_file, '.')
		dt[1] = dt[1]*"dwell_timesy"
		md"""
		Dwell times file: $(dwell_times_file = join(dt, '.'))
		"""
	end
end

# ╔═╡ c85c05ff-21e9-44de-b908-72a8864abec0
begin
	if data_type == "txt"
		data_file_path = path_data * data_file
		dwell_times_path = path_dwell_times * dwell_times_file
	else
		data_file_path = path_data * data_file
		dwell_times_path = ""
	end
end

# ╔═╡ 8d1de076-505d-4b09-abd4-579d9f3d91ba
data_file_path

# ╔═╡ edd4a010-6986-480b-9ec0-4d7fba52ac93
begin
	x, y = read_data(data_file_path, dwell_times_path)
	if data_type == "pickle"
		y = Δt .* y .* 1000
	end
end

# ╔═╡ 4349a657-0eac-40e7-8b84-8696aaa3e690
begin
	md"""
	Pick how many points to idealize (1000:$(length(x)))
	
	$(@bind data_size NumberField(1000:1000:length(x);default=50000))
	"""
end

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

# ╔═╡ 3f894fd5-2256-45a5-8ef9-039223f32014
begin
	data = get_specified_datapoints(x, y, Δt, UInt32(data_size))
	normalized_data = normalize_data(data)
	data["x"] = normalized_data
end

# ╔═╡ fada3c74-3f60-41fa-b011-5e24e112f1ee
begin
	method_outputs = []
	for method in (data_type == "pickle" ? methods_pickle : methods_txt)
		push!(method_outputs, calculate_method(data["x"], method, Δt))
	end
	method_outputs = Vector{MethodOutput}(method_outputs)
end

# ╔═╡ e6f2282a-fe7c-4467-ae54-6a439a875ac8
function plot_idealization_for_methods(data::Dict{String, Vector{Float32}}, method_outputs::Vector{MethodOutput}, T_left::Float32, T_right::Float32, Δt::Float32)
    @assert T_left <= T_right "N_left must be less or equal to N_right"
    @assert T_right <= (length(data["x"]) - 1) * Δt "T_right exceeds data duration"

    N_left = max(1, Int(round(T_left / Δt)) + 1)
    N_right = min(length(data["x"]), Int(round(T_right / Δt)) + 1)
    N = N_right - N_left + 1

    time = range(T_left, T_left + Δt*(N-1), length=N)

	color_blue = RGB([51,24,252] ./ 255.0 ...)  # a clear, less saturated blue
	color_dotted = RGB([244,53,63] ./ 255.0 ...)
  # lighter, semi-transparent orange

	c_dpi = 1000
	t_font = font(7)
	
	plots = []
	for method_output in method_outputs
	    if typeof(method_output) <: MikaMethodOutput
			vals = sort(unique(method_output.idealized_data))
			mapped = (method_output.idealized_data .== vals[2])
	        y2 = mapped[N_left:N_right]
			push!(plots, plot(time, y2, color=color_blue, legend=false, title="Mika Method", titlefont=t_font, dpi=c_dpi; yticks = 0:1))
	    elseif typeof(method_output) <: MeanDeviationMethodOutput
			y2 = method_output.idealized_data[N_left:N_right]
			push!(plots, plot(time, y2, color=color_blue, legend=false, title="Mean Deviation Method", titlefont=t_font, dpi=c_dpi; yticks = 0:1))
		elseif typeof(method_output) <: DeepChannelMethodOutput
			y2 = method_output.idealized_data[N_left:N_right]
			push!(plots, plot(time, y2, color=color_blue, legend=false, title="Deep Channel Method", titlefont=t_font, dpi=c_dpi; yticks = 0:1))
		elseif typeof(method_output) <: NaiveMethodOutput
			y2 = method_output.idealized_data[N_left:N_right]
			push!(plots, plot(time, y2, color=color_blue, legend=false, title="Naive Method", titlefont=t_font, dpi=c_dpi; yticks = 0:1))
			xlabel!("time [s]")
		elseif typeof(method_output) <: MDLMethodOutput
			y2 = method_output.idealized_data[N_left:N_right]
			push!(plots, plot(time, y2, color=color_blue, legend=false, title="MDL Method", titlefont=t_font, dpi=c_dpi; yticks = 0:1))
	    end
	end
    y1 = data["x"][N_left:N_right]
	plt1 = plot(time, y1, color=:green, legend=false, title="Ion channel current plot", titlefont=t_font, dpi=c_dpi)
	ylabel!("current [μA]")

	y3 = actual_idealize_data(data, what_first_dict, data_file, Δt)[N_left:N_right]
	plot(plt1, plots..., layout=grid(length(method_outputs) + 1, 1, heights=[0.4, 0.12, 0.12, 0.12, 0.12, 0.12]);  size=(680, 800))
	plot!(time, fill(y3, length(method_outputs) + 1), color=color_dotted, alpha=0.9, legend=false, titlefont=t_font, dpi=c_dpi, linestyle=:dot, lw=1.2)
	

end

# ╔═╡ dddebe29-e457-41f0-a548-6c31842b9953
plot_idealization_for_methods(data, method_outputs, T_left, T_right, Δt)

# ╔═╡ e2c14ce0-514d-4158-8cd8-f9c417790f00
begin
	accuracy_table = []
	
	for method in (data_type == "pickle" ? methods_pickle : methods_txt)
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

# ╔═╡ 207e7cb2-1f92-4907-9461-5eca67976314
begin
	global D = 0.0f0
	try
		D = parse.(Float32, join(split(split(data_file, "D")[end], ".")[1:2], "."))
	catch
		D = 1.0f0
	end
end

# ╔═╡ 3b62b303-920d-44ae-b940-1619a13ed286
begin
	dt = IonChannel.fit(IonChannel.UnitRangeTransform, data["x"])
	scaled_data = IonChannel.StatsBase.transform(dt, data["x"])
end

# ╔═╡ 10bef431-4c94-4f18-8bec-c7414f90ddf4
begin
	foo = calculate_probability_histogram(histogram_calculator(scaled_data))
	o = analyze_histogram_peaks(scaled_data)
	bar(o.edges, o.weights)
	vline!([o.edges[o.left_peak_index]], label="left peak", lw=2)
	vline!([o.edges[o.right_peak_index]], label="right peak", lw=2)
	vline!([o.edges[o.pmin_index]], label="minimum", lw=3)
	vline!([o.edges[o.midpoint]], label="midpoint")
end

# ╔═╡ c07ac23d-12e0-4a51-863c-fb0d9a6ae32a
# ╠═╡ disabled = true
#=╠═╡
begin
	struct NoiseAnalysis
		distance_between_peaks::Float32
		distance_between_middle::Float32
		min_height::Float32
	end
	
	function analyze_noise_levels(data_folder)
		what_first_file_path, data_paths, dwell_times_paths = read_all_file_paths(data_folder)
		all_paths = create_paths_dictionary(data_paths, dwell_times_paths)
		noise_analysis = Dict{String, Dict{String, Dict{String, NoiseAnalysis}}}()
		for (data_type, paths) in all_paths
			noise_analysis[data_type] = Dict{String, Dict{String, NoiseAnalysis}}()
			for d_type in keys(paths["data paths"])
				noise_analysis[data_type][d_type] = Dict{String, NoiseAnalysis}()
				N = length(paths["data paths"][d_type])
				# @info "Processing $d_type; $N files"
				for i in 1:N
					file = split(paths["data paths"][d_type][i], "/")[end]
					# @info "Processing file: $(file)"
					x = Float32[]
					y = Float32[]
					if data_type == "txt"
						continue
						# x, y = read_data(paths["data paths"][d_type][i], paths["dwell times paths"][d_type][i])
					else
						x, y = read_data(paths["data paths"][d_type][i])
					end
					data = get_specified_datapoints(x, y, Δt, UInt32(50000))
					dt = IonChannel.fit(IonChannel.UnitRangeTransform, data["x"])
					scaled_data = IonChannel.StatsBase.transform(dt, data["x"])
					# data["x"] = scaled_data
					
					# σ = IonChannel.std(scaled_data)
					# @info "Standard deviation: $(σ)"
					# output_stds[data_type][d_type][String(file)] = σ
	
					# find two peaks on a histogram and find distance between them
					# the further apart the peaks the better signal to noise ratio - maybe
	
					histogram_analysis::HistPeakAnalysis = analyze_histogram_peaks(scaled_data)
					edges = histogram_analysis.edges
					weights = histogram_analysis.weights
					distance_between_peaks::Float32 = abs(edges[histogram_analysis.right_peak_index] - edges[histogram_analysis.left_peak_index])
					distance_between_mid = 0.0f0
					if weights[histogram_analysis.left_peak_index] > weights[histogram_analysis.right_peak_index]
						distance_between_mid = abs(edges[histogram_analysis.left_peak_index] - 0.5f0)
					else
						distance_between_mid = abs(edges[histogram_analysis.right_peak_index] - 0.5f0)
					end
					noise_analysis[data_type][d_type][String(file)] = NoiseAnalysis(distance_between_peaks, distance_between_mid, histogram_analysis.pmin)			
					
					# file_dict[String(file)] = σ
				end
				# data_type_dict[d_type] = f
			end
		end
		noise_analysis
	end
end
  ╠═╡ =#

# ╔═╡ 39792c4e-4a3b-48ea-8442-f8b09e01c00b
#=╠═╡
noise_analysis = analyze_noise_levels(data_folder)
  ╠═╡ =#

# ╔═╡ 7f5fe22e-049f-4520-b6c6-91fe6c768f89
#=╠═╡
begin
	pickle_noise_analysis = noise_analysis["pickle"]
	pickle_noises_dict = Dict{String, Vector{Tuple{Float32, NoiseAnalysis}}}()
	for (type, noises) in pickle_noise_analysis
		temp_vec = Tuple{Float32, NoiseAnalysis}[]
		for (file_name, noise) in noises
			push!(temp_vec, (get_noise_level(file_name), noise))
		end
		sort!(temp_vec, by = x -> x[1])
		pickle_noises_dict[type] = temp_vec
	end
	pickle_noises_dict
end
  ╠═╡ =#

# ╔═╡ e4a0b0a8-922a-4d0d-bf04-c6aab14d1a39


# ╔═╡ e4c7bfd1-efea-44ec-aaae-57419db5bc5e
function plot_noise_levels(pickle_noises_dict, field::Symbol)
	vec_to_plot_p20 = [(D, getfield(noise, field)) for (D, noise) in pickle_noises_dict["p20"]]
	vec_to_plot_m20 = [(D, getfield(noise, field)) for (D, noise) in pickle_noises_dict["m20"]]
	plot(vec_to_plot_p20, label="p20")
	plot!(vec_to_plot_m20, label="m20")
	plot!(; xticks=0:5:100)
	# plot!(zipped_noises_for_plot[2:2:end])
end

# ╔═╡ f32f72fd-9cc3-4e5c-a4d2-846fcc436e4e
#=╠═╡
plot_noise_levels(pickle_noises_dict, field1)
  ╠═╡ =#

# ╔═╡ ea402fe3-c406-4a1d-923d-863b92065fcc
#=╠═╡
plot_noise_levels(pickle_noises_dict, field2)
  ╠═╡ =#

# ╔═╡ 06d9e137-dddf-47e4-ae3d-a94923b733c8
#=╠═╡
plot_noise_levels(pickle_noises_dict, field3)
  ╠═╡ =#

# ╔═╡ 961ef59a-ed64-4491-88f1-4dcec50bc100
begin
	field1 = :distance_between_middle
	field2 = :distance_between_peaks
	field3 = :min_height
	function get_ranges_of_noise_levels(pickle_noises_dict, field::Symbol)
		min_distances_for_d_levels = Dict()
		for (p_type, values) in pickle_noises_dict
			min_distances_for_d_levels[p_type] = Dict("VL" => [1.0f0, 0.0f0], "L" => [1.0f0, 0.0f0], "M" => [1.0f0, 0.0f0], "H" => [1.0f0, 0.0f0], "VH" => [1.0f0, 0.0f0])
			for (noise_level, n_analysis) in values
				if noise_level <= 1.0f0
					min_distances_for_d_levels[p_type]["VL"][1] = getfield(n_analysis, field) < min_distances_for_d_levels[p_type]["VL"][1] ? getfield(n_analysis, field) : min_distances_for_d_levels[p_type]["VL"][1]
					
					min_distances_for_d_levels[p_type]["VL"][2] = getfield(n_analysis, field) > min_distances_for_d_levels[p_type]["VL"][2] ? getfield(n_analysis, field) : min_distances_for_d_levels[p_type]["VL"][2]
				elseif noise_level <= 2.0f0
					min_distances_for_d_levels[p_type]["L"][1] = getfield(n_analysis, field) < min_distances_for_d_levels[p_type]["L"][1] ? getfield(n_analysis, field) : min_distances_for_d_levels[p_type]["L"][1]
					
					min_distances_for_d_levels[p_type]["L"][2] = getfield(n_analysis, field) > min_distances_for_d_levels[p_type]["L"][2] ? getfield(n_analysis, field) : min_distances_for_d_levels[p_type]["L"][2]
				elseif noise_level <= 5.0f0
					min_distances_for_d_levels[p_type]["M"][1] = getfield(n_analysis, field) < min_distances_for_d_levels[p_type]["M"][1] ? getfield(n_analysis, field) : min_distances_for_d_levels[p_type]["M"][1]
					
					min_distances_for_d_levels[p_type]["M"][2] = getfield(n_analysis, field) > min_distances_for_d_levels[p_type]["M"][2] ? getfield(n_analysis, field) : min_distances_for_d_levels[p_type]["M"][2]
				elseif noise_level <= 10.0f0
					min_distances_for_d_levels[p_type]["H"][1] = getfield(n_analysis, field) < min_distances_for_d_levels[p_type]["H"][1] ? getfield(n_analysis, field) : min_distances_for_d_levels[p_type]["H"][1]
					
					min_distances_for_d_levels[p_type]["H"][2] = getfield(n_analysis, field) > min_distances_for_d_levels[p_type]["H"][2] ? getfield(n_analysis, field) : min_distances_for_d_levels[p_type]["H"][2]
				else
					min_distances_for_d_levels[p_type]["VH"][1] = getfield(n_analysis, field) < min_distances_for_d_levels[p_type]["VH"][1] ? getfield(n_analysis, field) : min_distances_for_d_levels[p_type]["VH"][1]
					
					min_distances_for_d_levels[p_type]["VH"][2] = getfield(n_analysis, field) > min_distances_for_d_levels[p_type]["VH"][2] ? getfield(n_analysis, field) : min_distances_for_d_levels[p_type]["VH"][2]
				end	
			end
		end
		min_distances_for_d_levels
	end
end

# ╔═╡ 71437c8b-50e7-4461-a105-b1b74fed68c9
#=╠═╡
distances_for_d_levels = get_ranges_of_noise_levels(pickle_noises_dict, field3)
  ╠═╡ =#

# ╔═╡ d6aea140-fd7c-473d-8b59-2989c09e4023
#=╠═╡
begin
	plt1 = plot(;title="p20")
	for (noise_level, value) in distances_for_d_levels["p20"]
		vspan!(value, label=noise_level, alpha=0.5)
	end
	plt1
end
  ╠═╡ =#

# ╔═╡ faeddaa0-63b7-4f45-81d1-e96ee03056a4
#=╠═╡
begin
	plt2 = plot(; title="m20")
	for (noise_level, value) in distances_for_d_levels["m20"]
		vspan!(value, label=noise_level, alpha=0.5)
	end
	plt2
end
  ╠═╡ =#

# ╔═╡ Cell order:
# ╟─91ef147a-729a-11f0-1157-03caaf19ff7b
# ╠═dbd814ae-a166-4096-a3bc-69a169aa1e5a
# ╠═da4ca310-423a-4dec-89e4-c508f143b225
# ╠═471b5827-a5cc-4cf0-9134-62b77750d473
# ╟─da9e0281-ea18-4f4d-84fd-2022826cdb12
# ╠═234d65fb-c411-45a4-8cc3-3a32b4df8ec2
# ╟─bbb40b9b-d7df-40e6-af32-57c03420583e
# ╟─07583db4-1eba-4084-9a65-0b31173f37d9
# ╠═960fb302-6252-4b90-bf1f-f5d744c4c2ac
# ╠═1599b194-dfb4-4f3e-bc4b-c62857f74a54
# ╠═61740328-3bbe-4a40-ad22-5dc1b9e1eb76
# ╠═c85c05ff-21e9-44de-b908-72a8864abec0
# ╠═8d1de076-505d-4b09-abd4-579d9f3d91ba
# ╠═edd4a010-6986-480b-9ec0-4d7fba52ac93
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
# ╠═18d40559-432e-43d9-9027-7efcbc681a7d
# ╠═d3d9a860-61e4-48aa-8bb2-f7adde600f51
# ╠═946919c5-db55-4357-a57a-ae205e7f05e4
# ╠═4e86a3c1-39ac-4075-98a5-2f23af247b0b
# ╠═642fd764-2e35-435e-b553-d6f2afc6265a
# ╠═45496c71-959e-451a-aa2f-b91263535cd9
# ╠═178e295d-fb7a-4e09-bc81-eccde9e8c386
# ╠═6159272b-fd43-4dc9-bb2f-bceba7329003
# ╠═8b5726e5-1de9-44e1-bac7-68af7b899b24
# ╠═efe8289b-4b3d-4a8d-9f82-cdd8e8a8dd34
# ╠═fcdf4e9d-d2bc-45d9-989c-d0f38b8f5774
# ╠═4a6f6699-1b84-4329-846d-08ae252cf762
# ╠═fada3c74-3f60-41fa-b011-5e24e112f1ee
# ╠═e6f2282a-fe7c-4467-ae54-6a439a875ac8
# ╠═dddebe29-e457-41f0-a548-6c31842b9953
# ╟─18a20f62-284b-42ca-bad3-ebef333cfda8
# ╠═a61b3892-d50e-46ed-8a04-fcbe1f11e43e
# ╟─2cc62dc0-75d2-46a3-90c2-c3dd08c3a3e2
# ╟─f558c3bb-a336-455b-8397-8d8e8c6707c7
# ╠═e2c14ce0-514d-4158-8cd8-f9c417790f00
# ╠═1d0925c7-032a-4c2c-a4d8-8e4a25547024
# ╟─8f012c3e-d5f1-4267-a686-ee458a833a2a
# ╟─9f241a6f-27f2-44b5-887d-81d469294550
# ╟─e336e5e0-94a9-4093-8969-276ae308988f
# ╟─207e7cb2-1f92-4907-9461-5eca67976314
# ╠═3b62b303-920d-44ae-b940-1619a13ed286
# ╠═10bef431-4c94-4f18-8bec-c7414f90ddf4
# ╠═c07ac23d-12e0-4a51-863c-fb0d9a6ae32a
# ╠═39792c4e-4a3b-48ea-8442-f8b09e01c00b
# ╠═7f5fe22e-049f-4520-b6c6-91fe6c768f89
# ╠═e4a0b0a8-922a-4d0d-bf04-c6aab14d1a39
# ╠═e4c7bfd1-efea-44ec-aaae-57419db5bc5e
# ╠═f32f72fd-9cc3-4e5c-a4d2-846fcc436e4e
# ╠═ea402fe3-c406-4a1d-923d-863b92065fcc
# ╠═06d9e137-dddf-47e4-ae3d-a94923b733c8
# ╠═961ef59a-ed64-4491-88f1-4dcec50bc100
# ╠═71437c8b-50e7-4461-a105-b1b74fed68c9
# ╠═d6aea140-fd7c-473d-8b59-2989c09e4023
# ╠═faeddaa0-63b7-4f45-81d1-e96ee03056a4
