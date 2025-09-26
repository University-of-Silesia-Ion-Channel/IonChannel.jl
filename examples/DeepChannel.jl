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

# ╔═╡ e2b4a3b8-5965-430e-bc5a-4c9719c628bf
begin
	include("../src/IonChannel.jl")
	using .IonChannel, PyCall
	import .IonChannel: plot, plot!, title!
	import Pkg
	Pkg.activate(".")
	using PlutoUI
end

# ╔═╡ 0de62389-9084-4025-af37-9bd80abbbb48
md"""
# Machine learning method
"""

# ╔═╡ 168eebad-8721-4bdd-935d-480c8de65a52
md"""
### Reading data from files
"""

# ╔═╡ b8faaa19-5175-42be-9147-bd595131c200
project_directory_files = cd(readdir, pwd())

# ╔═╡ 0e270caf-2f13-41dc-9e1a-aa4c09e6fca1
md"""
Pick data folder (has to be within the notebooks directory)

$(@bind data_folder Select(project_directory_files))
"""

# ╔═╡ 5db88048-45c9-46b0-8df0-5fd976a9d235
md"""
Pick data type:

$(@bind data_type Select(["txt", "pickle"]))
"""

# ╔═╡ f5e97166-deba-458a-9eab-6cb505edb00a
if data_type == "txt"
	data_names= cd(readdir, pwd() * "/$(data_folder)/sampling/");
elseif data_type == "pickle"
	data_names = cd(readdir, pwd() * ("/$(data_folder)/pickles/"));
end;

# ╔═╡ 78e887c8-1624-4935-afaf-31aefceafe39
if data_type == "txt"
	md"""
	Pick membrane voltage: $(@bind voltage Select(data_names))
	"""
else
	md"""
	Pick pickle type: $(@bind pickle_sub Select(data_names))
	"""
end

# ╔═╡ cc7abcec-6015-478f-8aac-28e19eae526f
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

# ╔═╡ 7f0c1415-f82c-41cb-a309-c7c5331c9c2b
md"""
Data file: $(@bind data_file Select(data_filenames))
"""

# ╔═╡ c440bfa0-dfa0-4611-9b50-0b62424ec36d
begin
	if data_type == "txt"
	local dt = split(data_file, '.')
	dt[1] = dt[1]*"dwell_timesy"
	md"""
	Dwell times file: $(dwell_times_file = join(dt, '.'))
	"""
	end
end

# ╔═╡ 7363bda1-58f8-494f-bf74-ed200d459ad2
begin
	if data_type == "txt"
		data_file_path = path_data * data_file
		dwell_times_path = path_dwell_times * dwell_times_file
	else
		data_file_path = path_data * data_file
		dwell_times_path = ""
	end
end

# ╔═╡ 84805655-173a-4155-a798-68ce2b5db049
md"""
### Standardizing data
"""

# ╔═╡ 04b70190-6432-4ff5-9081-4970b4ee7644
Δt::Float32 = 1e-4

# ╔═╡ 8332068f-6869-40ba-98ea-682edad0d856
begin
	x, y = read_data(data_file_path, dwell_times_path)
	if data_type == "pickle"
		y = Δt .* y .* 1000
	end
end

# ╔═╡ ab4517f4-7c65-40b9-926d-bcc142fc390a
begin
	md"""
	Pick how many points to idealize (1000:$(length(x)))
	
	$(@bind data_size NumberField(1000:1000:length(x);default=50000))
	"""
end

# ╔═╡ a671e284-396d-497e-90d6-f00454c028f1
begin
	data = get_specified_datapoints(x, y, Δt, UInt32(data_size))
	normalized_data = normalize_data(data)
	data["x"] = normalized_data
end

# ╔═╡ b50a1e47-c544-4670-a31b-e0eb14d7ebcf
models_path = pwd() * "/models/"

# ╔═╡ 121ae240-56c1-4681-a647-753190041954
models = cd(readdir, models_path)

# ╔═╡ e34612f2-3b54-4a8e-8318-3b71607468d9
md"""
Choose model: $(@bind model_file Select(models))
"""

# ╔═╡ b7c9c99d-3b8f-47d4-99d5-470c8ac8506d
begin
	keras = pyimport("tensorflow.keras")
	model = keras.models.load_model(models_path * model_file)
end

# ╔═╡ ac1caaad-1e8e-47d7-8816-413f4e963db1
prediction = deep_channel_method(data["x"], Δt, DeepChannelMethod(model))

# ╔═╡ c7fac3b8-953b-45dd-8f6a-c2f4740b9db6
data["dwell times"]

# ╔═╡ 0ed9c483-8023-4891-ac65-6d239674b05c
error, h1, h2 = calculate_mean_square_error(data, prediction.dwell_times_approx)

# ╔═╡ e640bc2a-a823-488d-b4de-8561525bdc6f
plot([h1, h2])

# ╔═╡ 5cf2cf5e-c867-4461-bb18-7437ad14d525
begin
	md"""
	Left range index $(@bind N_left Slider(0:data_size; default=0.0, show_value=true))
	"""
end

# ╔═╡ 9c1e5d2c-7745-4e7a-8927-c3baeb9b1bea
begin
	md"""
	Right range index $(@bind N_right Slider(N_left:data_size-1; default=N_left+500, show_value=true))
	"""
end

# ╔═╡ e904b7ce-a027-4174-8740-0bb15554b345
T_left = trunc(N_left * Δt; digits=4)

# ╔═╡ a6b856f0-345b-4905-91f2-caf9edcf5494
T_right = trunc(N_right * Δt ;digits=4)

# ╔═╡ 0f42c978-80a3-425b-91f2-b1b4e53d4ef2
show_approx_on_plot(data, prediction, T_left, T_right, Δt)

# ╔═╡ f59dac24-8e8f-45e1-b5f1-aa295b7bfa1c
# ╠═╡ disabled = true
#=╠═╡
mean_error(DeepChannelMethod(model), Δt, UInt32(225000), true)
  ╠═╡ =#

# ╔═╡ Cell order:
# ╠═e2b4a3b8-5965-430e-bc5a-4c9719c628bf
# ╟─0de62389-9084-4025-af37-9bd80abbbb48
# ╟─168eebad-8721-4bdd-935d-480c8de65a52
# ╠═b8faaa19-5175-42be-9147-bd595131c200
# ╠═0e270caf-2f13-41dc-9e1a-aa4c09e6fca1
# ╠═5db88048-45c9-46b0-8df0-5fd976a9d235
# ╠═f5e97166-deba-458a-9eab-6cb505edb00a
# ╠═78e887c8-1624-4935-afaf-31aefceafe39
# ╠═cc7abcec-6015-478f-8aac-28e19eae526f
# ╠═7f0c1415-f82c-41cb-a309-c7c5331c9c2b
# ╠═c440bfa0-dfa0-4611-9b50-0b62424ec36d
# ╠═7363bda1-58f8-494f-bf74-ed200d459ad2
# ╠═8332068f-6869-40ba-98ea-682edad0d856
# ╠═ab4517f4-7c65-40b9-926d-bcc142fc390a
# ╠═84805655-173a-4155-a798-68ce2b5db049
# ╠═04b70190-6432-4ff5-9081-4970b4ee7644
# ╠═a671e284-396d-497e-90d6-f00454c028f1
# ╠═b50a1e47-c544-4670-a31b-e0eb14d7ebcf
# ╠═121ae240-56c1-4681-a647-753190041954
# ╠═e34612f2-3b54-4a8e-8318-3b71607468d9
# ╠═b7c9c99d-3b8f-47d4-99d5-470c8ac8506d
# ╠═ac1caaad-1e8e-47d7-8816-413f4e963db1
# ╠═c7fac3b8-953b-45dd-8f6a-c2f4740b9db6
# ╠═0ed9c483-8023-4891-ac65-6d239674b05c
# ╠═e640bc2a-a823-488d-b4de-8561525bdc6f
# ╟─5cf2cf5e-c867-4461-bb18-7437ad14d525
# ╠═9c1e5d2c-7745-4e7a-8927-c3baeb9b1bea
# ╟─e904b7ce-a027-4174-8740-0bb15554b345
# ╟─a6b856f0-345b-4905-91f2-caf9edcf5494
# ╠═0f42c978-80a3-425b-91f2-b1b4e53d4ef2
# ╠═f59dac24-8e8f-45e1-b5f1-aa295b7bfa1c
