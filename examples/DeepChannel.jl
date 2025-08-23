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

# ╔═╡ d33efcba-7c16-11f0-3f37-534c0931ce92
begin
	import Pkg
	
	Pkg.activate("../")
	ENV["PYTHON_JL_RUNTIME_PYTHON"] = Sys.which("python3")
	ENV["PYTHON"] = Sys.which("python3")
	Pkg.add(["PyCall", "PlutoUI", "Plots", "OpenSSL_jll"])
	Pkg.build("PyCall")
	using PyCall
	using OpenSSL_jll, PlutoUI, Plots
end

# ╔═╡ e2b4a3b8-5965-430e-bc5a-4c9719c628bf
begin
	include("../src/IonChannel.jl")
	using .IonChannel
	import .IonChannel: plot, plot!, title!
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

# ╔═╡ 5b955299-4e73-455b-b049-4c19c0af2bab
md"""
Pick data folder (has to be within the notebooks directory)

$(@bind data_folder Select(project_directory_files))
"""

# ╔═╡ e989640d-cc49-456b-9033-5d49375e7b9b
voltage_names= cd(readdir, pwd() * "/$(data_folder)/sampling/");

# ╔═╡ 2e6bdf37-5340-4ddd-a4ef-f0252bc3834f
md"""
Pick membrane voltage: $(@bind voltage Select(voltage_names))
"""

# ╔═╡ cf12577c-4dfe-4fd9-a1fa-4091cf1a83df
begin
	path_data = pwd() * "/$(data_folder)/sampling/$(voltage)/";
	path_dwell_times = pwd() * "/$(data_folder)/dwell_times/$(voltage)/";
end;

# ╔═╡ 9ab4502a-4555-4931-8536-0e7d6a495c8c
begin
	data_filenames = cd(readdir, path_data);
	clean_filenames = filter(
        fname -> occursin(r"^ce\d+\.txt$", fname),
        data_filenames
    )
end

# ╔═╡ e2532e7b-2654-4b78-9786-ba83db187ef6
md"""
Data file: $(@bind data_file Select(clean_filenames))
"""

# ╔═╡ 01c4d515-2b61-4e61-9560-ad9bb8f2d2be
begin
	local dt = split(data_file, '.')
	dt[1] = dt[1]*"dwell_timesy"
	md"""
	Dwell times file: $(dwell_times_file = join(dt, '.'))
	"""
end

# ╔═╡ c4009667-2e93-44b2-82e0-81a03d9747df
begin
	data_file_path = path_data * data_file
	dwell_times_path = path_dwell_times * dwell_times_file
end;

# ╔═╡ a16d98ce-e7a5-4f6d-a015-c7130df3e4dc
x, y = read_data(data_file_path, dwell_times_path)

# ╔═╡ ab4517f4-7c65-40b9-926d-bcc142fc390a
md"""
Pick how many points to idealize (1000:$(length(x)))

$(@bind data_size NumberField(1000:1000:length(x);default=50000))
"""

# ╔═╡ 84805655-173a-4155-a798-68ce2b5db049
md"""
### Standardizing data
"""

# ╔═╡ 04b70190-6432-4ff5-9081-4970b4ee7644
Δt = 1e-4

# ╔═╡ a671e284-396d-497e-90d6-f00454c028f1
begin
	data = get_specified_datapoints(x, y, Δt, data_size)
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
	Right range index $(@bind N_right Slider(N_left:data_size; default=N_left+500, show_value=true))
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
# ╠═d33efcba-7c16-11f0-3f37-534c0931ce92
# ╠═e2b4a3b8-5965-430e-bc5a-4c9719c628bf
# ╟─0de62389-9084-4025-af37-9bd80abbbb48
# ╠═168eebad-8721-4bdd-935d-480c8de65a52
# ╠═b8faaa19-5175-42be-9147-bd595131c200
# ╠═5b955299-4e73-455b-b049-4c19c0af2bab
# ╠═e989640d-cc49-456b-9033-5d49375e7b9b
# ╠═2e6bdf37-5340-4ddd-a4ef-f0252bc3834f
# ╠═cf12577c-4dfe-4fd9-a1fa-4091cf1a83df
# ╠═9ab4502a-4555-4931-8536-0e7d6a495c8c
# ╠═e2532e7b-2654-4b78-9786-ba83db187ef6
# ╠═01c4d515-2b61-4e61-9560-ad9bb8f2d2be
# ╠═c4009667-2e93-44b2-82e0-81a03d9747df
# ╠═a16d98ce-e7a5-4f6d-a015-c7130df3e4dc
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
# ╟─9c1e5d2c-7745-4e7a-8927-c3baeb9b1bea
# ╟─e904b7ce-a027-4174-8740-0bb15554b345
# ╟─a6b856f0-345b-4905-91f2-caf9edcf5494
# ╠═0f42c978-80a3-425b-91f2-b1b4e53d4ef2
# ╠═f59dac24-8e8f-45e1-b5f1-aa295b7bfa1c
