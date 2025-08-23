### A Pluto.jl notebook ###
# v0.20.16

using Markdown
using InteractiveUtils

# ╔═╡ 6621c8ee-8008-11f0-3f5b-7dcaa4dac102
begin
	all_data_files = []
	all_dwell_times_files = []
	for voltage in cd(readdir, "./data/sampling")
		data_filenames = cd(readdir, "./data/sampling/$voltage")
		clean_filenames = filter(
	        fname -> occursin(r"^ce\d+\.txt$", fname),
	        data_filenames
	    )
		push!(all_data_files, "./data/sampling/$voltage/" .* clean_filenames)
		dwell_times_filenames = cd(readdir, "./data/dwell_times/$voltage")
		clean_dt_filenames = filter(
	        fname -> occursin(r".*.txt$", fname),
	        dwell_times_filenames
	    )
		push!(all_dwell_times_files, "./data/dwell_times/$voltage/" .* clean_dt_filenames)
	end
	
	idealized_big_data = Dict{String, Vector{Int8}}([])
	local what_first_i = 1
	for voltage_i in 1:length(all_data_files)
		for sample_i in 1:length(all_data_files[voltage_i])
			@info "$(all_data_files[voltage_i][sample_i])"
			x, y = read_data(all_data_files[voltage_i][sample_i], all_dwell_times_files[voltage_i][sample_i])
			file_data = get_specified_datapoints(x, y, Δt)
			filename = split(all_data_files[voltage_i][sample_i], '/')[end]
			@info "$(filename)"
			what_first = what_first_dict[filename]
			idealized_value = what_first
			@info "First value of idalization: $(what_first)"
			idealized_values = actual_idealize_data(file_data, what_first_dict, filename, Δt)
			idealized_big_data[split(all_data_files[voltage_i][sample_i], "/")[end]] = idealized_values
		end
	end
end

# ╔═╡ 24a6f52f-1b2c-4025-bec8-a31ca247b96f
open("idealizations.txt", "w") do file
    for (key, values) in idealized_big_data
        write(file, key * "\n")
        values_str = join(string.(values), ",")
        write(file, values_str * "\n")
    end
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.6"
manifest_format = "2.0"
project_hash = "da39a3ee5e6b4b0d3255bfef95601890afd80709"

[deps]
"""

# ╔═╡ Cell order:
# ╠═6621c8ee-8008-11f0-3f5b-7dcaa4dac102
# ╠═24a6f52f-1b2c-4025-bec8-a31ca247b96f
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
