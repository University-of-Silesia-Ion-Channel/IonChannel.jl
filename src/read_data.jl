using StatsBase
using Normalization

"""
read_data(data_file_path::String, dwell_times_path::String) -> Tuple{Vector{Float32}, Vector{Float32}}

Read numerical data from two text files and return them as vectors of `Float32`.

# Arguments
- `data_file_path::String`: Path to a text file containing numerical values (one per line)
representing the primary data set.
- `dwell_times_path::String`: Path to a text file containing numerical values (one per line)
corresponding to dwell times.

# Returns
A tuple `(x, y)` where:
- `x::Vector{Float32}`: Values read from the first file.
- `y::Vector{Float32}`: Values read from the second file.

# Example
```
x, y = read_data("data.txt", "dwell_times.txt")
```

Both files must contain one floating-point number per line, with optional whitespace.
"""
function read_data(data_file_path::String, dwell_times_path::String) :: 
Tuple{Vector{Float32}, Vector{Float32}}
    # reading data to RAM
    x = open(data_file_path) do f; parse.(Float32, strip.(readlines(f))); end
    y = open(dwell_times_path) do f; parse.(Float32, strip.(readlines(f))); end
    x, y
end

"""
read_all_file_paths(data_folder::String) -> Tuple{String, Vector{String}, Vector{String}}

Collect and return file paths for experimental data files and their corresponding dwell time files
from a structured folder hierarchy.

# Arguments
- `data_folder::String`: The path (relative or absolute) to the main data folder.  
The expected directory structure inside `data_folder` is:

data_folder/
    sampling/
        <voltage1>/
            file1.txt
            file2.txt
            ...
        <voltage2>/
            ...
    dwell_times/
        <voltage1>/
            file1_dwell_timesy.txt
            ...
        <voltage2>/
            ...

# Returns
A tuple:
1. `what_first_path::String` — Path to the first special file found at the top level of `data_folder`
(taken from the 3rd entry in its directory listing).
2. `data_file_paths::Vector{String}` — Full paths to all sampled data files across all voltages.
3. `dwell_times_file_paths::Vector{String}` — Full paths to the corresponding dwell time files,
in the same order as `data_file_paths`.

# File matching logic
- Data files are taken from `/sampling/<voltage>/`, selecting every second file starting at index 2.
- Dwell time files are taken from `/dwell_times/<voltage>/`, selecting every second file starting at index 1.
- The dwell time file names are derived from the data file names by replacing the base name
with `<basename>dwell_timesy` and preserving the original extension.

# Example
```
what_first, data_paths, dwell_paths = read_all_file_paths("experiment_data")
println("First reference file: ", what_first)
println("Number of data files: ", length(data_paths))
println("Number of dwell time files: ", length(dwell_paths))
```
# Notes
- Assumes a specific directory and file naming convention.
- File order consistency is crucial for correctly matching data with dwell times.
"""
function read_all_file_paths(data_folder::String) :: Tuple{String, Vector{String}, Vector{String}}
    voltage_names= cd(readdir, pwd() * "/$(data_folder)/sampling/")
    data_file_paths = []
    dwell_times_file_paths = []
    for voltage in voltage_names
        path_data = pwd() * "/$(data_folder)/sampling/$(voltage)/"
        path_dwell_times = pwd() * "/$(data_folder)/dwell_times/$(voltage)/"
        data_filenames = cd(readdir, path_data)
            clean_filenames = filter(
            fname -> occursin(r"^ce\d+\.txt$", fname),
            data_filenames
        )
        for data_file in clean_filenames
            data_file_path = path_data * data_file
            data_file_path = path_data * data_file
            dwell_times_file = split(data_file, '.')
            dwell_times_file[1] = dwell_times_file[1] * "dwell_timesy"
            dwell_times_file = join(dwell_times_file, '.')
            dwell_times_path = path_dwell_times * dwell_times_file
            dwell_times_path = path_dwell_times * dwell_times_file
            push!(data_file_paths, data_file_path)
            push!(dwell_times_file_paths, dwell_times_path)
        end
    end
    # what_first_file_name = cd(readdir, pwd() * "/$(data_folder)/")[3]
    what_first_path = pwd() * "/$(data_folder)/first.txt"
    what_first_path, data_file_paths, dwell_times_file_paths
end

"""
    get_specified_datapoints(x::Vector{Float32}, y::Vector{Float32}, Δt::Float32, data_size=-1) 
        -> Dict{String, Vector{Float32}}

Extract a segment of the `x` data vector and the corresponding `dwell times` segment 
up to a specified number of data points or total recording time.

# Arguments
- `x::Vector{Float32}`: Time series or measurement values.
- `y::Vector{Float32}`: Corresponding dwell times between events or state changes.
- `Δt::Float32`: Sampling interval in seconds for the `x` data.
- `data_size::UInt32` (optional, default=`0`):  
Number of samples to include in the result.  
- If `0`, the entire dataset is used.
- Otherwise, selects the first `data_size` samples of `x` and the dwell times covered by them.

# Returns
`Dict{String, Vector{Float32}}`:
- `"x"` → The truncated `x` vector containing the first `data_size` samples (or all samples if `data_size=-1`).
- `"dwell times"` → A truncated version of `y` containing only those dwell time segments whose cumulative sum
does not exceed `max_time = data_size * Δt`.

# Method
1. Determine the number of samples to include (`data_size`), using the full length of `x` if `data_size==-1`.
2. Compute the `max_time` in seconds corresponding to the selected number of points.
3. Include only those dwell time segments from `y` whose cumulative sum is less than or equal to `max_time`.
4. Return both the truncated `x` and the matching truncated `y` in a dictionary.

# Example
```
x, y = read_data("data.txt", "dwell_times.txt")
data = get_specified_datapoints(x, y, Δt, 50)
println(data["x"])
println(data["dwell times"])
```

# Notes
- The dwell times are selected based on cumulative duration, **not** index count.
- This function preserves the original order of the dwell times.
- It assumes that `x` and `y` represent compatible datasets in terms of recording sequence.

"""
function get_specified_datapoints(x::Vector{Float32}, y::Vector{Float32}, Δt::Float32, data_size::UInt32=UInt32(0)) :: Dict{String, Vector{Float32}}
    N = length(x)
    data_size = data_size == 0 || data_size > N ? N : data_size
    max_time = data_size*Δt
    Y = y[findall(t -> t <= max_time, cumsum(y))]
    data = Dict("x" => x[1:data_size], "dwell times" => Y)
    data
end

"""
    normalize_data(data::Dict{String, Vector{Float32}}) -> Vector{Float32}

Normalize the `"x"` values in a data dictionary to have zero mean and unit variance
using z-score normalization.

# Arguments
- `data::Dict{String, Vector{Float32}}`:  
A dictionary containing at least the key `"x"` mapped to a vector of floating-point values
(such as raw measurement data).  
Other keys (e.g., `"dwell times"`) may be present but are ignored.

# Returns
- `Vector{Float32}`:  
A new vector of the same length as `data["x"]`, where each element has been normalized:

z_i = (x_i - μ) / σ

where μ is the mean of `x`, and σ is its standard deviation.

# Method
1. Fit a [`ZScore`] scaling model
to the `"x"` values using `fit(ZScore, data["x"])`.
2. Apply the `normalize` function from `StatsBase` to transform the data into z-scores.
3. Return the transformed vector.

# Example
```
x, y = read_data("data.txt", "dwell_times.txt")
data = get_specified_datapoints(x, y, Δt, 50)
normalized_x = normalize_data(data)
```

# Notes
- Requires the **StatsBase.jl** and **Normalization.jl** package.
- The `"x"` vector must not be empty and must contain finite real values.
- This function does not modify the original dictionary; it returns a new normalized vector.
"""
function normalize_data(data::Dict{String, Vector{Float32}}) :: Vector{Float32}
    N = fit(ZScore, data["x"])
    normalized_data = normalize(data["x"], N)
    normalized_data
end

"""
    combine_time_with_data(data::Vector{Float32}, Δt::Float32, batch_size=1) -> Vector{Tuple{Float32, Float32}}

Create a time-stamped version of raw data, pairing each value with its time in the sampled sequence.

# Arguments
- `data::Vector{Float32}`  
The signal or measurement data to be time-stamped.
- `Δt::Float32`  
The sampling interval (seconds) between consecutive data points.
- `batch_size::UInt8` (optional, default = 1)  
The step size for batch-wise processing; typically leave as 1 for full sequence.

# Returns
- `Vector{Tuple{Float32, Float32}}`  
A vector of `(time, value)` pairs. The time values run from 0 to the end in steps of `Δt * batch_size`, each paired with the matching data value.

# Description
Pairs each data point with its corresponding timestamp, supporting batch-wise access for processing algorithms that operate on downsampled or chunked data.

# Example
```
data = [0.1, 0.2, 0.3, 0.4]
Δt = 0.01
pairs = combine_time_with_data(data, Δt)
```
"""
function combine_time_with_data(data::Vector{Float32}, Δt::Float32, batch_size::UInt8=UInt8(1)) :: Vector{Tuple{Float32, Float32}}
    data_to_process = data[1:batch_size:end]
    data_with_times = collect(zip(0:Δt*batch_size:length(data_to_process), data_to_process))
    data_with_times
end