using StatsBase
using Distributions

"""
create_idealizations(data_folder::String, Δt::Float32) -> Dict{String, Vector{Int8}}

Generate idealized binary traces from raw experimental data and dwell times.

This function reads data files and their associated dwell time files from a given folder
(using [`read_all_file_paths`](@ref)), then produces idealized sequences of `0` and `1` values
based on the dwell times and initial state information.

# Arguments
- `data_folder::String`: Path to the main data directory.  
Must follow the folder structure expected by [`read_all_file_paths`](@ref), with subfolders:
- `sampling/<voltage>/`
- `dwell_times/<voltage>/`
and an additional file at the top level used to determine the starting state for each trace.
- `Δt::Float32` (optional, default=1e-4): Time step in seconds used to convert dwell times to point counts.

# Returns
- `Dict{String, Vector{Int8}}`:  
A dictionary where:
- Keys are **file names** (without path) of the raw data files.
- Values are idealized traces (`Vector{Int8}`), where elements are `0` or `1`.

# Method
1. Calls [`read_all_file_paths`](@ref) to get:
- Path to the "first state" file (`what_first_path`)
- List of raw data file paths
- List of dwell times file paths
2. Reads the initial states for each trace from the "first state" file.
3. For each dataset:
- Reads raw signal (`x`) and dwell times (`y`).
- Constructs a binary vector that alternates between `0` and `1` according to the dwell times.
- Ensures the idealized vector length matches the raw data length by trimming or padding as necessary.
4. Stores the idealized vector in the result dictionary under the base file name of the dataset.

# Example
```
idealized = create_idealizations("experiment_data")
println(keys(idealized)) # List of processed file names
println(idealized["trace01.dat"][1:20]) # First 20 idealized points of a trace
```

# Notes
- The time step `Δt` is fixed at `1e-4` seconds.
- Dwell times are converted to integer point counts via `round(Int, dt/Δt)`.
- If a computed segment length is `0`, it is replaced with a segment of length `1`.
- If the total length exceeds the raw signal length, the idealized trace is truncated;  
otherwise, zeros or ones are appended to match the length.
- Switching between `0` and `1` starts from the value given in the `what_first_path` file
for the specific trace.
"""
function create_idealizations(data_folder::String, Δt::Float32=Float32(1e-4)) :: Dict{String, Vector{Int8}}

    what_first_path, data_paths, dwell_times_paths = read_all_file_paths(data_folder)
    N = length(data_paths)    

    what_first_files = open(what_first_path) do f; parse.(Int8, strip.(readlines(f))); end

    idealized_big_data = Dict{String, Vector{Int8}}([])

    for i in 1:N
    x = open(data_paths[i]) do f; parse.(Float32, strip.(readlines(f))); end
    y = open(dwell_times_paths[i]) do f; parse.(Float32, strip.(readlines(f))); end

    data = Dict("x" => x, "dwell times" => y)
    # idealize data
    what_first = what_first_files[i]
    idealized_value = what_first
    idealized_values = Vector{Int8}([])
    for dt in data["dwell times"]
        how_many = round(Int, dt/Δt)
        append!(idealized_values, idealized_value * ones(how_many != 0 ? how_many : 1))
        idealized_value = idealized_value == 0 ? 1 : 0
    end
    if length(idealized_values) > length(data["x"])
        idealized_values = idealized_values[1:length(data["x"])]
    else
        append!(idealized_values, idealized_value * ones(length(data["x"]) - length(idealized_values)))
    end
    idealized_big_data[split(data_paths[i], "/")[end]] = idealized_values
    end
    idealized_big_data
end

"""
histogram_calculator(data::Vector{Float32}, bins::Int16=-1) -> Histogram

Compute a histogram of the given data vector with Freedman-Diaconis binning or a specified number of bins.

# Arguments
- `data::Vector{Float32}`: A vector of floating-point numbers representing the data to histogram.
- `bins::Int16` (optional, default=-1): Number of bins to divide the data range into. If set to -1, the Freedman-Diaconis rule is used to determine the optimal number of bins.

# Returns
- `Histogram`: A `Histogram` object (from `StatsBase.jl`) representing the frequency distribution
of the data across the specified bins.

# Method
1. Determine the minimum (`min_data`) and maximum (`max_data`) values in the data.
2. Divide the interval [`min_data`, `max_data`] into `bins` equal-width bins.
3. Use `StatsBase.fit(Histogram, data, edges)` to compute the histogram counts and bin edges.
4. Return the histogram object.

# Example
```
using StatsBase

data = randn(1000) # 1000 samples from a normal distribution
hist = histogram_calculator(data, UInt16(50))
println(hist.weights) # Counts per bin
println(hist.edges) # Bin edges
```

# Notes
- Requires `StatsBase.jl` for the `Histogram` type and `fit` function.
- The bins are equally spaced between the minimum and maximum data values.
- The returned `Histogram` object contains bin edges and counts, suitable for further analysis or plotting.
"""
function histogram_calculator(data::Vector{Float32}, nbins::Int16=Int16(-1)) :: Histogram
    extrema_of_data = extrema(data)
    min_data = extrema_of_data[1]
    max_data = extrema_of_data[2]
    if nbins > 0
        edges = range(min_data, stop=max_data, length=nbins+1)
        histogram_of_data = fit(Histogram, data, edges)
        return histogram_of_data
    end
    IQR::Float32 = iqr(data)
	n::UInt32 = length(data)
	bin_width::Float32 = 2.0 * (IQR/∛n)
	number_of_bins = round(Int, (max_data - min_data) / bin_width)
	histogram_of_data = fit(Histogram, data, nbins= number_of_bins)
    histogram_of_data
end

"""
calculate_probability_histogram(histogram::Histogram) -> Histogram

Convert the counts in a histogram to probabilities, producing a probability histogram.

# Arguments
- `histogram::Histogram`  
A histogram object (from `StatsBase`) containing bin edges and weighted counts.

# Returns
- `Histogram`  
A new histogram object with the same bin edges as the input, but with weights normalized to sum to 1, representing probabilities.

# Description
This function takes a histogram of counts or weights and converts it to a probability histogram by dividing each bin’s weight by the total sum of all bin weights.

# Example
```
data = randn(1000)
hist = fit(Histogram, data, 50)
prob_hist = calculate_probability_histogram(hist)

println(sum(prob_hist.weights)) # Should print 1.0 (or very close due to floating point)
```

# Notes 
- TO BE DEPRECATED in the future;
"""
function calculate_probability_histogram(histogram::Histogram) :: Histogram
    normalize(histogram, mode=:pdf)
end

"""
    analyze_histogram_peaks(prob_hist::Histogram) -> HistPeakAnalysis

Analyze a probability histogram to find the indices and values of major peaks and the trough between them.

# Arguments
- `prob_hist::Histogram`  
A probability histogram (from `StatsBase.Histogram`), typically with normalized weights.

# Returns
- [`HistPeakAnalysis`](@ref)  
A structure bundling bin edges, weights, indices and values for the two main peaks, the midpoint index, and the minimum value (trough) between the peaks.

# Description
This function examines the provided histogram to determine the location and values of:
- The left maximum (`pmax1`)
- Its left index (`pmax1_index`)
- The right maximum (`pmax2`)
- Its right index (`pmax2_index`)
- The midpoint index between the two maxima
- The minimum value (`pmin`) found between those peaks (used for thresholding)
All results are packed into a [`HistPeakAnalysis`](@ref) struct for downstream use.

# Example
```
data = randn(1000)
hist = fit(Histogram, data, 50)
prob_hist = calculate_probability_histogram(hist)
analysis = analyze_histogram_peaks(prob_hist)

println("First peak: ", analysis.left_peak, " at index ", analysis.left_peak_index)
println("Second peak: ", analysis.right_peak, " at index ", analysis.right_peak_index)
println("Minimum between peaks: ", analysis.pmin, " at index ", analysis.pmin_index)
```
"""
function analyze_histogram_peaks(prob_hist::Histogram) :: HistPeakAnalysis
    edges = collect(prob_hist.edges[1])
    weights = prob_hist.weights

    # First peak (absolute max)
    pmax1, pmax1_index = findmax(weights)
    half_hist = Int(round(length(edges) / 2))

    if pmax1_index < half_hist
        # Peak on left
        midpoint = floor(Int, (pmax1_index + length(edges)) / 2)
        pmax2, pmax2_index = findmax(weights[(midpoint+1):end])
        pmax2_index += midpoint

        # Min between peaks
        pmin, pmin_index = findmin(weights[pmax1_index:pmax2_index])
        pmin_index += pmax1_index - 1
    else
        # Peak on right
        midpoint = floor(Int, pmax1_index / 2)
        pmax2, pmax2_index = findmax(weights[1:midpoint])

        pmin, pmin_index = findmin(weights[pmax2_index:pmax1_index])
        pmin_index += pmax2_index - 1

        # Swap so pmax1 is always 'left' peak
        pmax1, pmax2 = pmax2, pmax1
        pmax1_index, pmax2_index = pmax2_index, pmax1_index
    end

    return HistPeakAnalysis(edges, weights,
        pmax1, pmax1_index,
        pmax2, pmax2_index,
        midpoint,
        pmin, pmin_index
    )
end

"""
    calculate_method(data::Vector{Float32}, c_method::IdealizationMethod, Δt::Float32)

Run the dwell-time estimation algorithm associated with the given method type.

# Arguments
- `data::Vector{Float32}` - The signal to analyse.
- `m::IdealizationMethod` - Parameters for the chosen method; algorithm resolved by `method(m)`.
- `Δt::Float32` - Sampling interval in seconds.

# Returns
- `Vector{Float32}` - Dwell times detected by the chosen method.

# Example
```
m = MeanDeviationMethod(0.05, 0.5)
dwell_times = calculate_dwell_times(signal, m, 1e-4)  
```
"""
function calculate_method(data::Vector{Float32}, c_method::IdealizationMethod, Δt::Float32)
    method_function(c_method)(data, Δt, c_method)
end