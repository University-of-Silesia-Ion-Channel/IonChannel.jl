using Distributions
using StatsBase
using Markdown
using PlutoUI

"""
    calculate_mean_square_error(
        data::Dict{String, Vector{Float64}}, 
        dwell_times_approx::Vector{Float64}, 
        dt_bins::Int=100
    ) -> Tuple{Float64, Histogram, Histogram}

Calculate the mean squared error (MSE) between histograms of actual and approximate dwell times, 
and return both histograms for further analysis.

# Arguments
- `data::Dict{String, Vector{Float64}}`  
A dictionary containing at least the key `"dwell times"` mapped to a vector of observed dwell times (Float64).
- `dwell_times_approx::Vector{Float64}`  
A vector of approximate dwell times to compare against the observed data.
- `dt_bins::Int` (optional, default=`100`)  
Number of bins to use when computing histograms for error calculation.

# Returns
`Tuple{Float64, Histogram, Histogram}` — `(mse, hist_data, hist_approx)` where:
1. `mse::Float64`  
Mean squared error between the histogram bin counts of actual and approximate dwell times.  
Computed as the average squared difference in bin weights.
2. `hist_data::Histogram`  
Histogram of actual dwell times from `data["dwell times"]`.
3. `hist_approx::Histogram`  
Histogram of `dwell_times_approx`.

# Description
1. Creates histograms for both **actual** (`data["dwell times"]`) and **approximate** dwell times using [`histogram_calculator`](@ref).
2. Fits exponential distributions to each dataset and computes differences in the fitted scale parameters (`θ` values).
3. Calculates the MSE between the histogram weights of the actual and approximate dwell times.
4. Returns the MSE along with both histogram objects for plotting or further analysis.

# Notes
- Uses `StatsBase` for histogram handling and `fit_mle` from `Distributions` for maximum likelihood estimation.
- The intermediate variables `breakpoints` and `accuracy` are currently computed but not returned.
- The returned histograms can be directly plotted with `Plots.jl` or analyzed further.

# Example
```
using Distributions, StatsBase

data = Dict("dwell times" => rand(Exponential(1.0), 1000))
dwell_times_approx = rand(Exponential(1.1), 1000)

mse, h_data, h_approx = calculate_mean_square_error(data, dwell_times_approx, 50)

println("Mean squared error: ", mse)
display(h_data)
display(h_approx)
```
"""
function calculate_mean_square_error(data::Dict{String, Vector{Float64}}, dwell_times_approx::Vector{Float64}, dt_bins=100) :: Tuple{Float64, Histogram, Histogram}
    breakpoints = cumsum(dwell_times_approx)
    h_dwell_times = histogram_calculator(data["dwell times"], dt_bins)
    h_dwell_times_approx = histogram_calculator(dwell_times_approx, dt_bins)
    θ_approx = fit_mle(Exponential, dwell_times_approx).θ
    θ_data = fit_mle(Exponential, data["dwell times"]).θ
    error = abs(θ_data - θ_approx)
    accuracy = round((θ_data - error)/θ_data; digits=8)
    mean²error = sum((h_dwell_times.weights .- h_dwell_times_approx.weights).^2) / dt_bins
    mean²error, h_dwell_times, h_dwell_times_approx
end

"""
    idealize_data(data::Vector{Float64}, dwell_times_approx::Vector{Float64}, hist_analysis::HistPeakAnalysis, Δt::Float64) -> Vector{Float64}

Construct an idealized signal from approximate dwell times and histogram peak analysis.

# Arguments
- `data::Vector{Float64}`  
The original raw data signal.
- `dwell_times_approx::Vector{Float64}`  
Approximate dwell times computed by an idealization method, in seconds.
- `hist_analysis::HistPeakAnalysis`  
Result of histogram peak analysis providing peak intensities and midpoint for thresholding.
- `Δt::Float64`  
Sampling interval of the data, in seconds.

# Returns
- `Vector{Float64}`  
Idealized signal reconstructed by alternating between the two states 0 and 1,
segmented according to the provided dwell times and aligned in length with the original data.

# Description
This function assigns each dwell segment one of two peak intensities depending on which side of the midpoint the starting value lies. It alternates between these intensities for consecutive dwell segments. The output length matches the input data by truncation or padding.

# Example
```
ideal_signal = idealize_data(data, dwell_times_approx, hist_analysis, Δt)
plot(data, label="Raw data")
plot!(ideal_signal, label="Idealized signal")
```
"""
function idealize_data(data::Vector{Float64}, dwell_times_approx::Vector{Float64}, hist_analysis::HistPeakAnalysis, Δt::Float64) :: Vector{Float64}
    Imid = hist_analysis.edges[hist_analysis.midpoint]
    Imax1 = hist_analysis.edges[hist_analysis.pmax1_index]
    Imax2 = hist_analysis.edges[hist_analysis.pmax2_index]
    Imax1, Imax2 = Imax1 > Imax2 ? (Imax2, Imax1) : (Imax1, Imax2)
    # println("Imid: ", Imid, ", Imax1: ", Imax1, ", Imax2: ", Imax2)
    idealized_value = data[1] < Imid ? Imax1 : Imax2
    # println("Initial idealized value: ", idealized_value)
    idealized_values = []
    # println("Dwell times [1]: ", dwell_times_approx[1]/Δt)
    for dt in dwell_times_approx
        how_many = round(Int, dt/Δt)
        append!(idealized_values, idealized_value * ones(how_many != 0 ? how_many : 1))
        idealized_value = idealized_value == Imax1 ? Imax2 : Imax1
    end

    if length(idealized_values) > length(data)
        idealized_values = idealized_values[1:length(data)]
    else
        append!(idealized_values, idealized_value * ones(length(data) - length(idealized_values)))
    end
    idealized_values
end

function actual_idealize_data(data::Dict{String, Vector{Float64}}, what_first_dict::Dict{String, Int64}, data_file_name::AbstractString, Δt::Float64) :: Vector{Int8}
	what_first = what_first_dict[data_file_name]
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
	idealized_values
end	

function accuracy_of_idealization(actual_idealization::Vector{Int8}, approx_idealization::Vector{Int8}) :: Float64
    # Calculate the accuracy as the proportion of matching states
    if actual_idealization[1] != approx_idealization[1]
        approx_idealization = 1 .- approx_idealization
    end
	sum(actual_idealization .== approx_idealization) / length(approx_idealization)
end

"""
    mean_error(method::IdealizationMethod, Δt::Float64, data_size::UInt32, ::Bool=false) -> Float64

Compute the **average mean squared error (MSE)** across multiple datasets,
using a specified idealization method to approximate dwell times.

# Arguments
- `method::IdealizationMethod`  
An instance of a concrete subtype of [`IdealizationMethod`](@ref),  
such as [`MeanDeviationMethod`](@ref), which stores:
    - The dwell-time estimation function.
    - Method parameters (e.g. δ, λ values).
- `Δt::Float64`  
Sampling interval (seconds) of the recordings.
- `data_size::UInt32=50` (optional, default=`50`)
Number of data points to include in each dataset for MSE calculation.
- `verbose::Bool=false`
If `true`, prints detailed processing information for each dataset.

# Returns
- `Float64`:  
The average MSE between actual dwell times and those estimated by `method`,
computed over all matching datasets found in the `"data"` folder.

# Description
1. Uses [`read_all_file_paths`](@ref) to find all raw data and dwell time files.
2. For each dataset:
- Loads raw data (`x`) and actual dwell times (`y`).
- Truncates data to desired length via [`get_specified_datapoints`](@ref).
- Normalizes the `"x"` signal using [`normalize_data`](@ref).
- Estimates dwell times by calling [`calculate_method`](@ref)  
    with the normalized data, `method`, and `Δt`.
- Computes the MSE between actual and estimated dwell times
    via [`calculate_mean_square_error`](@ref), taking only the MSE value
    (first element of its tuple return).
3. Averages the per-dataset MSE values.

# Notes
- This function assumes the folder `"data"` exists and has the expected  
subfolder structure required by [`read_all_file_paths`](@ref).
- The helper [`calculate_method`](@ref) is responsible for interpreting  
the `IdealizationMethod` instance and calling its stored method function.
- The `"x"` and `"dwell times"` vectors are assumed to be aligned and compatible.
- `calculate_mean_square_error` now returns a tuple; only the first element is used.

# Example
```
m = MeanDeviationMethod(deviation_from_mean_method, 0.05, 0.5)
Δt = 1e-4
avg_mse = mean_error(m, Δt)
println("Average MSE across datasets: ", avg_mse)
```
"""
function mean_error(method::IdealizationMethod, Δt::Float64, data_size::UInt32, verbose::Bool=false) :: MeanError
    # read all files
    what_first_file_path, data_paths, dwell_times_paths = read_all_file_paths("data")
    N = length(data_paths)
    sum_mean_squared_error = 0.0
    mean_squared_errors = []
    sum_accuracy = 0.0
    accuracies = []
    if verbose
        @info "$(what_first_file_path)"
    end
    what_first_dict = Dict(
	    String(split(line,',')[1]) => parse(Int, split(line,',')[2])
	    for line in eachline(what_first_file_path)
	)

    
    for i in 1:N
        # if verbose
        #     push!(messages, Markdown.md("----------------------------------------"))
        #     push!(messages, Markdown.md("Processing file $(data_paths[i])"))  
        # end
        if verbose
            @info "Processing file $(data_paths[i])"
        end
        x, y = read_data(data_paths[i], dwell_times_paths[i])
        data = get_specified_datapoints(x, y, Δt, data_size)
        normalized_data = normalize_data(data)
        data["x"] = normalized_data
        method_output = calculate_method(normalized_data, method, Δt)

        mse = calculate_mean_square_error(data, method_output.dwell_times_approx)[1]
        push!(mean_squared_errors, mse)

        actual_idealized_data = actual_idealize_data(data, what_first_dict, split(data_paths[i], '/')[end], Δt)
        if typeof(method_output) <: MikaMethodOutput
			vals = sort(unique(method_output.idealized_data))
			mapped = (method_output.idealized_data .== vals[2])
			approx_idealization = Vector{Float64}(mapped)
		else
			approx_idealization = method_output.idealized_data
		end
        acc = accuracy_of_idealization(actual_idealized_data, approx_idealization)
        sum_accuracy += acc
        push!(accuracies, acc)
        if verbose
            @info "Mean squared error: $mse"
            @info "Accuracy of idealization: $acc"
        end
        sum_mean_squared_error += mse
    end
    mean_squared_error = sum_mean_squared_error / N
    mean_accuracy = sum_accuracy / N
    MeanError(mean_squared_error, mean_accuracy, mean_squared_errors, accuracies)
end

