using Distributions
using StatsBase
using Markdown
using DataFrames

"""
    calculate_mean_square_error(
        data::Dict{String, Vector{Float32}}, 
        dwell_times_approx::Vector{Float32}, 
        dt_bins::UInt16=100
    ) -> Tuple{Float32, Histogram, Histogram}

Calculate the mean squared error (MSE) between histograms of actual and approximate dwell times, 
and return both histograms for further analysis.

# Arguments
- `data::Dict{String, Vector{Float32}}`  
A dictionary containing at least the key `"dwell times"` mapped to a vector of observed dwell times (Float32).
- `dwell_times_approx::Vector{Float32}`  
A vector of approximate dwell times to compare against the observed data.
- `dt_bins::UInt16` (optional, default=`100`)  
Number of bins to use when computing histograms for error calculation.

# Returns
`Tuple{Float32, Histogram, Histogram}` — `(mse, hist_data, hist_approx)` where:
1. `mse::Float32`  
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
function calculate_mean_square_error(data::Dict{String, Vector{Float32}}, dwell_times_approx::Vector{Float32}, dt_bins::UInt16=UInt16(100)) :: Tuple{Float32, Histogram, Histogram}
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
    idealize_data(data::Vector{Float32}, dwell_times_approx::Vector{Float32}, hist_analysis::HistPeakAnalysis, Δt::Float32) -> Vector{UInt8}

Construct an idealized signal from approximate dwell times and histogram peak analysis.

# Arguments
- `data::Vector{Float32}`  
The original raw data signal.
- `dwell_times_approx::Vector{Float32}`  
Approximate dwell times computed by an idealization method, in seconds.
- `hist_analysis::HistPeakAnalysis`  
Result of histogram peak analysis providing peak intensities and midpoint for thresholding.
- `Δt::Float32`  
Sampling interval of the data, in seconds.

# Returns
- `Vector{UInt8}`  
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
function idealize_data(data::Vector{Float32}, dwell_times_approx::Vector{Float32}, hist_analysis::HistPeakAnalysis, Δt::Float32) :: Vector{Float32}
    Imid = hist_analysis.edges[hist_analysis.midpoint]
    Imax1 = hist_analysis.edges[hist_analysis.pmax1_index]
    Imax2 = hist_analysis.edges[hist_analysis.pmax2_index]
    Imax1, Imax2 = Imax1 > Imax2 ? (Imax2, Imax1) : (Imax1, Imax2)
    # println("Imid: ", Imid, ", Imax1: ", Imax1, ", Imax2: ", Imax2)
    idealized_value = data[1] < Imid ? Imax1 : Imax2
    # println("Initial idealized value: ", idealized_value)
    idealized_values = Vector{Float32}([])
    # println("Dwell times [1]: ", dwell_times_approx[1]/Δt)
    for dt in dwell_times_approx
        how_many = round(UInt16, dt/Δt)
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

"""
    actual_idealize_data(
        data::Dict{String, Vector{Float32}},
        what_first_dict::Dict{String, Int64},
        data_file_name::AbstractString,
        Δt::Float32
    ) :: Vector{UInt8}

Reconstruct a binary idealized signal from dwell-time information by alternating
between states `0` and `1` with durations matching the provided dwell times.

# Arguments
- `data::Dict{String, Vector{Float32}}`: Dictionary containing at least:
    - `"dwell times"`: a vector of dwell time durations for consecutive states.
    - `"x"`: the original data trace (used to match output length).
- `what_first_dict::Dict{String, Int64}`: Dictionary specifying the starting
   state for each file, mapped by file name.
- `data_file_name::AbstractString`: Name of the dataset; used to select
   the initial state from `what_first_dict`.
- `Δt::Float32`: Sampling interval of the signal; used to convert dwell
   durations into sample counts.

# Returns
- `Vector{UInt8}`: A binary idealized trace of the same length as `data["x"]`,
   representing predicted open/closed states (`0`/`1`).

# Method
1. Retrieve the starting state (`what_first`) from `what_first_dict` for the
   current file.
2. Initialize `idealized_values` with that state.
3. For each dwell time in `data["dwell times"]`:
   - Convert dwell duration into a number of samples: `how_many = round(Int, dt / Δt)`.
   - Append that many samples of the current state.
   - Switch to the opposite state (`0 → 1` or `1 → 0`).
   - Ensure at least one sample is generated for very small dwell times.
4. Adjust the length of the result:
   - If longer than `data["x"]`, truncate.
   - If shorter, pad with the next state to match the signal length.
5. Return the finalized binary idealized trace.

# Notes
- The function enforces that the returned vector always matches
  the length of the original recorded trace (`data["x"]`).
- Alternating states assumes a two-state system
  (`0` and `1`), which is standard in patch-clamp idealization.
- Handles dwell durations shorter than `Δt` by ensuring at least one sample.
"""
function actual_idealize_data(data::Dict{String, Vector{Float32}}, what_first_dict::Dict{String, Int64}, data_file_name::AbstractString, Δt::Float32) :: Vector{UInt8}
	what_first = what_first_dict[data_file_name]
	idealized_value = what_first

	idealized_values = Vector{UInt8}([])
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

"""
    accuracy_of_idealization(
        actual_idealization::Vector{UInt8},
        approx_idealization::Vector{UInt8}
    ) :: Float32

Compute the accuracy of an approximate idealization compared to the ground-truth
(or reference) idealization.

# Arguments
- `actual_idealization::Vector{UInt8}`: The reference (true) idealized signal,
   represented as a binary vector (`0` and `1` states).
- `approx_idealization::Vector{UInt8}`: The test or approximated idealized trace
   to be evaluated.

# Returns
- `Float32`: The accuracy, computed as the proportion of matching states
  between the two sequences, in the range `[0.0, 1.0]`.

# Method
1. If the first state of `approx_idealization` does not match the first state of
   `actual_idealization`, then `approx_idealization` is inverted
   (`0 ↔ 1`) to ensure label consistency.
2. Accuracy is calculated as:
    ```
    sum(actual_idealization .== approx_idealization) / length(approx_idealization)
    ```
    i.e., the fraction of samples where both sequences agree.

# Notes
    - This function assumes both inputs are of equal length.
    - If a systematic state-label swap occurred (e.g., model outputs `1` for "closed"
    but reference uses `0`), the inversion step ensures a fair comparison.
    - Accuracy is a simple per-sample metric and may not capture temporal
    misalignments (e.g., small shifts in breakpoints).
"""
function accuracy_of_idealization(actual_idealization::Vector{UInt8}, approx_idealization::Vector{UInt8}) :: Float32
    # Calculate the accuracy as the proportion of matching states
    if actual_idealization[1] != approx_idealization[1]
        approx_idealization = 1 .- approx_idealization
    end
	sum(actual_idealization .== approx_idealization) / length(approx_idealization)
end

"""
    mean_error(method::IdealizationMethod, Δt::Float32, data_size::UInt32, ::Bool=false) -> Float32

Compute the **average mean squared error (MSE)** across multiple datasets,
using a specified idealization method to approximate dwell times.

# Arguments
- `method::IdealizationMethod`  
An instance of a concrete subtype of [`IdealizationMethod`](@ref),  
such as [`MeanDeviationMethod`](@ref), which stores:
    - The dwell-time estimation function.
    - Method parameters (e.g. δ, λ values).
- `Δt::Float32`  
Sampling interval (seconds) of the recordings.
- `data_size::UInt32`
Number of data points to include in each dataset for MSE calculation.
- `verbose::Bool=false`
If `true`, prints detailed processing information for each dataset.

# Returns
- `Float32`:  
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
the [`IdealizationMethod`](@ref) instance and calling its stored method function.
- The `"x"` and `"dwell times"` vectors are assumed to be aligned and compatible.
- [`calculate_mean_square_error`](@ref) now returns a tuple; only the first element is used.

# Example
```
m = MeanDeviationMethod(deviation_from_mean_method, 0.05, 0.5)
Δt = 1e-4
avg_mse = mean_error(m, Δt, UInt32(10000))
println("Average MSE across datasets: ", avg_mse)
```
"""
function mean_error(method::IdealizationMethod, Δt::Float32, data_size::UInt32, verbose::Bool=false)
    what_first_file_path, data_paths, dwell_times_paths = read_all_file_paths("data")
    data_paths_dict = create_paths_dictionary(data_paths, dwell_times_paths)

    mean_squared_errors = Dict{String, Dict{String, Union{Vector{Float32}, Float32}}}([])
    accuracies = Dict{String, Dict{String, Union{Vector{Float32}, Float32}}}([])

    table = Dict{String, Dict{String, Vector{Float32}}}(["errors" => Dict{String, Vector{Float32}}(), "accuracies" => Dict{String, Vector{Float32}}()])
    mean_error_dict = Dict{String, Float32}()
    mean_accuracy_dict = Dict{String, Float32}()
    if verbose
        @info "$(what_first_file_path)"
    end
    what_first_dict = Dict(
	    String(split(line,',')[1]) => parse(Int, split(line,',')[2])
	    for line in eachline(what_first_file_path)
	)

    for voltage in keys(data_paths_dict["data paths"])
        if verbose
            @info "Processing voltage $(voltage)"
        end
        # mean_squared_errors[voltage] = Dict{String, Union{Vector{Float32}, Float32}}("errors" => Vector{Float32}([]), "mean error" => 0.0f0)
        # accuracies[voltage] = Dict{String, Union{Vector{Float32}, Float32}}("accuracies" => Vector{Float32}([]), "mean accuracy" => 0.0f0)
        # table["errors"] = Dict{String, Vector{Float32}}(voltage => Float32[])
        # table["accuracies"] = Dict{String, Vector{Float32}}(voltage => Float32[])
        N = length(data_paths_dict["data paths"][voltage])
        temp_error = 0.0f0
        temp_acc = 0.0f0
        acc_table = Float32[]
        errors_table = Float32[]
        for i in 1:N
            x, y = read_data(data_paths_dict["data paths"][voltage][i], data_paths_dict["dwell times paths"][voltage][i])
            data = get_specified_datapoints(x, y, Δt, data_size)
            normalized_data = normalize_data(data)
            data["x"] = normalized_data
            method_output = calculate_method(normalized_data, method, Δt)

            mse = calculate_mean_square_error(data, method_output.dwell_times_approx)[1]
            # push!(mean_squared_errors[voltage]["errors"], mse)
            # push!(table["errors"][voltage], mse)
            temp_error += mse
            push!(errors_table, mse)
            actual_idealized_data = actual_idealize_data(data, what_first_dict, split(data_paths_dict["data paths"][voltage][i], '/')[end], Δt)
            if typeof(method_output) <: MikaMethodOutput
                vals = sort(unique(method_output.idealized_data))
                mapped = (method_output.idealized_data .== vals[2])
                approx_idealization = Vector{UInt8}(mapped)
            else
                approx_idealization = method_output.idealized_data
            end
            acc = accuracy_of_idealization(actual_idealized_data, approx_idealization)
            push!(acc_table, acc)
            # push!(accuracies[voltage]["accuracies"], acc)
            # push!(table["accuracies"][voltage], acc)
            temp_acc += acc
        end
        table["errors"][voltage] = errors_table
        table["accuracies"][voltage] = acc_table
        mean_error_dict[voltage] = temp_error / N
        mean_accuracy_dict[voltage] = temp_acc / N
        # mean_squared_errors[voltage]["mean error"] = temp_error / N
        # accuracies[voltage]["mean accuracy"] = temp_acc / N
    end
    # table, accuracies, mean_squared_errors
    table, mean_accuracy_dict, mean_error_dict
end

"""
    dicts_to_dataframes(table::Dict{String,Dict{String,Vector{Float32}}},
                        mean_accuracy_dict::Dict{String,Float32},
                        mean_error_dict::Dict{String,Float32})

Return (df_errors, df_accuracies, df_summary) as DataFrames.

- df_errors: columns per voltage with MSE values (missing padded).
- df_accuracies: columns per voltage with accuracy values (missing padded).
- df_summary: one row per voltage with mean_error and mean_accuracy.
"""
function dicts_to_dataframes(table::Dict{String,Dict{String,Vector{Float32}}},
                                mean_accuracy_dict::Dict{String,Float32},
                                mean_error_dict::Dict{String,Float32})

    # Helper to convert Dict{String,Vector{Float32}} -> DataFrame with missing padding
    function vector_dict_to_df(d::Dict{String,Vector{Float32}})
        # Determine maximum length among all vectors
        maxlen = isempty(d) ? 0 : maximum(length.(values(d)))
        # Build a NamedTuple of columns with element type Union{Missing,Float32}
        cols = (; (Symbol(k) => Union{Missing,Float32}[ i <= length(v) ? v[i] : missing
                                                for i in 1:maxlen ]
                    for (k, v) in d)...)
        DataFrame(cols)
    end

    df_errors = vector_dict_to_df(table["errors"])
    df_accuracies = vector_dict_to_df(table["accuracies"])

    # Align voltages using keys of mean_error_dict as the canonical set
    voltages = collect(keys(mean_error_dict))
    df_summary = DataFrame(
        voltage = voltages,
        mean_error = Float32[ mean_error_dict[v] for v in voltages ],
        mean_accuracy = Float32[ get(mean_accuracy_dict, v, NaN32) for v in voltages ],
    )

    return df_errors, df_accuracies, df_summary
end