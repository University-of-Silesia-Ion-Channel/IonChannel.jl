"""
    naive_method(
        data::Vector{Float32},
        Δt::Float32,
        method::NaiveMethod
    ) :: NaiveMethodOutput

Perform a simple threshold-based idealization of a time-series by detecting
crossings between two histogram-derived peaks.

This method:
1) builds a histogram of the signal,
2) converts it to a probability histogram,
3) analyzes peak structure to estimate a threshold between two modes,
4) detects state transitions (0 ↔ 1) when the signal crosses that threshold,
5) returns breakpoints, dwell times, and the binary idealized state sequence.

# Arguments
- `data::Vector{Float32}`: The raw time-series trace to be idealized.
- `Δt::Float32`: Sampling interval used to assign times to samples and compute dwell durations.
- `method::NaiveMethod`: Configuration object that must provide `number_of_histogram_bins`
  and supports the auxiliary functions used here (see Notes).

# Returns
- `NaiveMethodOutput`: A struct with fields:
  - `breakpoints::Vector{Float32}`: Times at which state transitions occur.
  - `dwell_times::Vector{Float32}`: Durations between successive breakpoints
    (first dwell time is `breakpoints[1]`).
  - `idealized_data::Vector{Int}`: Binary sequence of states per sample (`0` or `1`).

# Method
1. Compute histogram and probability histogram:
   - `histogram_of_data = histogram_calculator(data, method.number_of_histogram_bins)`
   - `prob_hist = calculate_probability_histogram(histogram_of_data)`
2. Estimate threshold from histogram peak analysis:
   - `hist_analysis = analyze_histogram_peaks(prob_hist)`
   - `threshold = hist_analysis.edges[hist_analysis.pmin_index]`
3. Combine times with data:
   - `data_with_times = combine_time_with_data(data, Δt)` (expected shape `N×2`, columns: time, value)
4. Initialize state from first sample relative to `threshold`.
5. Iterate through samples:
   - Detect upward crossings (0→1) when previous < threshold and current > threshold.
   - Detect downward crossings (1→0) when previous > threshold and current < threshold.
   - Record transition times in `breakpoints` and append `current_state` to `idealized_data`.
6. Compute dwell times:
   - `dwell_times = append!([breakpoints[1]], diff(breakpoints))`
7. Return `NaiveMethodOutput(breakpoints, dwell_times, idealized_data)`.

# Notes
- The following helper functions are expected to be available in scope:
  - `histogram_calculator(data::AbstractVector, nbins::Integer)`
  - `calculate_probability_histogram(histogram)`
  - `analyze_histogram_peaks(prob_hist)` returning at least `edges` and `pmin_index`
  - `combine_time_with_data(data, Δt)` returning a 2-column array `[time value]`
- State labeling convention: below threshold → `0`, above threshold → `1`.
- `idealized_data` is constructed sample-by-sample; breakpoints are continuous-time instants.
- The code assumes at least one crossing; if `breakpoints` is empty, computing
  `dwell_times` as written will error. Consider guarding this case in production.

# Example
```
trace = rand(Float32, 10_000)
Δt = 1f-4
method = NaiveMethod(number_of_histogram_bins = 100)

result = naive_method(trace, Δt, method)

println("Breakpoints: ", result.breakpoints)
println("Dwell times: ", result.dwell_times)
println("First 20 states: ", result.idealized_data[1:20])
```
"""
function naive_method(data::Vector{Float32}, Δt::Float32, method::NaiveMethod) :: NaiveMethodOutput
    # accessor functions for point for better readability
    value(point) = point[2]
    time(point) = point[1]

    histogram_of_data = histogram_calculator(data, method.number_of_histogram_bins)
    prob_hist = calculate_probability_histogram(histogram_of_data)
    hist_analysis = analyze_histogram_peaks(prob_hist)

    threshold = hist_analysis.edges[hist_analysis.pmin_index]
    data_with_times = combine_time_with_data(data, Δt)

    breakpoints = []
    previous_point = data_with_times[1]
    if value(previous_point) < threshold
        current_state = 0 # starting at the bottom
    else
        current_state = 1 # starting at the top
    end

    idealized_data = [current_state]
    for point in data_with_times[2:end, :]
        if current_state == 0
            if value(previous_point) < threshold && value(point) > threshold
                push!(breakpoints, time(point))
                current_state = 1
            end
        else
            if value(point) < threshold && value(previous_point) > threshold
                push!(breakpoints, time(point))
                current_state = 0
            end
        end
        # change the previous point to the next one
        push!(idealized_data, current_state)
        previous_point = point
    end
    dwell_times = append!([breakpoints[1]], diff(breakpoints))
    NaiveMethodOutput(breakpoints, dwell_times, idealized_data)
end

"""
    method_function(::NaiveMethod) -> Function

Return the algorithm function associated with a `NaiveMethod`.

This allows code like `calculate_method(data, m, Δt)` to work for any
`IdealizationMethod` subtype without changing the executor logic.
"""
method_function(::NaiveMethod) = naive_method