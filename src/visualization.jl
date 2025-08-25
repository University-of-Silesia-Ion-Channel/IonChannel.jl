using Plots
using StatsBase

"""
show_approx_on_plot(data::Dict{String, Vector{Float32}}, method_output::MethodOutput, T_left::Float32, T_right::Float32, Δt::Float32)

Plot raw data and overlay exact (ground truth) and approximated (idealized) breakpoints, with threshold line for Mika methods.

# Arguments
- `data::Dict{String, Vector{Float32}}`  
Dictionary containing `"x"` (raw data signal) and `"dwell times"` (vector of dwell segment durations).
- `method_output::MethodOutput`  
Output from an idealization algorithm, providing breakpoints and (for Mika) threshold information.
- `T_left::Float32`
Start time (seconds) for the plot interval.
- `T_right::Float32`  
End time (seconds) for the plot interval.
- `Δt::Float32`  
Sampling interval (seconds).

# Description
Plots the time series `"x"` for the duration `[T_left, T_right]` as a green line.  
- **Red vertical lines:** Mark true breakpoints, based on cumulative dwell times from `T_left` to `T_right`.
- **Blue vertical lines:** Indicate approximated breakpoints detected by the idealization method, up to from `T_left` to `T_right`.
- **Threshold line:** If the `method_output` is a [`MikaMethodOutput`](@ref), overlays a horizontal line at the chosen threshold.

Line opacities scale with the interval length.  
Additional features—such as band limits or idealized data overlay—can be enabled by uncommenting code sections.

# Example
```
show_approx_on_plot(data, result, 0.5, 0.9, 1e-4)
```  
"""
function show_approx_on_plot(data::Dict{String, Vector{Float32}}, method_output::MethodOutput, T_left::Float32, T_right::Float32, Δt::Float32)
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

    approx_indices = findall(t -> t >= T_left && t <= T_right, method_output.breakpoints)
    approx_breakpoints_to_draw = method_output.breakpoints[approx_indices]

    alpha_val = min(0.5, 1 / T_diff)
    vline!(breakpoints_to_draw; alpha=alpha_val, label="Exact", color="red")

    if typeof(method_output) <: MikaMethodOutput
        hline!([method_output.threshold], label="threshold", lw=2)
    end

    vline!(approx_breakpoints_to_draw, label="Approximated", alpha=alpha_val, color="blue")

end


"""
    plot_idealization_representation(data::Dict{String, Vector{Float32}}, method_output::MethodOutput, T_left::Float32, T_right::Float32, Δt::Float32)

Visualize ion channel data and its idealization using stacked subplots.  
Supports [`MikaMethodOutput`](@ref) and [`MeanDeviationMethodOutput`](@ref).

# Arguments
- `data::Dict{String, Vector{Float32}}`  
Dictionary containing `"x"` (raw signal data) and associated fields as needed.
- `method_output::MethodOutput`  
Output of an idealization algorithm. If it is [`MikaMethodOutput`](@ref), uses the built-in idealized data; if [`MeanDeviationMethodOutput`](@ref), reconstructs the idealization using dwell times and histogram analysis.
- `T_left::Float32`
Start time (seconds) for the plot interval.
- `T_right::Float32`  
Total plot duration (seconds).
- `Δt::Float32`  
Sampling interval (seconds).

# Description
Creates a two-row subplot:
- **Top plot:** Raw ion channel current signal (`data["x"]`) in green.
- **Bottom plot:** Idealized signal in blue.
    - For Mika methods: uses `method_output.idealized_data`.
    - For Mean Deviation methods: reconstructs idealization using dwell times and peak analysis.

Axes are labeled by time and current.  
Legend is suppressed for clarity.  
Other idealization outputs can be supported by extending the function.

# Example
```
# For Mika idealization
plot_idealization_representation(data, mika_output, 0.5, 1.0, 1e-4)

# For mean deviation idealization
plot_idealization_representation(data, mean_deviation_output, 0.5, 1.0, 1e-4)
```

# Note
For [`MeanDeviationMethodOutput`](@ref), histogram calculation and peak analysis are rerun to reconstruct the idealized trace.
"""
function plot_idealization_representation(data::Dict{String, Vector{Float32}}, method_output::MethodOutput, T_left::Float32, T_right::Float32, Δt::Float32)
    @assert T_left <= T_right "N_left must be less or equal to N_right"
    @assert T_right <= (length(data["x"]) - 1) * Δt "T_right exceeds data duration"

    N_left = max(1, Int(round(T_left / Δt)) + 1)
    N_right = min(length(data["x"]), Int(round(T_right / Δt)) + 1)
    N = N_right - N_left + 1

    time = range(T_left, T_left + Δt*(N-1), length=N)

    y1 = data["x"][N_left:N_right]

    y2 = method_output.idealized_data[N_left:N_right]

    plt1 = plot(time, y1, color=:green, legend=false, title="Ion channel current plot", dpi=200)
    plt2 = plot(time, y2, color=:blue, legend=false, title="Idealization of ion channel", dpi=200)

    plot(plt1, plt2, layout=grid(2, 1, heights=[0.75, 0.25]))
    xlabel!("time [s]")
    ylabel!("current [pA]")

end

"""
    show_threshold_on_plot(data_histogram::Histogram, histogram_analysis::HistPeakAnalysis, method_output::MikaMethodOutput)

Plot the data histogram and visualize threshold-related features and detected peaks.

# Arguments
- `data_histogram::Histogram`  
The histogram of original data, typically created with [`histogram_calculator`](@ref).
- `histogram_analysis::HistPeakAnalysis`  
Results from peak and trough detection, providing bin indices for annotated visualization.
- `method_output::MikaMethodOutput`  
Output of the Mika method containing the chosen threshold index and other results.

# Description
Creates a bar plot of the data histogram and overlays vertical lines to indicate:
- the primary and secondary peaks,
- the histogram midpoint,
- the selected threshold location,
- and the minimum (trough) between peaks.

This visualization helps evaluate the effectiveness and placement of the threshold and key features used by the idealization algorithm.

# Example
```
show_threshold_on_plot(data_histogram, analysis, result)
```
"""
function show_threshold_on_plot(data_histogram::Histogram, histogram_analysis::HistPeakAnalysis, method_output::MikaMethodOutput)
    bar(data_histogram, label="Histogram (density)", alpha=0.5)
    vline!([histogram_analysis.edges[histogram_analysis.pmax1_index]], label="pmax1")
    vline!([histogram_analysis.edges[histogram_analysis.pmax2_index]], label="pmax2")
    vline!([histogram_analysis.edges[histogram_analysis.midpoint]], label="half")
    vline!([histogram_analysis.edges[method_output.threshold_index]], label="threshold", lw=3)
    vline!([histogram_analysis.edges[histogram_analysis.pmin_index]], label="minimum")
    # r1 = collect(prob_hist_edges_list[pmax1_index - 1]:0.1:prob_hist_edges_list[pmin_index + 2])
    # plot!(r1, line1[1].*r1 .+ line1[2], label="Line 1")
    # r2 = collect(prob_hist_edges_list[pmin_index-1]:0.1:prob_hist_edges_list[pmax2_index+1])
    # plot!(r2, line2[1].*r2 .+ line2[2], label="Line 2")
    # vline!([x₁], label="x1", linewidth=2)
    # vline!([x₂], label="x2", linewidth=2)
end

