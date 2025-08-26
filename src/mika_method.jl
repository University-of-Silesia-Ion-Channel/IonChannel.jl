using Statistics
using Distributions
using StatsBase
using HypothesisTests

"""
    method_function(::MikaMethod) -> Function

Return the algorithm function associated with a [`MikaMethod`](@ref).

This allows code like `calculate_method(data, m, Δt)` to work for any
[`IdealizationMethod`](@ref) subtype without changing the executor logic.
"""
method_function(::MikaMethod) = mika_method

"""
    point(hist::HistPeakAnalysis, indexfield::Symbol, valuefield::Symbol) -> Point

Extract a point (x, y) from a histogram analysis result using field names.

# Arguments
- `hist::HistPeakAnalysis`  
A structurecontaining edges and value fields.
- `indexfield::Symbol`  
The field name indicating the index (e.g., `:pmax1_index`).
- `valuefield::Symbol`  
The field name indicating the value (e.g., `:pmax1`).

# Returns
- [`Point`](@ref)  
A [`Point`](@ref) instance where `x` is taken from `hist.edges` at the specified index,
and `y` is the field value accessed from `hist`.

# Description
This function provides a generic way to extract coordinates from a histogram analysis result, 
allowing flexible selection of peak or trough points for plotting or further calculations.

# Example
```
analysis = analyze_histogram_peaks(prob_hist)
peak_pt = point(analysis, :pmax1_index, :pmax1)
println(peak_pt.x, ", ", peak_pt.y)
```
"""
function point(hist::HistPeakAnalysis, indexfield::Symbol, valuefield::Symbol)::Point
    x = hist.edges[getfield(hist, indexfield)]
    y = getfield(hist, valuefield)
    Point(x, y)
end

"""
    line(point1::Point, point2::Point) -> Line

Construct a [`Line`](@ref) (y = a*x + b) passing through two points.

# Arguments
- `point1::Point`  
The first point `(x₁, y₁)` through which the line will pass.
- `point2::Point`  
The second point `(x₂, y₂)` through which the line will pass.

# Returns
- [`Line`](@ref)  
A line in slope-intercept form (`y = a*x + b`), where `a` is the slope and `b` is the intercept.

# Description
Computes the slope (`a`) and y-intercept (`b`) for the line passing through two supplied points.

# Example
```
p1 = Point(1.0, 2.0)
p2 = Point(3.0, 5.0)
l = line(p1, p2)
println("y = \$(l.a)x + \$(l.b)")
```
"""
function line(point1::Point, point2::Point)::Line
    a = (point1.y - point2.y) / (point1.x - point2.x)
    b = point1.y - a * point1.x
    Line(a, b)
end

"""
    get_threshold_width(hist_analysis::HistPeakAnalysis, ϵ::Float32) -> ThresholdWidth

Compute the threshold band used for state discrimination in idealization, based on histogram analysis and a weighting parameter.

# Arguments
- `hist_analysis::HistPeakAnalysis`  
Result of peak analysis on a histogram, containing peak and minimum bin indices and values.
- `ϵ::Float32`  
Weighting parameter that adjusts the positions of threshold bounds between minimum and peak values.

# Returns
- [`ThresholdWidth`](@ref)  
Structure detailing the central threshold and its lower (`x₁`) and upper (`x₂`) bounds.

# Description
Calculates a threshold band (region between two values) by linearly interpolating between the minimum and each peak, weighted by ϵ. This band is used to classify points in the signal for idealization methods.

# Example
```
thr_width = get_threshold_width(hist_analysis, 0.1)
println("Threshold center: ", thr_width.threshold_centre)
println("Lower bound: ", thr_width.x₁)
println("Upper bound: ", thr_width.x₂)
```
"""
function get_threshold_width(hist_analysis::HistPeakAnalysis, ϵ::Float32) :: ThresholdWidth
    max1_val, max2_val, min_val = hist_analysis.edges[hist_analysis.pmax1_index], hist_analysis.edges[hist_analysis.pmax2_index], hist_analysis.edges[hist_analysis.pmin_index]
    point_max1 = point(hist_analysis, :pmax1_index, :pmax1)
	point_min = point(hist_analysis, :pmin_index, :pmin)
	point_max2 = point(hist_analysis, :pmax2_index, :pmax2)
    line_max1 = line(point_max1, point_min) # line from max1 to min
    line_max2 = line(point_max2, point_min) # line from max2 to min
	
	s1 = abs(line_max1.a)
    s2 = abs(line_max2.a)

    # Fallback to symmetric if slopes are unusable
    denom = s1 + s2
    if !isfinite(s1) || !isfinite(s2) || denom == 0
        ε1 = ϵ
        ε2 = ϵ
    else
        # weight gentler side (smaller slope) more
        w1 = s1 / denom
        w2 = s2 / denom
        ε1 = ϵ * Float32(w1)
        ε2 = ϵ * Float32(w2)
    end

    x₁ = (1 - ε1) * min_val + ε1 * max1_val
    x₂ = (1 - ε2) * min_val + ε2 * max2_val
	ThresholdWidth(min_val, x₁, x₂)
end;

"""
    calculate_approximation(data_with_times::Vector{Tuple{Float32, Float32}}, threshold::ThresholdWidth) 
        -> Tuple{Vector{Float32}, Vector{Float32}}

Estimate breakpoints and dwell times in a time-stamped signal using a threshold band.

# Arguments
- `data_with_times::Vector{Tuple{Float32, Float32}}`  
Vector of `(time, value)` pairs for the signal, e.g. output from [`combine_time_with_data`](@ref).
- `threshold::ThresholdWidth`  
Threshold band object defining the central threshold and its lower/upper bounds.

# Returns
- `Tuple{Vector{Float32}, Vector{Float32}}`  
A tuple containing:
    1. `breakpoints::Vector{Float32}` — estimated transition times (seconds)
    2. `dwell_times::Vector{Float32}` — computed dwell times (seconds) between transitions

# Description
This algorithm tracks signal transitions using a threshold band, finding:
- Intervals where the signal is within the threshold band
- Median times of such intervals are treated as breakpoints
- Dwell times are computed as time differences between breakpoints

The method handles both edge and intermediate threshold crossing logic, smoothing via medians if needed, and properly catching transitions even when noise fluctuates.

# Example
```
pairs = combine_time_with_data(data, Δt)
thr_band = get_threshold_width(hist_analysis, ϵ)
breaks, dwell_times = calculate_approximation(pairs, thr_band)
println("First dwell time: ", dwell_times)
```
"""
function calculate_approximation(data_with_times::Vector{Tuple{Float32,Float32}}, threshold::ThresholdWidth)::Tuple{Vector{Float32},Vector{Float32}}
    # accessor functions for point for better readability
    value(point) = point[2]
    time(point) = point[1]

    breakpoints = []
    previous_point = data_with_times[1]
    x1, x2 = threshold.x₁ <= threshold.x₂ ? (threshold.x₁, threshold.x₂) : (threshold.x₂, x₁)
    if value(previous_point) < threshold.threshold_centre
        current_state = 0 # starting at the bottom
    else
        current_state = 1 # starting at the top
    end

    temp_time_list = []
    for point in data_with_times[2:end, :]
        # check if point is in the band
        if x1 < value(point) < x2
            # add that time to the temporary list
            push!(temp_time_list, time(point))
        elseif !isempty(temp_time_list)
            if current_state == 0
                if value(point) > x2
                    # add the breakpoint
                    push!(breakpoints, median(temp_time_list))
                    # change current state
                    current_state = 1
                end
            else
                if value(point) < x1
                    # add the breakpoint
                    push!(breakpoints, median(temp_time_list))
                    # change current state
                    current_state = 0
                end
            end
            # reinitialize the temporary list
            temp_time_list = []
        else
            # naive portion of the algorithm (when ϵ small or 0)
            if current_state == 0
                if value(previous_point) < x1 && value(point) > x2
                    push!(breakpoints, time(point))
                    current_state = 1
                end
            else
                if value(point) < x1 && value(previous_point) > x2
                    push!(breakpoints, time(point))
                    current_state = 0
                end
            end
        end
        # change the previous point to the next one
        previous_point = point
    end
    dwell_times = append!([breakpoints[1]], diff(breakpoints))
    breakpoints, dwell_times
end

"""
    noise_test(noise::Noise) :: Float32

Evaluate the normality of noise samples by batching the data and averaging
Shapiro-Wilk test p-values.

This function splits the noise sequence into fixed-size batches, runs a
`ShapiroWilkTest` on each batch, and returns the mean p-value as a summary
normality score.

# Arguments
- `noise::Noise`: A noise container with fields:
    - `ξ::Vector{Float32}`: The raw noise samples.
    - `μ::Float32`: Mean of the noise (metadata; not used in this function).
    - `σ::Float32`: Standard deviation of the noise (metadata; not used here).

# Returns
- `Float32`: The average p-value from Shapiro-Wilk tests across all complete
  batches. Higher values suggest better agreement with normality.

# Details
- Uses a fixed `batch_size = 50`.
- Only complete batches are tested:
  `num_batches = div(length(noise_data(noise)), batch_size)`.
  Any remainder samples are ignored.
- For each batch, computes `ShapiroWilkTest(batch)` and stores `pvalue(test)`.
- Returns `mean(pvals)`.

# Notes
- Requires `HypothesisTests` (for `ShapiroWilkTest`) and `Statistics` (for `mean`).
- If the total number of samples is less than `batch_size`, the function returns
  `NaN` (since there are no batches/p-values to average). Consider guarding
  against this in calling code.
- The function assumes the presence of `noise_data(noise)`, which should return
  the vector of samples to test (e.g., `noise.ξ`). If not defined, replace
  `noise_data(noise)` with `noise.ξ`.
"""
function noise_test(noise::Noise)::Float32
    # data_1 = rand(Normal(0, 1), 50000)
    batch_size = 50
    num_batches = div(length(noise_data(noise)), batch_size)
    pvals = Float32[]

    for i in 1:num_batches
        batch = noise_data(noise)[(i-1)*batch_size+1:i*batch_size]
        test = ShapiroWilkTest(batch)
        push!(pvals, pvalue(test))
    end

    mean_pval = mean(pvals)
    mean_pval
end

"""
    mika_method(data::Vector{Float32}, Δt::Float32, method::MikaMethod) -> MikaMethodOutput

Apply the Mika idealization algorithm to a signal, using histogram-based thresholding and noise optimization.

# Arguments
- `data::Vector{Float32}`  
The original raw signal to be idealized.
- `Δt::Float32`  
Sampling interval in seconds.
- `method::MikaMethod`  
Parameters for the Mika method, including histogram bin count.

# Returns
- [`MikaMethodOutput`](@ref)  
Comprehensive structure containing breakpoints, dwell times, idealized signal, noise statistics, and optimized thresholds.

# Description
This method:
1. Computes a signal histogram and probability histogram.
2. Analyzes histogram peaks to set up threshold bands.
3. Estimates initial breakpoints and dwell times via the threshold region.
4. Constructs the idealized signal alternating between two peak levels.
5. Assesses noise and fits a normal distribution for MSE optimization.
6. Iteratively adjusts (optimizes) the threshold to minimize noise MSE.
7. Returns all outcomes in a [`MikaMethodOutput`](@ref) struct for further analysis or plotting.

# Example
```
params = MikaMethod(ϵ=0.05, number_of_histogram_bins=80)
result = mika_method(data, Δt, params)
println("Breakpoints: ", breakpoints(result))
println("Noise MSE: ", noise_mse(result))
```
"""
function mika_method(data::Vector{Float32}, Δt::Float32, method::MikaMethod)::MikaMethodOutput
    histogram_of_data = histogram_calculator(data, method.number_of_histogram_bins)
    prob_hist = calculate_probability_histogram(histogram_of_data)
    hist_analysis = analyze_histogram_peaks(prob_hist)
    # calculate initial mse
    data_with_times = combine_time_with_data(data, Δt)

    threshold::ThresholdWidth = get_threshold_width(hist_analysis, Float32(0.0))

    breakpoints, dwell_times_approx = calculate_approximation(data_with_times, threshold)

    idealized_data = idealize_data(data, dwell_times_approx, hist_analysis, Δt)
    best_idealized_data = idealized_data

    # threshold = threshold_width.threshold_centre
    # threshold_index = hist_analysis.pmin_index

    noise_ = noise(data, idealized_data)
    best_noise = noise_
    noise_mse = noise_test(noise_)
    best_noise_mse = noise_mse

    # store initial values to revert back to them
    initial_breakpoints = breakpoints
    initial_dwell_times_approx = dwell_times_approx
    initial_idealized_data = idealized_data
    initial_noise = best_noise
    initial_threshold = threshold
    initial_best_noise_mse = best_noise_mse

    max_left = sum(hist_analysis.weights[hist_analysis.pmax1_index:hist_analysis.pmin_index])
    max_right = sum(hist_analysis.weights[hist_analysis.pmin_index:hist_analysis.pmax2_index])

    step = max_left > max_right ? -1 : 1

    best_threshold = threshold
    best_centre_index = hist_analysis.pmin_index
    # optimize the threshold
    previous_noise_mse = noise_mse
    # @info "Initial noise MSE: $noise_mse"
    for min_ind in hist_analysis.pmin_index+step:step:hist_analysis.pmax1_index
        hist_analysis.pmin_index = min_ind
        threshold = get_threshold_width(hist_analysis, Float32(0.0))

        temp_breakpoints, temp_dwell_times_approx = calculate_approximation(data_with_times, threshold)

        idealized_data = idealize_data(data, temp_dwell_times_approx, hist_analysis, Δt)
        noise_ = noise(data, idealized_data)
        noise_mse = noise_test(noise_)

        if noise_mse > best_noise_mse
            best_centre_index = min_ind
            best_threshold = threshold
            best_noise = noise_
            previous_noise_mse = noise_mse
            best_noise_mse = noise_mse
            best_idealized_data = idealized_data
            breakpoints, dwell_times_approx = temp_breakpoints, temp_dwell_times_approx
        end
    end

    # checking only centre because other values are the same
    if best_threshold.threshold_centre == threshold.threshold_centre
        best_centre_index = hist_analysis.pmin_index
        best_threshold = initial_threshold
        breakpoints = initial_breakpoints
        dwell_times_approx = initial_dwell_times_approx
        best_idealized_data = initial_idealized_data
        best_noise = initial_noise
        best_noise_mse = initial_best_noise_mse
        # return MikaMethodOutput(initial_breakpoints, initial_dwell_times_approx, initial_idealized_data, initial_noise, initial_threshold, initial_threshold_index, initial_best_noise_mse)
    end

    # MikaMethodOutput(breakpoints, dwell_times_approx, best_idealized_data, best_noise, threshold, threshold_index, best_noise_mse)
    hist_analysis.pmin_index = best_centre_index
    # @info "Starting fine-tuning of ϵ"
    # @info "Best noise MSE after histogram edge optimization: $best_noise_mse at threshold centre = $(best_threshold.threshold_centre)"
    for ϵ::Float32 in 0.01:0.01:0.2
        
        threshold = get_threshold_width(hist_analysis, ϵ)
        temp_breakpoints, temp_dwell_times_approx = calculate_approximation(data_with_times, threshold)

        idealized_data = idealize_data(data, temp_dwell_times_approx, hist_analysis, Δt)
        noise_ = noise(data, idealized_data)
        noise_mse = noise_test(noise_)
        # @info "Noise MSE: $noise_mse at ϵ = $ϵ"
        if noise_mse > best_noise_mse
            # @info "New best noise MSE: $noise_mse at ϵ = $ϵ"
            # threshold = hist_analysis.edges[min_ind]
            # threshold_index = min_ind
            best_threshold = threshold
            best_noise = noise_
            previous_noise_mse = noise_mse
            best_noise_mse = noise_mse
            best_idealized_data = idealized_data
            breakpoints, dwell_times_approx = temp_breakpoints, temp_dwell_times_approx
        end
    end
    MikaMethodOutput(breakpoints, dwell_times_approx, best_idealized_data, best_noise, best_threshold, best_noise_mse)
end