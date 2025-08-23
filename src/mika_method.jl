using Statistics
using Distributions
using StatsBase
using HypothesisTests

"""
    method_function(::MikaMethod) -> Function

Return the algorithm function associated with a `MikaMethod`.

This allows code like `calculate_method(data, m, Δt)` to work for any
`IdealizationMethod` subtype without changing the executor logic.
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
- `Point`  
A `Point` instance where `x` is taken from `hist.edges` at the specified index,
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
function point(hist::HistPeakAnalysis, indexfield::Symbol, valuefield::Symbol) :: Point
    x = hist.edges[getfield(hist, indexfield)]
    y = getfield(hist, valuefield)
    Point(x, y)
end

"""
    line(point1::Point, point2::Point) -> Line

Construct a `Line` (y = a*x + b) passing through two points.

# Arguments
- `point1::Point`  
The first point `(x₁, y₁)` through which the line will pass.
- `point2::Point`  
The second point `(x₂, y₂)` through which the line will pass.

# Returns
- `Line`  
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
function line(point1::Point, point2::Point) :: Line
    a = (point1.y - point2.y)/(point1.x - point2.x)
    b = point1.y - a * point1.x
    Line(a, b)
end

"""
    calculate_approximation(data_with_times::Vector{Tuple{Float32, Float32}}, threshold::ThresholdWidth) 
        -> Tuple{Vector{Float32}, Vector{Float32}}

Estimate breakpoints and dwell times in a time-stamped signal using a threshold band.

# Arguments
- `data_with_times::Vector{Tuple{Float32, Float32}}`  
Vector of `(time, value)` pairs for the signal, e.g. output from `combine_time_with_data`.
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
function calculate_approximation(data_with_times::Vector{Tuple{Float32, Float32}}, threshold::ThresholdWidth) :: Tuple{Vector{Float32}, Vector{Float32}}
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
    fit_normal_to_noise(noise::Noise) -> Tuple{Histogram, Vector{Float32}}

Fit a normal distribution to noise data and compute its PDF alongside the normalized noise histogram.

# Arguments
- `noise::Noise`  
A `Noise` struct containing residuals and their statistical parameters (mean and standard deviation).

# Returns
- `Tuple{Histogram, Vector{Float32}}`  
A tuple containing:
    1. A normalized histogram of the noise data representing the empirical PDF.
    2. A vector of PDF values from the fitted normal distribution evaluated at the histogram bin edges.

# Description
This function fits a Normal distribution model (`Normal(μ, σ)`) to the noise residuals,
then produces a normalized histogram (PDF mode) of the noise data for comparison.
It also computes theoretical PDF values of the fitted distribution at corresponding bin edges.

# Example
```
noise_obj = noise(data, idealized_values)
hist_pdf, fitted_pdf = fit_normal_to_noise(noise_obj)
plot(hist_pdf, label="Noise histogram PDF")
plot!(hist_pdf.edges, fitted_pdf, label="Fitted normal PDF")
```
"""
function fit_normal_to_noise(noise::Noise) :: Tuple{Histogram, Vector{Float32}}
    fitted_dist = Normal(μ(noise), σ(noise))
    noise_histogram_pdf = fit(Histogram, noise_data(noise), nbins=100)
    noise_histogram_pdf = normalize(noise_histogram_pdf, mode=:pdf)
    fitted_dist_pdf = pdf.(fitted_dist, noise_histogram_pdf.edges)[1]
    noise_histogram_pdf, fitted_dist_pdf
end

"""
    fit_mse(noise::Histogram, fit::Vector{Float32}) -> Float32

Calculate the mean squared error (MSE) between the histogram bin weights and a fitted model vector.

# Arguments
- `noise::Histogram`  
Histogram representing the observed noise distribution (bin weights).
- `fit::Vector{Float32}`  
Vector of fitted values corresponding to PDF or model estimates evaluated at histogram bins.

# Returns
- `Float32`  
The mean squared error computed as the average squared difference between the histogram weights and fitted values (excluding the last fitted value).

# Description
Computes a goodness-of-fit metric quantifying how closely the fitted vector approximates the observed histogram weights.

# Example
```
mse = fit_mse(histogram_noise, fitted_pdf_vector)
println("Mean squared error: ", mse)
```
"""
fit_mse(noise::Histogram, fit::Vector{Float32}) :: Float32 = sum((noise.weights .- fit[1:end-1]) .^ 2) / length(noise.weights)


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
function noise_test(noise::Noise) :: Float32
	# data_1 = rand(Normal(0, 1), 50000)
	batch_size = 50
	num_batches = div(length(noise_data(noise)), batch_size)
	pvals = Float32[]
	
	for i in 1:num_batches
	    batch = noise_data(noise)[(i-1)*batch_size+1 : i*batch_size]
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
Parameters for the Mika method, including histogram bin count and threshold mixing factor.

# Returns
- `MikaMethodOutput`  
Comprehensive structure containing breakpoints, dwell times, idealized signal, noise statistics, and optimized thresholds.

# Description
This method:
1. Computes a signal histogram and probability histogram.
2. Analyzes histogram peaks to set up threshold bands.
3. Estimates initial breakpoints and dwell times via the threshold region.
4. Constructs the idealized signal alternating between two peak levels.
5. Assesses noise and fits a normal distribution for MSE optimization.
6. Iteratively adjusts (optimizes) the threshold to minimize noise MSE.
7. Returns all outcomes in a `MikaMethodOutput` struct for further analysis or plotting.

# Example
```
params = MikaMethod(ϵ=0.05, number_of_histogram_bins=80)
result = mika_method(data, Δt, params)
println("Breakpoints: ", breakpoints(result))
println("Noise MSE: ", noise_mse(result))
```
"""
function mika_method(data::Vector{Float32}, Δt::Float32, method::MikaMethod) :: MikaMethodOutput
    histogram_of_data = histogram_calculator(data, method.number_of_histogram_bins)
    prob_hist = calculate_probability_histogram(histogram_of_data)
    hist_analysis = analyze_histogram_peaks(prob_hist)
    # calculate initial mse
    data_with_times = combine_time_with_data(data, Δt)

    threshold_width = get_threshold_width(hist_analysis, method.ϵ)

    breakpoints, dwell_times_approx = calculate_approximation(data_with_times, threshold_width)

    idealized_data = idealize_data(data, dwell_times_approx, hist_analysis, Δt)
    best_idealized_data = idealized_data

    threshold = threshold_width.threshold_centre
    threshold_index = hist_analysis.pmin_index

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
    initial_threshold_index = threshold_index
    initial_best_noise_mse = best_noise_mse

    max_left = sum(hist_analysis.weights[hist_analysis.pmax1_index:hist_analysis.pmin_index])
    max_right = sum(hist_analysis.weights[hist_analysis.pmin_index:hist_analysis.pmax2_index])
    
    step = max_left > max_right ? -1 : 1
    
    # optimize the threshold
    previous_noise_mse = noise_mse
    # @info "Initial noise MSE: $noise_mse"
    for min_ind in hist_analysis.pmin_index+step:step:hist_analysis.pmax1_index
        hist_analysis.pmin_index = min_ind
        threshold_width = get_threshold_width(hist_analysis, method.ϵ)
        
        temp_breakpoints, temp_dwell_times_approx = calculate_approximation(data_with_times, threshold_width)
        
        idealized_data = idealize_data(data, temp_dwell_times_approx, hist_analysis, Δt)
        noise_ = noise(data, idealized_data)
        noise_mse = noise_test(noise_)
        
        if noise_mse > best_noise_mse
            threshold = hist_analysis.edges[min_ind]
            threshold_index = min_ind
            best_noise = noise_
            previous_noise_mse = noise_mse
            best_noise_mse = noise_mse
            best_idealized_data = idealized_data
            breakpoints, dwell_times_approx = temp_breakpoints, temp_dwell_times_approx
        end
    end

    if threshold_index == hist_analysis.pmax1_index
        return MikaMethodOutput(initial_breakpoints, initial_dwell_times_approx, initial_idealized_data, initial_noise, initial_threshold, initial_threshold_index, initial_best_noise_mse)
    end

    MikaMethodOutput(breakpoints, dwell_times_approx, best_idealized_data, best_noise, threshold, threshold_index, best_noise_mse)
end