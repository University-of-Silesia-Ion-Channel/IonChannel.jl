using PyCall
"""
    IdealizationMethod

Abstract supertype for all idealization method parameter types.

Each concrete subtype of [`IdealizationMethod`](@ref) stores only the parameters
specific to that method. The associated computation algorithm is defined
separately and linked to the subtype via a `method(::YourMethodType)` function.

# Usage
To implement a new method:
1. Create a new struct subtype of [`IdealizationMethod`](@ref) that stores parameters.
2. Write the algorithm function:

    my_algorithm(data::Vector{Float32}, Δt::Float32, m::MyMethodType) -> Vector{Float32}

3. Link it to the type with:

    method(::MyMethodType) = my_algorithm
"""
abstract type IdealizationMethod end

"""
    MethodOutput

Abstract supertype for all output types produced by idealization methods in [`IonChannel`](@ref).

All result types describing the outcome of an idealization algorithm—such as dwell times, breakpoints, noise metrics, and other analysis products—should subtype this.

This enables consistent downstream handling, generic dispatch, and unified documentation for results from different algorithms.
"""
abstract type MethodOutput end

"""
    HistPeakAnalysis

Structure holding the results of peak and trough analysis on a histogram.

# Fields
- `edges::Vector{Float32}`  
The bin edges of the histogram.
- `weights::Vector{Float32}`  
The counts or weights for each histogram bin.
- `pmax1::Float32`  
The height (weight) of the first (primary) maximum peak.
- `pmax1_index::Int`  
The bin index of the first maximum peak.
- `pmax2::Float32`  
The height (weight) of the second maximum peak.
- `pmax2_index::Int`  
The bin index of the second maximum peak.
- `midpoint::Int`  
The index representing the midpoint between the two main peaks.
- `pmin::Float32`  
The height (weight) of the minimum (trough) between the two main peaks.
- `pmin_index::Int`  
The bin index of the minimum (trough) between the two main peaks.

# Description
[`HistPeakAnalysis`](@ref) encapsulates detailed information about the histogram peaks and troughs necessary for further analysis such as threshold determination.

It is typically produced by the function [`analyze_histogram_peaks`](@ref).

# Example
```
hist_analysis = analyze_histogram_peaks(prob_hist)
println("Primary peak at bin ", hist_analysis.pmax1_index, " with weight ", hist_analysis.pmax1)
println("Minimum between peaks at bin ", hist_analysis.pmin_index)
```
"""
mutable struct HistPeakAnalysis
    edges::Vector{Float32}
    weights::Vector{Float32}
    pmax1::Float32
    pmax1_index::Int
    pmax2::Float32
    pmax2_index::Int
    midpoint::Int
    pmin::Float32
    pmin_index::Int
end

"""
    ThresholdWidth

Structure representing the threshold range for an idealization method.

# Fields
- `threshold_centre::Float32`  
The central threshold value used for detecting state transitions.
- `x₁::Float32`  
The lower bound of the threshold band.
- `x₂::Float32`  
The upper bound of the threshold band.

# Description
[`ThresholdWidth`](@ref) encapsulates the main threshold and its bounds within which the signal values are considered for state transition analysis.

It is typically computed during histogram-based threshold optimization and used for breakpoint detection.

# Example
```
thr_width = ThresholdWidth(0.5, 0.3, 0.7)
println("Center threshold: ", thr_width.threshold_centre)
println("Lower bound: ", thr_width.x₁)
println("Upper bound: ", thr_width.x₂)
```
"""
struct ThresholdWidth
    threshold_centre::Float32
    x₁::Float32
    x₂::Float32
end

"""
    Point

A simple 2D point with `x` and `y` coordinates.

# Fields
- `x::Float32` — The horizontal coordinate.
- `y::Float32` — The vertical coordinate.

# Example
```
p = Point(1.0, 2.0)
println(p.x) # 1.0
println(p.y) # 2.0
```
"""
struct Point
    x::Float32
    y::Float32
end

"""
    Line

Represents a 2D line defined by the equation `y = a * x + b`.

# Fields
- `a::Float32` — The slope of the line.
- `b::Float32` — The y-intercept of the line.

# Example
```
line = Line(2.0, 1.0) # y = 2x + 1
```
"""
struct Line
    a::Float32
    b::Float32
end

"""
    Noise

Mutable structure holding noise characteristics computed as differences between raw data and idealized values.

# Fields
- `ξ::Vector{Float32}`  
The vector of residuals (noise), computed as `data - idealized_values`.
- `μ::Float32`  
The mean of the noise vector.
- `σ::Float32`  
The standard deviation of the noise vector.

# Description
[`Noise`](@ref) encapsulates statistical information about deviation of the observed data from an idealized signal, useful for assessing fit quality and noise properties.

---

`noise_data(noise::Noise) -> Vector{Float32}`  
Accessor function returning the raw noise vector `ξ`.

`μ(noise::Noise) -> Float32`  
Accessor function returning the mean of the noise.

`σ(noise::Noise) -> Float32`  
Accessor function returning the standard deviation of the noise.

---

`noise(data::Vector{Float32}, idealized_values::Vector{Float32}) -> Noise`  
Computes and returns a [`Noise`](@ref) object representing the noise between raw data and idealized values.

# Arguments
- `data::Vector{Float32}`  
The original data vector.
- `idealized_values::Vector{Float32}`  
The idealized or reconstructed data vector.

# Returns
- [`Noise`](@ref)  
A [`Noise`](@ref) struct containing the residual vector and its statistical properties.

# Example
```
data = [1.0, 2.1, 2.9, 4.1]
idealized = [1.0, 2.0, 3.0, 4.0]
n = noise(data, idealized)
println("Noise mean: ", μ(n))
println("Noise std dev: ", σ(n))
```
"""
mutable struct Noise
    ξ::Vector{Float32}
    μ::Float32
    σ::Float32
end

"""
    noise_data(noise::Noise) -> Vector{Float32}

Accessor function returning the raw noise vector `ξ`.
"""
noise_data(noise::Noise) = noise.ξ
"""
    μ(noise::Noise) -> Float32

Accessor function returning the mean of the noise.
"""    
μ(noise::Noise) = noise.μ

"""
    σ(noise::Noise) -> Float32
    
Accessor function returning the standard deviation of the noise.
"""
σ(noise::Noise) = noise.σ

"""
Computes and returns a [`Noise`](@ref) object representing the noise between raw data and idealized values.
See more in the [`Noise`](@ref) struct documentation.
"""
function noise(data::Vector{Float32}, idealized_values::Vector{Float32}) :: Noise
    ξ = data .- idealized_values
    μ = mean(ξ)
    σ = std(ξ)
    Noise(ξ, μ, σ)
end

"""
    MikaMethodOutput <: MethodOutput

Structure holding the results of the Mika idealization method.

# Fields

- `breakpoints::Vector{Float32}`  
Vector of cumulative breakpoint times (in seconds) indicating detected state changes.

- `dwell_times_approx::Vector{Float32}`  
Estimated dwell times (in seconds) between detected transitions.

- `idealized_data::Vector{Float32}`  
The reconstructed idealized signal corresponding to the original data.

- `noise::Noise`  
A [`Noise`](@ref) struct capturing residuals between the original and idealized data and associated statistics.

- `threshold::ThresholdWidth`  
The threshold and its bounds used for state transition detection.

- `noise_mse::Float32`  
The mean squared error between the noise and its fitted normal distribution.

# Description

[`MikaMethodOutput`](@ref) bundles all significant outputs of the Mika idealization procedure, enabling comprehensive downstream analysis and visualization.

# Accessor Functions

Each field has a corresponding accessor function for convenient retrieval:

- `breakpoints(optimized_data::MikaMethodOutput)`  
- `dwell_times_approx(optimized_data::MikaMethodOutput)`  
- `idealized_data(optimized_data::MikaMethodOutput)`  
- `noise(optimized_data::MikaMethodOutput)`  
- `threshold(optimized_data::MikaMethodOutput)`  
- `noise_mse(optimized_data::MikaMethodOutput)`  

These allow accessing parts of the output without field syntax.

# Example

```
result = mika_method(data, Δt, MikaMethod(120))

println("Threshold used: ", threshold(result))
println("Number of breakpoints detected: ", length(breakpoints(result)))
```
"""
struct MikaMethodOutput <: MethodOutput
    breakpoints::Vector{Float32}
    dwell_times_approx::Vector{Float32}
    idealized_data::Vector{Float32}
    noise::Noise
    threshold::ThresholdWidth
    noise_mse::Float32
end

breakpoints(optimized_data::MikaMethodOutput) = optimized_data.breakpoints
dwell_times_approx(optimized_data::MikaMethodOutput) = optimized_data.dwell_times_approx
idealized_data(optimized_data::MikaMethodOutput) = optimized_data.idealized_data
noise(optimized_data::MikaMethodOutput) = optimized_data.noise
threshold(optimized_data::MikaMethodOutput) = optimized_data.threshold
noise_mse(optimized_data::MikaMethodOutput) = optimized_data.noise_mse  

"""
    MikaMethod <: IdealizationMethod

Parameters and a method function for the Mika idealization method.

# Fields

- `number_of_histogram_bins::UInt16`  
The number of bins used in histogram computations during idealization.

# Description

[`MikaMethod`](@ref) stores parameters essential for controlling the thresholding and histogram binning in the Mika method. The associated algorithm is linked via `method_function(::MikaMethod)`.

# Example
```
m = MikaMethod(100)
result = calculate_method(data, m, Δt)
```
"""
mutable struct MikaMethod <: IdealizationMethod
    number_of_histogram_bins::UInt16
end

"""
    MeanDeviationMethod <: IdealizationMethod

Parameters for the *deviation-from-running-mean* idealization method.

# Fields
- `δ::Float32` - Deviation offset subtracted from the absolute deviation before thresholding.

# Description
[`MeanDeviationMethod`](@ref) stores only these two numeric parameters.
The actual algorithm is implemented in [`deviation_from_mean_method`](@ref)
and linked via:

    method(::MeanDeviationMethod) = deviation_from_mean_method

# Example
```m = MeanDeviationMethod(Float32(0.0))

println(method_function(m)) # shows the stored function
println(δ(m)) # 0.00
```
"""
mutable struct MeanDeviationMethod <: IdealizationMethod
    δ::Float32
end

"""
    δ(m::MeanDeviationMethod) -> Float32

Return the `δ` (delta) parameter from the given [`MeanDeviationMethod`](@ref) instance.
"""
δ(m::MeanDeviationMethod) = m.δ

"""
    MeanDeviationMethodOutput <: MethodOutput

Type describing the output of the mean deviation idealization method.

# Fields

- `breakpoints::Vector{Float32}`
    The cumulative breakpoints (timepoints, in seconds) at which state changes occur in the idealized trace, as computed by the method.

- `dwell_times_approx::Vector{Float32}`
    The sequence of estimated dwell times between state transitions (in seconds) produced by the algorithm.

# Description

This struct bundles the main outputs from the mean deviation idealization routine.
It allows downstream analysis, comparison, or visualization of both raw and processed results.

Typically constructed and returned by [`deviation_from_mean_method`](@ref) or generic idealization runners.

# Example

```
result = deviation_from_mean_method(normalized_data, Δt, MeanDeviationMethod(δ, λ))
result.breakpoints # Vector of breakpoints (seconds)
result.dwell_times_approx # Vector of dwell times (seconds)
```  
"""
struct MeanDeviationMethodOutput <: MethodOutput
    breakpoints::Vector{Float32}
    dwell_times_approx::Vector{Float32}
    idealized_data::Vector{UInt8}
end


"""
    DeepChannelMethod <: IdealizationMethod

Wrapper for a Python-based deep learning model used to idealize ion-channel
time-series data into discrete states.

This type stores a reference to a Python model (e.g., a Keras/PyTorch model
accessed via `PyCall.PyObject`) that supports per-sample classification
(prediction) over the input trace after appropriate preprocessing.

# Fields
- `model::PyObject`: The underlying Python model object. It is expected to
  expose a prediction API compatible with the calling code (e.g., a `.predict`
  method that accepts a 4D tensor shaped like `(N, 1, 1, 1)` and returns
  class probabilities per sample).

# Usage
- Typically passed into `deep_channel_method(data, Δt, c_method)` to produce:
  - per-sample predicted class labels,
  - transition breakpoints,
  - approximate dwell times.

# Notes
- Ensure PyCall is properly initialized and the Python environment includes all
  dependencies required by the model (e.g., TensorFlow/Keras or PyTorch).
- Input preprocessing (e.g., scaling with `UnitRangeTransform` and reshaping to
  `(N,1,1,1)`) is handled by [`deep_channel_method`](@ref).
- The model is expected to output a probability distribution over states for
  each sample; downstream code uses `argmax` to select the most likely class.
"""
struct DeepChannelMethod <: IdealizationMethod
    model::PyObject
end


"""
    DeepChannelMethodOutput <: MethodOutput

Container for the outputs of the deep learning–based idealization (DeepChannel method).

Holds the approximate dwell times, transition times (breakpoints), and the full
per-sample idealized state sequence.

# Fields
- `dwell_times_approx::Vector{Float32}`: Approximate durations spent in each
  detected state segment, typically computed from successive transition times.
- `breakpoints::Vector{Float32}`: Transition times (relative to the start of
  the trace) at which the state changes.
- `idealized_data::Vector{UInt8}`: Per-sample predicted states (e.g., `0`/`1`)
  produced by the model after preprocessing and decoding.

# Notes
- Units of `dwell_times_approx` and `breakpoints` depend on the sampling
  interval used upstream (e.g., `Δt` in seconds or milliseconds).
- `idealized_data` provides the full-resolution state assignment, while
  `dwell_times_approx` and `breakpoints` summarize state changes.
- Designed to be returned by [`deep_channel_method`](@ref).
"""
struct DeepChannelMethodOutput <: MethodOutput
    dwell_times_approx::Vector{Float32}
    breakpoints::Vector{Float32}
    idealized_data::Vector{UInt8}
end

"""
    NaiveMethod <: IdealizationMethod

Configuration for a simple histogram-threshold-based idealization method.

This method estimates a decision threshold from the empirical distribution of
signal amplitudes (via histogram analysis) and detects state transitions when
the signal crosses that threshold.

# Fields
- `number_of_histogram_bins::UInt16`: Number of bins to use when building the
  amplitude histogram. A larger value can capture finer structure but may be
  noisier; typical ranges are 50-200 depending on data length and noise.

# Usage
- Pass an instance to `naive_method(data, Δt, method)` to produce:
  - transition `breakpoints`,
  - `dwell_times`,
  - per-sample binary `idealized_data` (0/1).

# Notes
- The chosen bin count influences the stability of the estimated threshold.
- Works best for bimodal amplitude distributions with reasonably separated modes.
"""
struct NaiveMethod <: IdealizationMethod
    number_of_histogram_bins::UInt16
end

"""
    NaiveMethodOutput <: MethodOutput

Container for outputs produced by the histogram-threshold-based idealization
pipeline ([`naive_method`](@ref)).

# Fields
- `dwell_times_approx::Vector{Float32}`: Approximate durations spent in each
  state segment, typically computed from differences between successive
  breakpoints (first element usually equals the first breakpoint time).
- `breakpoints::Vector{Float32}`: Transition times (relative to the start of
  the trace) at which the signal crosses the estimated threshold and the state
  flips.
- `idealized_data::Vector{UInt8}`: Per-sample binary state assignments (`0` or `1`)
  derived by tracking threshold crossings over the time series.

# Notes
- Units of `dwell_times_approx` and `breakpoints` are determined by the sampling
  interval used upstream (e.g., `Δt`).
- `idealized_data` provides the full-length state sequence, while
  `dwell_times_approx` and `breakpoints` summarize state changes.
- Designed to be returned by [`naive_method`](@ref).
"""
struct NaiveMethodOutput <: MethodOutput
    dwell_times_approx::Vector{Float32}
    breakpoints::Vector{Float32}
    idealized_data::Vector{UInt8}
end


"""
    MeanError

Aggregate metrics for evaluating idealization or prediction performance across
multiple segments, traces, or runs.

# Fields
- `mean_squared_error::Float32`: The overall mean of squared errors, typically
  averaged across all samples and/or batches.
- `mean_accuracy::Float32`: The overall accuracy aggregated across evaluations,
  commonly expressed as a fraction in `[0.0, 1.0]`.
- `mean_squared_errors::Vector{Float32}`: Per-segment or per-batch MSE values
  used to compute `mean_squared_error`.
- `accuracies::Vector{Float32}`: Per-segment or per-batch accuracy values used
  to compute `mean_accuracy`.

# Notes
- This struct is a convenient container when running repeated evaluations (e.g.,
  cross-validation folds, multiple traces, or bootstrapped subsets).
- `mean_squared_error` is typically computed from `mean(mean_squared_errors)`,
  and `mean_accuracy` from `mean(accuracies)`, but storing all components keeps
  downstream analysis flexible (e.g., computing variance or confidence intervals).
"""
struct MeanError
    mean_squared_error::Float32
    mean_accuracy::Float32
    mean_squared_errors::Vector{Float32}
    accuracies::Vector{Float32}
end

"""
    MDLMethod <: IdealizationMethod

Configuration for an idealization approach based on the Minimum Description Length (MDL) principle.

The MDL method selects a piecewise-constant segmentation of the signal that minimizes
the total description length (model complexity + data fit). It balances the number of
segments (penalizing over-segmentation) against how well segments explain the data.

# Fields
- `min_seg::UInt16`: Minimum allowed segment length (in samples). Prevents overly short
  segments that are likely due to noise.
- `threshold::Float32`: Amplitude or penalty threshold used within the MDL search/criteria.
  Its interpretation depends on the specific implementation (e.g., penalty weight, merge/split
  decision threshold).
- `number_of_histogram_bins::UInt16`: Number of bins for any histogram-based auxiliary steps
  (e.g., estimating noise statistics or amplitude modes) used by the MDL routine.

# Usage
- Construct and pass to an MDL-based idealization function, e.g.:
  `mdl_method(data, Δt, method::MDLMethod) -> MDLMethodOutput` (function and output
  type names may vary based on your implementation).

# Notes
- The MDL criterion typically trades off data fidelity and model complexity;
  tuning `min_seg` and `threshold` impacts the balance between false positives
  (spurious segments) and missed transitions.
- Works well when the signal is piecewise constant with sparse change points and
  approximately stationary noise.
"""
struct MDLMethod <: IdealizationMethod
    min_seg::UInt16
    threshold::Float32
    number_of_histogram_bins::UInt16
end


"""
    MDLMethodOutput <: MethodOutput

Container for the results of an MDL-based idealization.

Holds the detected change-point times (breakpoints), the corresponding
approximate dwell durations, and the full per-sample idealized state
sequence.

# Fields
- `breakpoints::Vector{Float32}`: Times (relative to the start of the trace)
  at which the signal is estimated to change state according to the MDL
  segmentation.
- `dwell_times_approx::Vector{Float32}`: Approximate durations spent in each
  state segment, typically computed from successive `breakpoints` (the first
  dwell time usually equals the first breakpoint time).
- `idealized_data::Vector{UInt8}`: Per-sample state assignments (e.g., `0`/`1`)
  produced by the MDL idealization over the entire trace length.

# Notes
- Units of `breakpoints` and `dwell_times_approx` depend on the sampling
  interval used upstream (e.g., `Δt` in seconds or milliseconds).
- `idealized_data` provides the full-resolution labeling, while
  `breakpoints` and `dwell_times_approx` summarize the segmentation.
- Designed to be returned by an MDL idealization routine (e.g., [`mdl_method`](@ref)).
"""
struct MDLMethodOutput <: MethodOutput
    breakpoints::Vector{Float32}
    dwell_times_approx::Vector{Float32}
    idealized_data::Vector{UInt8}
end