using PyCall
"""
    IdealizationMethod

Abstract supertype for all idealization method parameter types.

Each concrete subtype of `IdealizationMethod` stores only the parameters
specific to that method. The associated computation algorithm is defined
separately and linked to the subtype via a `method(::YourMethodType)` function.

# Usage
To implement a new method:
1. Create a new struct subtype of `IdealizationMethod` that stores parameters.
2. Write the algorithm function:

    my_algorithm(data::Vector{Float64}, Δt::Float64, m::MyMethodType) -> Vector{Float64}

3. Link it to the type with:

    method(::MyMethodType) = my_algorithm
"""
abstract type IdealizationMethod end

"""
    MethodOutput

Abstract supertype for all output types produced by idealization methods in IonChannel.

All result types describing the outcome of an idealization algorithm—such as dwell times, breakpoints, noise metrics, and other analysis products—should subtype this.

This enables consistent downstream handling, generic dispatch, and unified documentation for results from different algorithms.
"""
abstract type MethodOutput end

"""
    HistPeakAnalysis

Structure holding the results of peak and trough analysis on a histogram.

# Fields
- `edges::Vector{Float64}`  
The bin edges of the histogram.
- `weights::Vector{Float64}`  
The counts or weights for each histogram bin.
- `pmax1::Float64`  
The height (weight) of the first (primary) maximum peak.
- `pmax1_index::Int`  
The bin index of the first maximum peak.
- `pmax2::Float64`  
The height (weight) of the second maximum peak.
- `pmax2_index::Int`  
The bin index of the second maximum peak.
- `midpoint::Int`  
The index representing the midpoint between the two main peaks.
- `pmin::Float64`  
The height (weight) of the minimum (trough) between the two main peaks.
- `pmin_index::Int`  
The bin index of the minimum (trough) between the two main peaks.

# Description
`HistPeakAnalysis` encapsulates detailed information about the histogram peaks and troughs necessary for further analysis such as threshold determination.

It is typically produced by the function [`analyze_histogram_peaks`](@ref).

# Example
```
hist_analysis = analyze_histogram_peaks(prob_hist)
println("Primary peak at bin ", hist_analysis.pmax1_index, " with weight ", hist_analysis.pmax1)
println("Minimum between peaks at bin ", hist_analysis.pmin_index)
```
"""
mutable struct HistPeakAnalysis
    edges::Vector{Float64}
    weights::Vector{Float64}
    pmax1::Float64
    pmax1_index::Int
    pmax2::Float64
    pmax2_index::Int
    midpoint::Int
    pmin::Float64
    pmin_index::Int
end

"""
    ThresholdWidth

Structure representing the threshold range for an idealization method.

# Fields
- `threshold_centre::Float64`  
The central threshold value used for detecting state transitions.
- `x₁::Float64`  
The lower bound of the threshold band.
- `x₂::Float64`  
The upper bound of the threshold band.

# Description
`ThresholdWidth` encapsulates the main threshold and its bounds within which the signal values are considered for state transition analysis.

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
    threshold_centre::Float64
    x₁::Float64
    x₂::Float64
end

"""
    Point

A simple 2D point with `x` and `y` coordinates.

# Fields
- `x::Float64` — The horizontal coordinate.
- `y::Float64` — The vertical coordinate.

# Example
```
p = Point(1.0, 2.0)
println(p.x) # 1.0
println(p.y) # 2.0
```
"""
struct Point
    x::Float64
    y::Float64
end

"""
    Line

Represents a 2D line defined by the equation `y = a * x + b`.

# Fields
- `a::Float64` — The slope of the line.
- `b::Float64` — The y-intercept of the line.

# Example
```
line = Line(2.0, 1.0) # y = 2x + 1
```
"""
struct Line
    a::Float64
    b::Float64
end

"""
    Noise

Mutable structure holding noise characteristics computed as differences between raw data and idealized values.

# Fields
- `ξ::Vector{Float64}`  
The vector of residuals (noise), computed as `data - idealized_values`.
- `μ::Float64`  
The mean of the noise vector.
- `σ::Float64`  
The standard deviation of the noise vector.

# Description
`Noise` encapsulates statistical information about deviation of the observed data from an idealized signal, useful for assessing fit quality and noise properties.

---

`noise_data(noise::Noise) -> Vector{Float64}`  
Accessor function returning the raw noise vector `ξ`.

`μ(noise::Noise) -> Float64`  
Accessor function returning the mean of the noise.

`σ(noise::Noise) -> Float64`  
Accessor function returning the standard deviation of the noise.

---

`noise(data::Vector{Float64}, idealized_values::Vector{Float64}) -> Noise`  
Computes and returns a `Noise` object representing the noise between raw data and idealized values.

# Arguments
- `data::Vector{Float64}`  
The original data vector.
- `idealized_values::Vector{Float64}`  
The idealized or reconstructed data vector.

# Returns
- `Noise`  
A `Noise` struct containing the residual vector and its statistical properties.

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
    ξ::Vector{Float64}
    μ::Float64
    σ::Float64
end


noise_data(noise) = noise.ξ
μ(noise) = noise.μ
σ(noise) = noise.σ

function noise(data::Vector{Float64}, idealized_values::Vector{Float64}) :: Noise
    ξ = data .- idealized_values
    μ = mean(ξ)
    σ = std(ξ)
    Noise(ξ, μ, σ)
end

"""
    MikaMethodOutput <: MethodOutput

Structure holding the results of the Mika idealization method.

# Fields

- `breakpoints::Vector{Float64}`  
Vector of cumulative breakpoint times (in seconds) indicating detected state changes.

- `dwell_times_approx::Vector{Float64}`  
Estimated dwell times (in seconds) between detected transitions.

- `idealized_data::Vector{Float64}`  
The reconstructed idealized signal corresponding to the original data.

- `noise::Noise`  
A `Noise` struct capturing residuals between the original and idealized data and associated statistics.

- `threshold::Float64`  
The optimized threshold value used for state discrimination.

- `threshold_index::Int64`  
The index within histogram bins corresponding to the selected threshold.

- `noise_mse::Float64`  
The mean squared error between the noise and its fitted normal distribution.

# Description

`MikaMethodOutput` bundles all significant outputs of the Mika idealization procedure, enabling comprehensive downstream analysis and visualization.

# Accessor Functions

Each field has a corresponding accessor function for convenient retrieval:

- `breakpoints(optimized_data::MikaMethodOutput)`  
- `dwell_times_approx(optimized_data::MikaMethodOutput)`  
- `idealized_data(optimized_data::MikaMethodOutput)`  
- `noise(optimized_data::MikaMethodOutput)`  
- `threshold(optimized_data::MikaMethodOutput)`  
- `threshold_index(optimized_data::MikaMethodOutput)`  
- `noise_mse(optimized_data::MikaMethodOutput)`  

These allow accessing parts of the output without field syntax.

# Example

```
result = mika_method(data, Δt, MikaMethod(0.1, 120))

println("Threshold used: ", threshold(result))
println("Number of breakpoints detected: ", length(breakpoints(result)))
```
"""
struct MikaMethodOutput <: MethodOutput
    breakpoints::Vector{Float64}
    dwell_times_approx::Vector{Float64}
    idealized_data::Vector{Float64}
    noise::Noise
    threshold::Float64
    threshold_index::Int64
    noise_mse::Float64
end

breakpoints(optimized_data::MikaMethodOutput) = optimized_data.breakpoints
dwell_times_approx(optimized_data::MikaMethodOutput) = optimized_data.dwell_times_approx
idealized_data(optimized_data::MikaMethodOutput) = optimized_data.idealized_data
noise(optimized_data::MikaMethodOutput) = optimized_data.noise
threshold(optimized_data::MikaMethodOutput) = optimized_data.threshold
threshold_index(optimized_data::MikaMethodOutput) = optimized_data.threshold_index
noise_mse(optimized_data::MikaMethodOutput) = optimized_data.noise_mse  

"""
    MikaMethod <: IdealizationMethod

Parameters and a method function for the Mika idealization method.

# Fields

- `ϵ::Float64`  
A smoothing/noise parameter influencing threshold optimization.

- `number_of_histogram_bins::UInt8`  
The number of bins used in histogram computations during idealization.

# Description

`MikaMethod` stores parameters essential for controlling the thresholding and histogram binning in the Mika method. The associated algorithm is linked via `method_function(::MikaMethod)`.

# Example
```
m = MikaMethod(0.1, 100)
result = calculate_method(data, m, Δt)
```
"""
mutable struct MikaMethod <: IdealizationMethod
    ϵ::Float64
    number_of_histogram_bins::UInt8
end

"""
    MeanDeviationMethod <: IdealizationMethod

Parameters for the *deviation-from-running-mean* idealization method.

# Fields
- `δ::Float64` - Deviation offset subtracted from the absolute deviation before thresholding.
- `λ::Float64` - Threshold value above which a deviation from the mean indicates a state change.

# Description
`MeanDeviationMethod` stores only these two numeric parameters.
The actual algorithm is implemented in [`deviation_from_mean_method`](@ref)
and linked via:

    method(::MeanDeviationMethod) = deviation_from_mean_method

# Example
```m = MeanDeviationMethod(0.05, 0.5)

println(method_function(m)) # shows the stored function
println(δ(m)) # 0.05
println(λ(m)) # 0.5
```
"""
mutable struct MeanDeviationMethod <: IdealizationMethod
    δ::Float64
    λ::Float64
end

"""
    δ(m::MeanDeviationMethod) -> Float64

Return the `δ` (delta) parameter from the given `MeanDeviationMethod` instance.
"""
δ(m::MeanDeviationMethod) = m.δ

"""
    λ(m::MeanDeviationMethod) -> Float64

Return the `λ` (lambda) parameter from the given `MeanDeviationMethod` instance.
"""
λ(m::MeanDeviationMethod) = m.λ

"""
    MeanDeviationMethodOutput <: MethodOutput

Type describing the output of the mean deviation idealization method.

# Fields

- `breakpoints::Vector{Float64}`
    The cumulative breakpoints (timepoints, in seconds) at which state changes occur in the idealized trace, as computed by the method.

- `dwell_times_approx::Vector{Float64}`
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
    breakpoints::Vector{Float64}
    dwell_times_approx::Vector{Float64}
    idealized_data::Vector{Int8}
end

struct DeepChannelMethod <: IdealizationMethod
    model::PyObject
end

struct DeepChannelMethodOutput <: MethodOutput
    dwell_times_approx::Vector{Float64}
    breakpoints::Vector{Float64}
    idealized_data::Vector{Int8}
end

struct NaiveMethod <: IdealizationMethod
    number_of_histogram_bins::UInt16
end

struct NaiveMethodOutput <: MethodOutput
    dwell_times_approx::Vector{Float64}
    breakpoints::Vector{Float64}
    idealized_data::Vector{Int8}
end

struct MeanError
    mean_squared_error::Float64
    mean_accuracy::Float64
    mean_squared_errors::Vector{Float64}
    accuracies::Vector{Float64}
end

struct MDLMethod <: IdealizationMethod
    min_seg::UInt16
    threshold::Float64
    number_of_histogram_bins::UInt16
end

struct MDLMethodOutput <: MethodOutput
    breakpoints::Vector{Float64}
    dwell_times_approx::Vector{Float64}
    idealized_data::Vector{Int8}
end