"""
    module IonChannel

A Julia module for advanced idealization and analysis of ion channel single-molecule data.

# Overview

`IonChannel` provides type-safe, highly modular tools for the idealization, segmentation, and analysis of ion channel current traces. It implements and integrates methods based on running-mean deviation, histogram-based segmentation, threshold optimization, MDL change-point detection, and deep-learning-based state detection, supporting workflows in quantitative biophysics and single-molecule electrophysiology.

## Main Features

- Idealization Methods:
  - Abstract parameter types (`IdealizationMethod`) and result types (`MethodOutput`) for extensible algorithms.
  - Implementations: deviation-from-mean (`MeanDeviationMethod`), Mika method (`MikaMethod`), Naive histogram thresholding (`NaiveMethod`), MDL segmentation (`MDLMethod`), and DeepChannel (Python model via `PyCall`, `DeepChannelMethod`).

- Data Handling:
  - Read paired raw traces and dwell times from files: `read_data`, or enumerate dataset folders with `read_all_file_paths`.
  - Extract subsets, align dwell windows, and normalize signals: `get_specified_datapoints`, `normalize_data`, `combine_time_with_data`.

- Histogram Analysis:
  - Build histograms and probability histograms: `histogram_calculator`, `calculate_probability_histogram`.
  - Peak/trough and midpoint detection for bimodal signals: `analyze_histogram_peaks` → `HistPeakAnalysis`.

- Noise and Fit Evaluation:
  - Residual modeling and statistics with `Noise` and `noise_data`.
  - Dwell-time distribution comparison via `calculate_mean_square_error` (returns MSE and histograms).
  - End-to-end benchmark across datasets with `mean_error` (returns `MeanError`).

- Visualization:
  - Overlay exact vs. approximated breakpoints and thresholds: `show_approx_on_plot`.
  - Stacked plots of raw vs. idealized sequences: `plot_idealization_representation`.
  - Threshold and peaks on histogram (Mika): `show_threshold_on_plot`.

## Typical Workflow

1. Reading Data:
   Use `read_data` or `read_all_file_paths` to load trace data and associated dwell times.

2. Idealization:
   Choose an idealization method struct (e.g. `MeanDeviationMethod`, `MikaMethod`, `NaiveMethod`, `MDLMethod`, `DeepChannelMethod`) and run the algorithm (e.g., `calculate_method`, `mika_method`, `naive_method`, `mdl_method`, `deep_channel_method`) to obtain a `<: MethodOutput` with `breakpoints`, `dwell_times_approx`, and an `idealized_data` sequence.

3. Histogram and Threshold:
   Generate and analyze histograms with `histogram_calculator` and `analyze_histogram_peaks`. Convert counts to probabilities via `calculate_probability_histogram`.

4. Visualization:
   Visualize results using `show_approx_on_plot`, `show_threshold_on_plot`, and `plot_idealization_representation`.

5. Noise and Evaluation:
   Quantify residuals with `noise_data` and evaluate dwell-time fit via `calculate_mean_square_error`; assess full-run performance with `mean_error`.

## Extending the Module

To add a new idealization method:
1. Create a new parameter struct `<: IdealizationMethod` (store parameters only).
2. Implement your algorithm `my_method(data::Vector{Float32}, t::Float32, m::MyMethod) -> <: MethodOutput`.
3. Link your method via `method_function(::MyMethod) = my_method` so it works with the generic `calculate_method`.

## Key Types

- Abstractions: `IdealizationMethod`, `MethodOutput`
- Methods and outputs: `MeanDeviationMethod`, `MeanDeviationMethodOutput`; `MikaMethod`, `MikaMethodOutput`; `NaiveMethod`, `NaiveMethodOutput`; `MDLMethod`, `MDLMethodOutput`; `DeepChannelMethod`, `DeepChannelMethodOutput`
- Analysis helpers: `HistPeakAnalysis`, `ThresholdWidth`, `Noise`, `Point`, `Line`
- Evaluation aggregates: `MeanError`

## Selected API (canonical names)

- Data I/O and prep:
  - `read_data(data_file_path::String, dwell_times_path::String) -> (x::Vector{Float32}, y::Vector{Float32})`
  - `read_all_file_paths(data_folder::String) -> (what_first_path::String, data_file_paths::Vector{String}, dwell_times_file_paths::Vector{String})`
  - `get_specified_datapoints(x::Vector{Float32}, y::Vector{Float32}, t::Float32, data_size::UInt32=0) -> Dict("x"=>..., "dwell times"=>...)`
  - `normalize_data(data::Dict{String,Vector{Float32}}) -> Vector{Float32}`
  - `combine_time_with_data(x::Vector{Float32}, t::Float32; batch_size::UInt8=1) -> Vector{Tuple{Float32,Float32}}`

- Histogram and peaks:
  - `histogram_calculator(x::Vector{Float32}, bins::UInt16=100) -> StatsBase.Histogram`
  - `calculate_probability_histogram(h::Histogram) -> Histogram`
  - `analyze_histogram_peaks(h::Histogram) -> HistPeakAnalysis`

- Idealization methods:
  - `deviation_from_mean_method(x::Vector{Float32}, t::Float32, m::MeanDeviationMethod) -> MeanDeviationMethodOutput`
  - `mika_method(x::Vector{Float32}, t::Float32, m::MikaMethod) -> MikaMethodOutput`
  - `naive_method(x::Vector{Float32}, t::Float32, m::NaiveMethod) -> NaiveMethodOutput`
  - `mdl_method(x::Vector{Float32}, t::Float32, m::MDLMethod) -> MDLMethodOutput`
  - `deep_channel_method(x::Vector{Float32}, t::Float32, m::DeepChannelMethod) -> DeepChannelMethodOutput`
  - `method_function(::IdealizationMethod)` → returns the concrete callable for `calculate_method`
  - `calculate_method(x::Vector{Float32}, m::IdealizationMethod, t::Float32) -> MethodOutput`

- Evaluation and reconstruction:
  - `calculate_mean_square_error(data::Dict, dwell_times_approx::Vector{Float32}, dt_bins::UInt16=100) -> (mse::Float32, hist_data::Histogram, hist_approx::Histogram)`
  - `accuracy_of_idealization(actual::Vector{UInt8}, approx::Vector{UInt8}) -> Float32`
  - `mean_error(method::IdealizationMethod, t::Float32, data_size::UInt32, verbose::Bool=false) -> MeanError`
  - `actual_idealize_data(data::Dict, what_first_dict::Dict{String,Int64}, file_name::AbstractString, t::Float32) -> Vector{UInt8}`
  - `idealize_data(data::Vector{Float32}, dwell_times_approx::Vector{Float32}, hist::HistPeakAnalysis, t::Float32) -> Vector{Float32}`
  - `create_idealizations(data_folder::String, t::Float32=1e-4f0) -> Dict{String,Vector{Int8}}`

- MDL internals (exposed helpers):
  - `_mdl(segment::Vector{Float32}, BP::Vector{UInt32}) -> Float32`
  - `_test_breakpoint(segment::Vector{Float32}, candidate::Vector{UInt32}) -> Bool`
  - `detect_single_breakpoint(data::Vector{Float32}, min_seg::UInt16=300) -> Vector{UInt32}`
  - `detect_double_breakpoint(data::Vector{Float32}, min_seg::UInt16=300) -> Vector{UInt32}`
  - `detect_breaks_mdl(segment::Vector{Float32}, method::AbstractString, min_seg::UInt16=300) -> Vector{UInt32}`
  - `stepstat_mdl(data::Vector{Float32}, BP::Vector{UInt32}, threshold::Float32=0.8f0) -> (filtered::Vector{UInt32}, stepvalues::Vector{Float32})`

# Example
```
using IonChannel

x, y = read_data("trace.txt", "dwell.txt")

t = 1e-4f0
xnorm = normalize_data(Dict("x" => x))
params = MikaMethod(0.1f0, UInt16(120))
result = calculate_method(xnorm, params, t)

show_approx_on_plot(Dict("x"=>x, "dwell times"=>y), result, 0.5f0, 0.9f0, t)
```
# Authors

Created and maintained by [Piotr Mika](https://github.com/p-j-o-t-e-r).
"""
module IonChannel

  import Pkg
  Pkg.activate("../")
  packages = ["StatsBase", "Distributions", "Plots", "Normalization"]
  Pkg.add(packages)

	ENV["PYTHON_JL_RUNTIME_PYTHON"] = Sys.which("python3")
	ENV["PYTHON"] = Sys.which("python3")
  @info "Using Python at: $(ENV["PYTHON"])"
	Pkg.add(["PyCall", "PlutoUI", "OpenSSL_jll", "HypothesisTests"])
	Pkg.build("PyCall")
	using PyCall
	using OpenSSL_jll, PlutoUI, Plots

  include("./types.jl")
  include("./read_data.jl")
  include("./auxiliary.jl")
  include("./mika_method.jl")
  include("./mean_deviation_method.jl")
  include("./deep_channel_method.jl")
  include("./naive_method.jl")
  include("./method_mse.jl")
  include("./mdl.jl")
  include("./visualization.jl")
  
  export
    read_data,
    read_all_file_paths,
    create_idealizations,
    get_specified_datapoints,
    normalize_data,
    histogram_calculator,
    calculate_mean_square_error,
    show_approx_on_plot,
    mean_error,
    deviation_from_mean_method,
    IdealizationMethod,
    MeanDeviationMethod,
    method_function,
    δ,
    λ,
    calculate_method,
    MeanDeviationMethodOutput,
    MethodOutput,
    HistPeakAnalysis,
    ThresholdWidth,
    Point,
    Line,
    Noise,
    noise_data,
    μ,
    σ,
    noise,
    breakpoints,
    dwell_times_approx,
    idealized_data,
    threshold,
    threshold_index,
    noise_mse,
    MikaMethod,
    MikaMethodOutput,
    calculate_probability_histogram,
    analyze_histogram_peaks,
    point,
    line,
    combine_time_with_data,
    get_threshold_width,
    calculate_approximation,
    idealize_data,
    fit_normal_to_noise,
    fit_mse,
    mika_method,
    show_threshold_on_plot,
    plot_idealization_representation,
    noise_test,
    DeepChannelMethod,
    DeepChannelMethodOutput,
    deep_channel_method,
    NaiveMethod,
    NaiveMethodOutput,
    naive_method,
    MeanError,
    accuracy_of_idealization,
    actual_idealize_data,
    detect_single_breakpoint,
    detect_double_breakpoint,
    detect_breaks_mdl,
    _test_breakpoint,
    _mdl,
    stepstat_mdl,
    MDLMethod,
    MDLMethodOutput,
    mdl_method

end