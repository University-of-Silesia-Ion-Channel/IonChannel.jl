"""
    module IonChannel

A Julia module for advanced idealization and analysis of ion channel single-molecule data.

# Overview

`IonChannel` provides type-safe, highly modular tools for the idealization, segmentation, and analysis of ion channel current traces. It implements and integrates methods based on running-mean deviation, histogram-based segmentation, threshold optimization, and noise fitting, supporting workflows in quantitative biophysics and single-molecule electrophysiology.

## Main Features

- **Idealization Methods:**  
  - Abstract parameter types and output structures for extensible idealization algorithms.
  - Implementation of classic deviation-from-mean and Mika algorithm with optimized thresholding.

- **Data Handling:**  
  - Automated reading of raw and dwell-time data files with directory structure support.
  - Utility to extract, segment, and normalize time series data.

- **Histogram Analysis:**  
  - Probability histogram generation and peak/trough analysis.
  - Automated detection of major states, threshold bands, and signal features.

- **Noise and Fit Evaluation:**  
  - Construction of noise models from data vs. idealized results.
  - Statistical fitting, MSE evaluation, and retrieval of core metrics.

- **Visualization:**  
  - Interactive and batch plotting functions to visualize data, thresholds, and segmentation results.

## Typical Workflow

1. **Reading Data:**  
   Use `read_data` or `read_all_file_paths` to load trace data and associated dwell times.

2. **Idealization:**  
   Choose an idealization method struct (e.g. `MeanDeviationMethod`, `MikaMethod`) and call `calculate_method` to obtain an output struct with breakpoints and dwell times.

3. **Histogram and Threshold:**  
   Generate and analyze histograms with `histogram_calculator` and `analyze_histogram_peaks`. Compute optimal threshold bands with `get_threshold_width`.

4. **Visualization:**  
   Visualize results using plotting functions like `show_approx_on_plot`,  `show_threshold_on_plot` and `plot_idealization_representation`.

5. **Noise Analysis:**  
   Quantify and fit noise using the `Noise` struct and `fit_normal_to_noise`, and evaluate fit quality with `fit_mse`.

## Extending the Module

To add a new idealization method:
1. Create a new parameter struct `<: IdealizationMethod`.
2. Implement the algorithm returning a `<: MethodOutput` output.
3. Link your method via `method_function(::YourMethodType) = your_algorithm`.

## Key Types

- `IdealizationMethod`, `MethodOutput`
- `MeanDeviationMethod`, `MikaMethod`, and associated output types
- `HistPeakAnalysis`, `ThresholdWidth`, `Noise`

## References & Help

All public-facing types and major functions are documented with Julia docstrings—use the help system (`?function_or_type`) in Julia or Pluto to browse API details and usage examples.

# Example
```
using IonChannel

# Load data
data, dwell = read_data("trace.txt", "dwell.txt")

# Normalize and idealize
params = MikaMethod(ϵ=0.1, number_of_histogram_bins=120)
result = calculate_method(data, params, 1e-4)

# Visualize results
show_approx_on_plot(Dict("x"=>data, "dwell times"=>dwell), result, 0.5, 1e-4)
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