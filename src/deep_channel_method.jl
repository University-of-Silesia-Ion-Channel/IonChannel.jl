using PyCall

"""
    method_function(::DeepChannelMethod) -> Function

Return the algorithm function associated with a [`DeepChannelMethod`](@ref).

This allows code like `calculate_method(data, m, Δt)` to work for any
[`IdealizationMethod`](@ref) subtype without changing the executor logic.
"""
method_function(::DeepChannelMethod) = deep_channel_method

"""
    deep_channel_method(data::Vector{Float32}, Δt::Float32, c_method::DeepChannelMethod) 
        :: DeepChannelMethodOutput

Apply a deep learning–based ion-channel state detection method to time series data.

# Arguments
- `data::Vector{Float32}`: The raw time-series signal (e.g., ion-channel recording).
- `Δt::Float32`: Sampling interval of the signal. Used to convert sample indices into time units.
- `c_method::DeepChannelMethod`: A trained deep learning model encapsulated in a [`DeepChannelMethod`](@ref)
   object. The underlying model must support `.predict`.

# Returns
A [`DeepChannelMethodOutput`](@ref) containing:
- `dwell_times_approx::Vector{Float32}`: Approximate dwell times (time spent in each state).
- `breakpoints::Vector{Float32}`: Time points where state transitions occur.
- `class_predict_val::Vector{UInt8}`: Per-sample predicted state labels, starting at `0`.

# Method
1. The input data is scaled into `[0,1]` using `UnitRangeTransform`.
2. Data is reshaped into `(N,1,1,1)` to match the model’s expected input.
3. Predictions are computed with `c_method.model.predict`.
4. The most likely class per sample is extracted with `argmax`, producing
   `class_predict_val`.
5. Iterates through predictions to detect state transitions:
   - Accumulates dwell times per state.
   - Records breakpoints at transitions.
6. Returns dwell times, breakpoints, and the full predicted class sequence.

# Notes
- Predicted states are output as `UInt8` values starting from `0`.
- Breakpoints and dwell times are expressed in the same units as `Δt`.
- `data_augmentation` is currently inactive but reserved for possible extensions.

# Example
```
trace = rand(Float32, 10_000) # synthetic time-series
Δt = 0.1f0 # sampling interval
method = DeepChannelMethod(model) # previously constructed model

result = deep_channel_method(trace, Δt, method)

println(result.dwell_times_approx)
println(result.breakpoints)
println(result.class_predict_val[1:20])
```
"""
function deep_channel_method(data::Vector{Float32}, Δt::Float32, c_method::DeepChannelMethod) :: DeepChannelMethodOutput
    data_augmentation = 0
    N = length(data) + data_augmentation
    t = fit(UnitRangeTransform, data)
    scaled_data = StatsBase.transform(t, data)
    input_data = reshape(scaled_data, (N, 1, 1, 1))
    pred = c_method.model.predict(input_data, batch_size=16*1024, verbose=0)
    idxs = argmax(pred; dims=2)
    class_predict_val = Vector{UInt8}(vec(getindex.(idxs, 2)) .- 1)
    
    # create dwell times for the output of the model
    prev_value = class_predict_val[1]
    t = 0
    dwell_times_approx = Vector{Float32}([])
    breakpoints = Vector{Float32}([])
    
    for (i, value) in enumerate(class_predict_val[2:end])
        t += 1
        if value != prev_value
            push!(dwell_times_approx, t*Δt)
            push!(breakpoints, i*Δt)
            prev_value = value
            t = 0
        end
    end
    dwell_times_approx[1] -= data_augmentation * Δt
    DeepChannelMethodOutput(dwell_times_approx, breakpoints .- data_augmentation * Δt, class_predict_val)
end