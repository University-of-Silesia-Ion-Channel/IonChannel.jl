

"""
    method_function(::MeanDeviationMethod) -> Function

Return the algorithm function associated with a `MeanDeviationMethod`.

This allows code like `calculate_method(data, m, Δt)` to work for any
`IdealizationMethod` subtype without changing the executor logic.
"""
method_function(::MeanDeviationMethod) = deviation_from_mean_method

"""
    deviation_from_mean_method(data::Vector{Float32}, Δt::Float32, c_method::MeanDeviationMethod) -> Vector{Float32}

Estimate dwell times by detecting deviations from a running mean.

# Arguments
- `data::Vector{Float32}` - The (usually normalized) signal to analyse.
- `Δt::Float32` - Sampling interval in seconds.
- `m::MeanDeviationMethod` - Parameters container; provides `m.δ` and `m.λ`.

# Returns
- `Vector{Float32}` - Estimated dwell times in seconds.

# Algorithm
1. Start with the first value as the mean.
2. For each new sample:
- If *not* in "different state" mode, update the running mean.
- Compute `abs(sample - mean) - m.δ`.
- If above `m.λ`, trigger a state change and record dwell segment length.
- Otherwise, continue counting the current dwell segment.
3. Return the list of dwell segment durations.

# Example
```
signal = [0.0, 0.1, 0.2, 1.5, 1.6, 1.7, 0.2, 0.1, 0.05]

m = MeanDeviationMethod(0.05, 0.5)

Δt = 0.1 # 100 ms sampling interval

dwell_times_est = deviation_from_mean_method(signal, Δt, m)

println(dwell_times_est) # e.g., [0.3, 0.3, 0.3]
```

# Notes
- Produces only complete dwell segments; last partial segment is not appended.
- Works best on normalized or detrended signals to remove baseline drift.
- `δ(c_method)` suppresses detection of very small fluctuations.
- `λ(c_method)` controls the sensitivity: smaller values detect more frequent changes.
"""
function deviation_from_mean_method(data::Vector{Float32}, Δt::Float32, c_method::MeanDeviationMethod) :: MeanDeviationMethodOutput
    temporary_dwell_time = 1
    different_state_than_mean = false

    dwell_times_approx = Vector{Float32}([])

    sum_t = data[1]
    nr_means = 1
    mean_t = data[1]

    state = data[1] < mean(data) ? 0 : 1
    idealized_data = [state]
    for t in 2:length(data)
        # Step 1: Running mean
        if !different_state_than_mean
            sum_t += data[t]
            nr_means += 1
            mean_t = sum_t / nr_means
        end        

        # Step 2: Deviation
        deviation = abs(data[t] - mean_t) - δ(c_method)
        
        # Step 3: Check for change
        if deviation > λ(c_method)
            if different_state_than_mean
                temporary_dwell_time += 1
            else
                push!(dwell_times_approx, temporary_dwell_time * Δt)
                temporary_dwell_time = 1
                state = state == 0 ? 1 : 0
                
            end
            different_state_than_mean = true
            else
            if different_state_than_mean
                push!(dwell_times_approx, temporary_dwell_time * Δt)
                temporary_dwell_time = 1
                state = state == 0 ? 1 : 0
                
            else
                temporary_dwell_time += 1
            end
            different_state_than_mean = false
        end
        push!(idealized_data, state)
    end
    breakpoints = cumsum(dwell_times_approx)
    MeanDeviationMethodOutput(breakpoints, dwell_times_approx, idealized_data)
end