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

method_function(::NaiveMethod) = naive_method