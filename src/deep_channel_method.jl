using PyCall

method_function(::DeepChannelMethod) = deep_channel_method
const np = pyimport("numpy")

function deep_channel_method(data::Vector{Float64}, Δt::Float64, c_method::DeepChannelMethod)
    data_augmentation = 0
    N = length(data) + data_augmentation
    t = fit(UnitRangeTransform, data)
    scaled_data = StatsBase.transform(t, data)
    input_data = reshape(scaled_data, (N, 1, 1, 1))
    pred = c_method.model.predict(input_data, batch_size=16*1024, verbose=0)
    class_predict_val = np.argmax(pred, axis=-1)
    # class_predict_val = [ind[2] - 1 for ind in class_predict_val_ind]
    
    # create dwell times for the output of the model
    prev_value = class_predict_val[1]
    t = 0
    dwell_times_approx = Vector{Float64}([])
    breakpoints = Vector{Float64}([])
    
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