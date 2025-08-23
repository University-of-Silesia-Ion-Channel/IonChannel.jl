using StatsBase

function _mdl(segment::Vector{Float32}, BP::Vector{UInt32})::Float32
    N = length(segment)
    BPi = Vector{UInt32}(unique(vcat([1], BP, [N])))
    p = length(BPi) - 1
    RSS = 0.0
    CL = 0.0
    for k in eachindex(BPi)[2:end]
        seg = BPi[k-1]:(BPi[k]-1)
        if length(seg) == 0
            continue
        end
        mu = mean(segment[seg])
        RSS += sum((segment[seg] .- mu).^2)
        Nseg = length(seg)
        if Nseg > 0
            CL += log(Nseg)
        end
    end
    if RSS <= 0
        return Inf
    end
    p*log(N) + 0.5*CL + (N/2)*log(RSS/N)
end

function _test_breakpoint(segment::Vector{Float32}, candidate::Vector{UInt32}) :: Bool
    if length(candidate) == 0
        return false
    end
    mdl_no = _mdl(segment, Vector{UInt32}([]))
    mdl_yes = _mdl(segment, candidate)
    return mdl_no > mdl_yes
end

function detect_breaks_mdl(segment::Vector{Float32}, method::AbstractString, min_seg::UInt16=UInt16(300))
    if method == "full"
        candidate = detect_single_breakpoint(segment, min_seg)
    elseif method == "full_two_break"
        candidate = detect_double_breakpoint(segment, min_seg)
    else
        return Vector{UInt32}([])
    end
    if _test_breakpoint(segment, candidate)
        ret = candidate
        # @info "MDL detected $(ret) breakpoint"
    else
        ret = Vector{UInt32}([])    
    end
    ret
end

function detect_single_breakpoint(data::Vector{Float32}, min_seg::UInt16=UInt16(300))::Vector{UInt32}
    n = length(data)
    if n < 2 * min_seg
        return Vector{UInt32}([])
    end

    mean1 = mean(data[1:min_seg - 1])
    mean2 = mean(data[min_seg:end])
    logL1 = sum((data[1:min_seg - 1] .- mean1) .^ 2)
    logL2 = sum((data[min_seg:end] .- mean2) .^2)
    bestlog = logL1 + logL2
    best_idx = 0

    for i in min_seg + 1:n - min_seg - 1
        new_mean1 = ((i - 1)*mean1 + data[i]) / i
        diff1 = new_mean1 - mean1

        new_mean2 = ((n - i + 1)*mean2 - data[i])/(n - i)
        diff2 = new_mean2 - mean2

        newL1 = logL1 + (data[i] - mean1)^2 - i*(diff1^2)
        newL2 = logL2 - (data[i] - new_mean2)^2 + (n - i + 1)*(diff2^2)
        Nloglik = newL1 + newL2

        if Nloglik < bestlog
            bestlog = Nloglik
            best_idx = i
        end
        mean1 = new_mean1
        mean2 = new_mean2
        logL1 = newL1
        logL2 = newL2
    end
    if best_idx == 0
        return Vector{UInt32}([])
    end
    out = Vector{UInt32}([best_idx])
    out
end

function detect_double_breakpoint(data::Vector{Float32}, min_seg::UInt16=UInt16(300))
    n = length(data)
    if n < 3 * min_seg
        return Vector{Int32}([])
    end
    cumx = cumsum(data)
    cumz = cumsum(data .^ 2)
    cumx_end = cumx[end]
    cumz_end = cumz[end]

    bestlog = Inf
    best_i = 0
    best_j = 0

    for i in min_seg:n - 2 * min_seg - 1
        cumx_i = cumx[i]
        cumz_i = cumz[i]
        mu1 = cumx_i/(i+1)
        logL1 = cumz_i - 2*mu1*cumx_i + (i+1)*mu1^2

        for j in i + min_seg:n - min_seg - 1
            l2 = j - i
            l3 = n - j
            cumx_j = cumx[j]
            cumz_j = cumz[j]

            mu2 = (cumx_j - cumx_i)/l2
            logL2 = cumz_j - cumz_i - 2*mu2*(cumx_j - cumx_i) + l2*mu2^2

            mu3 = (cumx_end - cumx_j)/l3
            logL3 = cumz_end - cumz_j - 2*mu3*(cumx_end - cumx_j) + l3*mu3^2

            Nloglik = logL1 + logL2 + logL3
            if Nloglik < bestlog
                bestlog = Nloglik
                best_i = i
                best_j = j
            end
        end
    end
    if best_i == 0 || best_j == 0
        return Vector{UInt32}([])
    end
    out = Vector{UInt32}(undef, 2)
    out[1] = best_i
    out[2] = best_j
    out
end

function stepstat_mdl(data::Vector{Float32}, BP::Vector{UInt32}, threshold::Float32=Float32(0.8)) :: Tuple{Vector{UInt32}, Vector{Float32}}
    push!(BP, UInt32(length(data)))
    stepvalue = zeros(Float32, length(BP))
    skip::UInt32 = 1
    i0::UInt32 = BP[1]
    for k in eachindex(BP)
        start::UInt32 = i0 + skip
        stop::UInt32 = BP[k] - skip
        if stop < start
            start = BP[k]
            stop = BP[k]
        end
        indices = start:stop+1
        if length(indices) == 0
            indices = Vector{UInt32}([BP[k]])
        end
        stepvalue[k] = mean(data[indices])
        i0 = BP[k]
    end

    jumps = diff(stepvalue)
    filtered = BP[1:end - 1][abs.(jumps) .> threshold]
    filtered, stepvalue
end

function mdl_method(data::Vector{Float32}, Δt::Float32, c_method::MDLMethod) :: MDLMethodOutput
    
	start::UInt32 = 1
	end_::UInt32 = length(data)
    BP_local = Vector{UInt32}([end_])
    BPlast::UInt32 = end_
    t0::UInt32 = start
    currentBP::UInt32 = BPlast

    while t0 < end_
        while true
			# @info "Searching breakpoints at: $t0:$(currentBP-1)"
            current_segment = data[t0:(currentBP-1)]
            br = detect_breaks_mdl(current_segment, "full", c_method.min_seg)
            if isempty(br) && length(current_segment) > 3 * c_method.min_seg
                br = detect_breaks_mdl(current_segment, "full_two_break", c_method.min_seg)
            end

            if !isempty(br)
                loc = Vector{UInt32}(br .+ t0 .- 1)
                BP_local = vcat(BP_local, loc)
                currentBP = loc[1]
            else
                break
            end
        end

        sort!(BP_local)
        t0 = currentBP + 1
        if currentBP != BPlast && currentBP != end_
            BPlast = currentBP
            idx = findall(x -> x == currentBP, BP_local)
            if !isempty(idx) && idx[1] + 1 <= length(BP_local)
                currentBP = BP_local[idx[1] + 1]
            end
        end
    end
    
    breaks = sort(BP_local)[1:end-1]
	final_breaks, step_values = stepstat_mdl(data, breaks, c_method.threshold)
	breakpoints = final_breaks .* Δt

	histogram_of_data = histogram_calculator(data, c_method.number_of_histogram_bins)
    prob_hist = calculate_probability_histogram(histogram_of_data)
    hist_analysis = analyze_histogram_peaks(prob_hist)

	threshold = hist_analysis.edges[hist_analysis.pmin_index]
	if data[1] < threshold
        current_state = 0 # starting at the bottom
    else
        current_state = 1 # starting at the top
    end
	prev_br_idx = 1
	idealized_data = [current_state]
	for br_idx in final_breaks
		append!(idealized_data, fill(current_state, br_idx - prev_br_idx))
		prev_br_idx = br_idx
		current_state = current_state == 0 ? 1 : 0
	end
	append!(idealized_data, fill(current_state, length(data) - prev_br_idx))
	dwell_times = append!([breakpoints[1]], diff(breakpoints))
	MDLMethodOutput(breakpoints, dwell_times, idealized_data)
end

method_function(::MDLMethod) = mdl_method