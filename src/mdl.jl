using StatsBase

"""
    _mdl(segment::Vector{Float32}, BP::Vector{UInt32}) :: Float32

Compute the Minimum Description Length (MDL) criterion for a piecewise-constant segmentation.

Given a signal `segment` and a set of candidate breakpoint indices `BP`, this function
evaluates the MDL cost consisting of a model complexity term and a data-fit term
(residual sum of squares within segments). Lower values indicate better segmentations.

The segmentation is induced by the sorted, unique set `{1} ∪ BP ∪ {N}`, where `N = length(segment)`.

Arguments:
- segment::Vector{Float32}: The data segment to be evaluated.
- BP::Vector{UInt32}: Candidate breakpoint indices (1-based, strictly within 1..N).

Returns:
- Float32: The MDL value; lower is better. Returns `Inf` if the residual sum of squares (RSS) is non-positive.

Details:
- For each segment between consecutive breakpoints, the mean is estimated and RSS accumulated.
- Complexity term includes `p*log(N)` with `p = number_of_segments - 1`, and a local length penalty `0.5*Σ log(Nseg)`.
- Fit term is `(N/2)*log(RSS/N)`.

Notes:
- Empty subsegments are skipped when accumulating RSS and complexity.
"""
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

"""
    _test_breakpoint(segment::Vector{Float32}, candidate::Vector{UInt32}) :: Bool

Test whether adding the proposed breakpoints reduces the MDL criterion.

Computes MDL with and without the candidate breakpoints and returns true if the
candidate segmentation yields a strictly lower MDL score.

Arguments:
- segment::Vector{Float32}: Data segment to test.
- candidate::Vector{UInt32}: Proposed breakpoint indices.

Returns:
- Bool: `true` if `MDL(with candidate) < MDL(without)`, otherwise `false`.

Notes:
- Returns `false` when `candidate` is empty.
"""
function _test_breakpoint(segment::Vector{Float32}, candidate::Vector{UInt32}) :: Bool
    if length(candidate) == 0
        return false
    end
    mdl_no = _mdl(segment, Vector{UInt32}([]))
    mdl_yes = _mdl(segment, candidate)
    return mdl_no > mdl_yes
end

"""
    detect_breaks_mdl(
        segment::Vector{Float32},
        method::AbstractString,
        min_seg::UInt16=UInt16(300)
    ) :: Vector{UInt32}

Detect candidate breakpoint(s) using MDL-backed single or double-break search.

Depending on `method`, the function delegates to a single- or double-break detector
and then validates the proposed breakpoints with the MDL test.

Arguments:
- segment::Vector{Float32}: Data to search for breakpoints.
- method::AbstractString: Either `"full"` (single-break) or `"full_two_break"` (double-break).
- min_seg::UInt16: Minimum allowed segment length (in samples).

Returns:
- Vector{UInt32}: Validated breakpoint indices. Empty if none pass the MDL test or method is unknown.
"""
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

"""
    detect_single_breakpoint(
        data::Vector{Float32},
        min_seg::UInt16=UInt16(300)
    ) :: Vector{UInt32}

Find a single change point that minimizes within-segment squared error, subject to a minimum segment length.

Performs a linear scan, maintaining incremental means and within-segment sums of squares to
identify the index that best splits the data into two segments with minimal total squared error.

Arguments:
- data::Vector{Float32}: Input sequence.
- min_seg::UInt16: Minimum length for each side of the breakpoint.

Returns:
- Vector{UInt32}: A vector with one breakpoint index, or empty if no valid split is found.

Notes:
- If `length(data) < 2*min_seg`, no split is attempted and the result is empty.
"""
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

"""
    detect_double_breakpoint(
        data::Vector{Float32},
        min_seg::UInt16=UInt16(300)
    ) :: Vector{UInt32}

Search for two change points that jointly minimize the sum of within-segment squared errors.

Uses cumulative sums to evaluate candidate pairs `(i, j)` efficiently, enforcing a minimum segment
length on all three resulting segments.

Arguments:
- data::Vector{Float32}: Input signal.
- min_seg::UInt16: Minimum segment length for each of the three segments.

Returns:
- Vector{UInt32}: A 2-element vector `[i, j]` with the best breakpoints, or empty if none.

Notes:
- If `length(data) < 3*min_seg`, returns an empty vector without searching.
"""
function detect_double_breakpoint(data::Vector{Float32}, min_seg::UInt16=UInt16(300)) ::Vector{UInt32}
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

"""
    stepstat_mdl(
        data::Vector{Float32},
        BP::Vector{UInt32},
        threshold::Float32=0.8f0
    ) :: Tuple{Vector{UInt32}, Vector{Float32}}

Estimate step values per segment and filter breakpoints by jump magnitude.

Given breakpoints `BP`, this function:
- appends the end index to form closed segments,
- estimates the mean (`stepvalue`) for each segment,
- computes jumps between consecutive segment means,
- filters breakpoints whose absolute jump exceeds `threshold`.

Arguments:
- data::Vector{Float32}: Input signal.
- BP::Vector{UInt32}: Candidate breakpoints (1-based).
- threshold::Float32: Minimum absolute difference between consecutive step means to retain a breakpoint.

Returns:
- (filtered::Vector{UInt32}, stepvalue::Vector{Float32}):
  - `filtered`: Breakpoints surviving the jump threshold.
  - `stepvalue`: Estimated mean level for each (original) segment.

Notes:
- Ensures each segment has at least one index; if an interval collapses, it uses the breakpoint index.
"""
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

"""
    mdl_method(
        data::Vector{Float32},
        Δt::Float32,
        c_method::MDLMethod
    ) :: MDLMethodOutput

Perform MDL-based idealization by iteratively detecting and validating breakpoints,
filtering spurious changes by step size, and constructing the idealized state sequence.

Pipeline:
1. Iterative breakpoint search:
   - Starting from the full range, repeatedly detect single or double breakpoints
     (via [`detect_breaks_mdl`](@ref)) within the current segment, respecting `c_method.min_seg`.
   - Accept and insert any proposed breakpoints that pass the MDL test.
2. Sort and consolidate all local breakpoints.
3. Filter by jump magnitude:
   - Use `stepstat_mdl(data, breaks, c_method.threshold)` to remove small steps and
     estimate segment means.
4. Convert sample indices to time:
   - `breakpoints = final_breaks .* Δt`.
5. Determine initial state using a histogram-derived threshold:
   - Estimate amplitude threshold via [`histogram_calculator`](@ref), [`calculate_probability_histogram`](@ref),
     and [`analyze_histogram_peaks`](@ref), and set initial state to 0 if `data[1] < threshold`, else 1.
6. Build per-sample idealized sequence by alternating states across `final_breaks`.
7. Compute dwell times as `[breakpoints[1]; diff(breakpoints)]`.

Arguments:
- data::Vector{Float32}: Input trace to be idealized.
- Δt::Float32: Sampling interval used to convert indices to time.
- c_method::MDLMethod: Configuration with `min_seg`, `threshold`, and `number_of_histogram_bins`.

Returns:
- MDLMethodOutput: Contains `breakpoints` (Float32 times), `dwell_times_approx`, and `idealized_data` (UInt8 states).

Notes:
- Expects the following helpers in scope: [`histogram_calculator`](@ref), [`calculate_probability_histogram`](@ref),
  and [`analyze_histogram_peaks`](@ref) returning an object with fields `edges` and `pmin_index`.
- State alternation assumes a two-state model (0/1) switching at each retained breakpoint.
- If `final_breaks` is empty, ensure calling code handles empty dwell times accordingly.
"""
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

"""
    method_function(::MDLMethod)

Dispatch helper that maps an [`MDLMethod`](@ref) configuration to its execution function.

Returns a callable with signature `(data::Vector{Float32}, Δt::Float32, c_method::MDLMethod) -> MDLMethodOutput`,
typically used in higher-level code to select the appropriate idealization routine based on method type.

Example:
```
f = method_function(MDLMethod(300, 0.8f0, 100))
out = f(data, Δt, c_method) # calls mdl_method
```
"""
method_function(::MDLMethod) = mdl_method