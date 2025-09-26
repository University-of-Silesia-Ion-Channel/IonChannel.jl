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
        threshold::Float32=0.8f0,
        Δt::Float32
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
- Δt::Float32: Sampling interval used to convert indices to time.

Returns:
- (filtered::Vector{UInt32}, stepvalue::Vector{Float32}):
  - `filtered`: Breakpoints surviving the jump threshold.
  - `stepvalue`: Estimated mean level for each (original) segment.

Notes:
- Ensures each segment has at least one index; if an interval collapses, it uses the breakpoint index.
"""
function stepstat_mdl(data::Vector{Float32}, BP::Vector{UInt32}, threshold::Float32) :: Tuple{Vector{UInt32}, Vector{Float32}}
    stepvalue = Float32[]
    b_idxs = vcat(1, BP, length(data))
    prev_b_idx = b_idxs[1]
	for b_idx in b_idxs[2:end]
		push!(stepvalue, mean(data[prev_b_idx:b_idx]))
		prev_b_idx = b_idx
	end
    jumps = diff(stepvalue)
    filtered = BP[abs.(jumps) .> threshold]
    filtered, stepvalue
end

"""
    mdl_method_part(
        data::Vector{Float32},
        c_method::MDLMethod
    ) :: Vector{UInt32}

Recursively detect candidate breakpoints in a time series using the MDL criterion.

This function implements an iterative, segment-wise search for change points,
operating over the range `[1, length(data)]`, and checking within each segment
for additional split points using both single- and double-breakpoint detection.
Segmentations are only accepted if they significantly reduce the MDL score, and all
proposed breakpoints are consolidated, sorted, and returned as sample indices.

### Arguments
- `data::Vector{Float32}`: Input signal to segment and search for change points.
- `c_method::MDLMethod`: Configuration object containing at least the `min_seg` field
  (minimum segment length).

### Returns
- `Vector{UInt32}`: Detected breakpoint sample indices (1-based), sorted, not including the last index.

### Implementation notes
- Starts by initializing the entire range as a candidate segment.
- For each current segment, attempts a single breakpoint search via [`detect_breaks_mdl`](@ref) with `"full"`.
  If no breakpoint is found, and the segment is sufficiently long, tries a double-break search (`"full_two_break"`).
- Successfully validated breakpoints (by MDL) are accepted and the search continues recursively within new subsegments.
- Consolidates all candidate breakpoints found, sorts, and removes the terminal endpoint before returning.

### Usage example
```
breaks = mdl_method_part(data, MDLMethod(300, 0.8f0, 100))
```
"""
function mdl_method_part(data::Vector{Float32}, c_method::MDLMethod) :: Vector{UInt32}
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
    breaks
end

"""
    mdl_method(
        data::Vector{Float32},
        Δt::Float32,
        c_method::MDLMethod
    ) :: MDLMethodOutput

Segment and idealize a univariate time series using MDL-driven breakpoint detection.

This function performs bidirectional search for optimal breakpoints using
[`mdl_method_part`](@ref) in forward and reverse directions, merges all unique
breaks, and applies amplitude-jump filtering via [`stepstat_mdl`](@ref).
It then reconstructs a two-state (0/1) idealized sequence, parameters such as
breakpoint times and mean values, and estimates state dwell-times.

### Pipeline
1. **Breakpoint search**: Call [`mdl_method_part`](@ref) on `data` and its reverse.
2. **Breakpoint consolidation**: Merge, deduplicate, and sort all break indices.
3. **Filtering**: Remove breakpoints with sub-threshold amplitude jumps using [`stepstat_mdl`](@ref).
4. **Thresholding**: Estimate a classification amplitude threshold from the histogram of `data`.
5. **Idealization**: Generate a state sequence alternating between 0/1 at each surviving breakpoint, based on thresholding of the initial value.
6. **Dwell times**: Compute per-state dwell times as intervals between breakpoints (in time units).

### Arguments
- `data::Vector{Float32}`: Input sequence (trace).
- `Δt::Float32`: Sampling interval, used to express break times and dwell times in physical units.
- `c_method::MDLMethod`: Configuration containing at least fields `min_seg` (minimum segment length), `threshold` (minimum step), and number of histogram bins.

### Returns
- `MDLMethodOutput`: Struct containing:
  - `breakpoints::Vector{Float32}`: Surviving breakpoint times.
  - `dwell_times_approx::Vector{Float32}`: Duration of each state (seconds or chosen `Δt` units).
  - `idealized_data::Vector{UInt8}`: Idealized 0/1 state sequence.
  - `all_breaks::Vector{Float32}`: All candidate breakpoints (pre-filtering).
  - `step_values::Vector{Float32}`: Mean values for each segment.

### Notes
- Assumes a two-state system alternating at each detected break.
- Returns a single dwell if no breakpoints survive filtering.
- Depends on helper functions: [`mdl_method_part`](@ref), [`stepstat_mdl`](@ref), [`histogram_calculator`](@ref), [`calculate_probability_histogram`](@ref), and [`analyze_histogram_peaks`](@ref).

### Example
```
c_method = MDLMethod(300, 0.8f0, 100)
Δt = 1.0f0 # 1 sample per time unit
result = mdl_method(data, Δt, c_method)
println(result.breakpoints)
println(result.dwell_times_approx)
```
"""
function mdl_method(data::Vector{Float32}, Δt::Float32, c_method::MDLMethod) :: MDLMethodOutput
    breaks_forward = mdl_method_part(data, c_method)
    breaks_backward = (length(data) + 1) .- mdl_method_part(data[end:-1:1], c_method)
    all_breaks::Vector{UInt32} = sort(unique(vcat(breaks_forward, breaks_backward)))
    # @info "$(length(all_breaks)) breakpoints detected before step filtering"
    final_breaks, step_values = stepstat_mdl(data, all_breaks, c_method.threshold)
    # @info "$final_breaks breakpoints after step filtering"
	breakpoints::Vector{Float32} = final_breaks .* Δt
	histogram_of_data = histogram_calculator(data)
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
    if !(isempty(breakpoints))
	    dwell_times = vcat([breakpoints[1]], diff(breakpoints))
    else
        dwell_times = [length(data) * Δt]
    end
	MDLMethodOutput(breakpoints, dwell_times, idealized_data, all_breaks .* Δt, step_values)
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