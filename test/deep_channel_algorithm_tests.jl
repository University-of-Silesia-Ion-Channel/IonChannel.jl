using Test

"""
Additional unit tests for DeepChannel that focus on the algorithm logic
without requiring full module dependencies. These tests verify the core
computational patterns used in the deep_channel_method function.
"""

@testset "DeepChannel Algorithm Logic Patterns" begin
    # Test argmax logic similar to what's used in deep_channel_method
    @testset "Argmax Prediction Processing" begin
        # Simulate predictions from a binary classifier
        predictions = Float32[
            0.8 0.2;  # Sample 1: class 0 more likely 
            0.3 0.7;  # Sample 2: class 1 more likely
            0.9 0.1;  # Sample 3: class 0 more likely
            0.1 0.9   # Sample 4: class 1 more likely
        ]
        
        # Test argmax extraction (like in deep_channel_method)
        idxs = argmax(predictions; dims=2)
        class_predict_val = vec(getindex.(idxs, 2)) .- 1  # Convert to 0-based indexing
        
        @test class_predict_val == UInt64[0, 1, 0, 1]  # Expected state sequence
    end
    
    @testset "Dwell Time and Breakpoint Calculation" begin
        # Simulate the core logic from deep_channel_method
        class_predict_val = UInt8[0, 0, 0, 1, 1, 0, 0, 1]  # State sequence
        Δt = 0.1f0
        data_augmentation = 0  # As used in the actual implementation
        
        # Replicate the dwell time calculation logic exactly
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
        
        # Apply post-processing (only if there are dwell times)
        if !isempty(dwell_times_approx)
            dwell_times_approx[1] -= data_augmentation * Δt
        end
        breakpoints = breakpoints .- data_augmentation * Δt
        
        # Expected: [0,0,0] -> [1,1] -> [0,0] -> [1]
        # Transitions at: index 3 (3*Δt), index 5 (2*Δt), index 7 (2*Δt)
        # Note: the algorithm only records dwell times when transitions occur
        expected_dwells = Float32[0.3, 0.2, 0.2]  # Only records completed segments
        expected_breaks = Float32[0.3, 0.5, 0.7]  # At indices 3, 5, 7
        
        @test length(dwell_times_approx) == 3  # Only transitions create dwell times
        @test dwell_times_approx ≈ expected_dwells
        @test breakpoints ≈ expected_breaks
    end
    
    @testset "Data Scaling Patterns" begin
        # Test data preprocessing similar to UnitRangeTransform
        data = Float32[1.0, 3.0, 2.0, 5.0, 4.0]
        
        # Manual unit range transform (scales to [0,1])
        min_val = minimum(data)
        max_val = maximum(data)
        scaled_data = (data .- min_val) ./ (max_val - min_val)
        
        @test minimum(scaled_data) ≈ 0.0f0
        @test maximum(scaled_data) ≈ 1.0f0
        @test issorted(scaled_data[[1,3,2,5,4]])  # Check relative ordering preserved
        
        # Test reshaping to expected input format (N, 1, 1, 1)
        N = length(scaled_data)
        input_data = reshape(scaled_data, (N, 1, 1, 1))
        @test size(input_data) == (N, 1, 1, 1)
        @test input_data[:, 1, 1, 1] == scaled_data
    end
    
    @testset "Edge Cases in State Processing" begin
        # Test with no state transitions (all same state)
        class_predict_val = UInt8[0, 0, 0, 0]
        Δt = 0.1f0
        
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
        
        # No transitions should result in empty arrays
        @test isempty(dwell_times_approx)
        @test isempty(breakpoints)
    end
    
    @testset "Single Transition Pattern" begin
        # Test with exactly one transition
        class_predict_val = UInt8[0, 0, 1, 1]
        Δt = 0.05f0
        data_augmentation = 0
        
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
        
        # Apply post-processing
        if !isempty(dwell_times_approx)
            dwell_times_approx[1] -= data_augmentation * Δt
        end
        breakpoints = breakpoints .- data_augmentation * Δt
        
        # Should have one dwell time (for the first segment) and one breakpoint
        @test length(dwell_times_approx) == 1  # Only the completed first segment
        @test length(breakpoints) == 1
        @test dwell_times_approx[1] ≈ 2 * Δt  # 2 samples in state 0
        @test breakpoints[1] ≈ 2 * Δt  # Transition at index 2
    end
    
    @testset "Alternating States Pattern" begin
        # Test rapidly alternating states
        class_predict_val = UInt8[0, 1, 0, 1, 0]
        Δt = 0.1f0
        data_augmentation = 0
        
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
        
        # Apply post-processing
        if !isempty(dwell_times_approx)
            dwell_times_approx[1] -= data_augmentation * Δt
        end
        breakpoints = breakpoints .- data_augmentation * Δt
        
        # Every transition creates a single-sample dwell time
        @test length(dwell_times_approx) == 4  # 4 transitions = 4 completed segments
        @test length(breakpoints) == 4  # 4 transitions
        @test all(dt -> dt ≈ Δt, dwell_times_approx)  # All single-sample dwells
    end
end

@testset "DeepChannel Data Validation Patterns" begin
    @testset "Input Data Validation" begin
        # Test various input data scenarios
        
        # Normal case
        data = Float32[0.1, 0.5, 0.9, 0.3]
        @test all(x -> isa(x, Float32), data)
        @test length(data) > 0
        
        # Edge case: single sample
        single_data = Float32[0.5]
        @test length(single_data) == 1
        @test isa(single_data[1], Float32)
        
        # Edge case: very small values
        small_data = Float32[1e-6, 2e-6, 3e-6]
        min_val = minimum(small_data)
        max_val = maximum(small_data)
        @test max_val > min_val  # Ensure range is non-zero
    end
    
    @testset "Timing Parameter Validation" begin
        # Test different Δt values
        valid_deltas = Float32[0.001, 0.01, 0.1, 1.0]
        
        for Δt in valid_deltas
            @test Δt > 0.0f0
            @test isa(Δt, Float32)
            
            # Test timing calculations
            dwell_samples = 5
            expected_time = dwell_samples * Δt
            @test expected_time > 0.0f0
        end
    end
    
    @testset "Output Structure Validation" begin
        # Test expected output structure patterns
        
        # Simulate a typical output
        dwell_times = Float32[0.1, 0.2, 0.15]
        breakpoints = Float32[0.1, 0.3]
        idealized_data = UInt8[0, 0, 1, 1, 0]
        
        # Validate structure
        @test length(breakpoints) == length(dwell_times) - 1
        @test all(dt -> dt > 0, dwell_times)
        @test all(bp -> bp >= 0, breakpoints)
        @test issorted(breakpoints)  # Breakpoints should be in order
        @test all(state -> state in [0, 1], idealized_data)  # Binary states
    end
end