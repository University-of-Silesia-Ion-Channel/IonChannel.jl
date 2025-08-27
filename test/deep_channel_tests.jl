using Test

# Check if we can load the IonChannel module with PyCall
const MODULE_LOADABLE = try
    # This will test if we can import the module without PyCall errors
    eval(:(import Pkg; using PyCall))  # Try to load PyCall first
    true
catch e
    @warn "Cannot load PyCall dependency: $e"
    false
end

# Basic structural tests that don't require PyCall
@testset "DeepChannel Module Structure (Basic)" begin
    # Test that the deep_channel_method.jl file exists and has expected content
    deep_channel_file = joinpath(@__DIR__, "..", "src", "deep_channel_method.jl")
    @test isfile(deep_channel_file)
    
    content = read(deep_channel_file, String)
    @test occursin("deep_channel_method", content)
    @test occursin("DeepChannelMethod", content)
    @test occursin("method_function", content)
    
    # Test that the module exports the right symbols
    main_file = joinpath(@__DIR__, "..", "src", "IonChannel.jl")
    @test isfile(main_file)
    
    main_content = read(main_file, String)
    @test occursin("DeepChannelMethod", main_content)
    @test occursin("DeepChannelMethodOutput", main_content) 
    @test occursin("deep_channel_method", main_content)
end

# Only run the main tests if PyCall is available
if MODULE_LOADABLE
    # Mock PyObject for testing
    struct MockPyObject
        predictions::Matrix{Float32}
    end

    # Mock predict method that mimics the expected behavior
    function Base.getproperty(mock::MockPyObject, name::Symbol)
        if name == :predict
            return function(input_data; batch_size=nothing, verbose=nothing)
                # Return the stored predictions
                return getfield(mock, :predictions)
            end
        else
            return getfield(mock, name)
        end
    end

    @testset "DeepChannel Types and Structure" begin
        # Test DeepChannelMethodOutput type structure
        dwell_times = Float32[0.1, 0.2, 0.15]
        breakpoints = Float32[0.1, 0.3]
        idealized_data = UInt8[0, 0, 1, 1, 0]
        
        output = IonChannel.DeepChannelMethodOutput(dwell_times, breakpoints, idealized_data)
        
        @test isa(output, IonChannel.DeepChannelMethodOutput)
        @test isa(output, IonChannel.MethodOutput)
        @test output.dwell_times_approx == dwell_times
        @test output.breakpoints == breakpoints
        @test output.idealized_data == idealized_data
        
        # Test that DeepChannelMethod can be created with mock object
        mock_predictions = rand(Float32, 10, 2)
        mock_model = MockPyObject(mock_predictions)
        method = IonChannel.DeepChannelMethod(mock_model)
        @test isa(method, IonChannel.DeepChannelMethod)
        @test isa(method, IonChannel.IdealizationMethod)
        @test method.model === mock_model
    end

    @testset "DeepChannel Method Function Mapping" begin
        # Test that method_function returns the correct function
        mock_predictions = rand(Float32, 5, 2)
        mock_model = MockPyObject(mock_predictions)
        method = IonChannel.DeepChannelMethod(mock_model)
        
        func = IonChannel.method_function(method)
        @test func === IonChannel.deep_channel_method
    end

    @testset "DeepChannel Algorithm Logic with Mock Model" begin
        # Create synthetic test data
        data = Float32[0.2, 0.1, 0.15, 0.8, 0.9, 0.85, 0.7, 0.8]  # Two-state pattern
        Δt = 0.1f0
        
        # Create mock predictions for a clear state transition pattern
        N = length(data)
        predictions = zeros(Float32, N, 2)
        
        # Create a clear state transition pattern
        for i in 1:N
            if i <= N÷2
                predictions[i, 1] = 0.9f0  # state 0 confident
                predictions[i, 2] = 0.1f0
            else
                predictions[i, 1] = 0.1f0  
                predictions[i, 2] = 0.9f0  # state 1 confident
            end
        end
        
        mock_model = MockPyObject(predictions)
        method = IonChannel.DeepChannelMethod(mock_model)
        
        # Run the deep channel method
        result = IonChannel.deep_channel_method(data, Δt, method)
        
        # Test output structure
        @test isa(result, IonChannel.DeepChannelMethodOutput)
        @test length(result.idealized_data) == length(data)
        @test all(x -> x in [0, 1], result.idealized_data)  # Binary states
        
        # Test that we get at least one state transition (given our mock model)
        @test length(result.dwell_times_approx) >= 1
        @test length(result.breakpoints) == length(result.dwell_times_approx) - 1
        
        # Test timing consistency
        @test all(dt -> dt > 0, result.dwell_times_approx)  # All dwell times positive
        @test all(bp -> bp >= 0, result.breakpoints)  # All breakpoints non-negative
        
        # Test that idealized data has the expected pattern for our mock
        # First half should be state 0, second half should be state 1
        first_half = result.idealized_data[1:length(data)÷2]
        second_half = result.idealized_data[(length(data)÷2+1):end]
        
        @test all(x -> x == 0, first_half)  # First half should be state 0
        @test all(x -> x == 1, second_half)  # Second half should be state 1
    end

    @testset "DeepChannel Edge Cases" begin
        # Test with minimal data
        data = Float32[0.5, 0.6]
        Δt = 0.1f0
        
        # Mock predictions: all state 0
        N = length(data)
        predictions = zeros(Float32, N, 2)
        predictions[:, 1] .= 0.8f0  # All state 0
        predictions[:, 2] .= 0.2f0
        
        mock_model = MockPyObject(predictions)
        method = IonChannel.DeepChannelMethod(mock_model)
        result = IonChannel.deep_channel_method(data, Δt, method)
        
        @test isa(result, IonChannel.DeepChannelMethodOutput)
        @test length(result.idealized_data) == length(data)
        
        # With no state transitions, we should still get valid output
        @test length(result.dwell_times_approx) >= 0
        @test length(result.breakpoints) >= 0
    end

    @testset "DeepChannel Breakpoint and Dwell Time Calculation" begin
        # Test the core logic of breakpoint and dwell time calculation
        data = Float32[0.1, 0.2, 0.3, 0.4]  # 4 samples
        Δt = 0.1f0
        
        # Mock predictions that create a specific state pattern: [0, 0, 1, 1]
        N = length(data)
        predictions = zeros(Float32, N, 2)
        
        for i in 1:N
            if i <= 2
                predictions[i, 1] = 0.9f0  # state 0
                predictions[i, 2] = 0.1f0
            else
                predictions[i, 1] = 0.1f0
                predictions[i, 2] = 0.9f0  # state 1
            end
        end
        
        mock_model = MockPyObject(predictions)
        method = IonChannel.DeepChannelMethod(mock_model)
        result = IonChannel.deep_channel_method(data, Δt, method)
        
        # Expected: 2 samples in state 0, then 2 samples in state 1
        # Should have 2 dwell times and 1 breakpoint
        @test length(result.dwell_times_approx) == 2
        @test length(result.breakpoints) == 1
        
        # First dwell time should be 2 * Δt = 0.2
        @test result.dwell_times_approx[1] ≈ 2 * Δt
        # Second dwell time should be 2 * Δt = 0.2  
        @test result.dwell_times_approx[2] ≈ 2 * Δt
        
        # Breakpoint should be at index 2 * Δt = 0.2
        @test result.breakpoints[1] ≈ 2 * Δt
    end
else
    @testset "DeepChannel Tests Skipped" begin
        @test_skip true  # Skip these tests when PyCall is not available
        @info "DeepChannel functional tests skipped - PyCall dependency not available"
    end
end