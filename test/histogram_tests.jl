using Test
using StatsBase

@testset "histogram_calculator and related" begin
    # simple bimodal data: half around -1, half around +2
    data = Dict("x" => vcat(fill(-1.0f0, 50), fill(2.0f0, 50)))
    scaled_data = IonChannel.normalize_data(data)
    hist = IonChannel.histogram_calculator(scaled_data)

    @test isa(hist, Histogram)
    # edges should be bins+1
    # @test length(hist.edges[1]) == 11
    # total weight should be positive and no larger than number of samples
    @test 0 < sum(hist.weights) <= length(scaled_data)

    # probability histogram should sum to ~1.0
    prob = IonChannel.calculate_probability_histogram(hist)
    @test abs(sum(prob.weights) * prob.edges[1].step.hi - 1.0f0) < 1e-6

    # analyze peaks: expect two peaks roughly around the two clusters
    analysis = IonChannel.analyze_histogram_peaks(scaled_data)
    @test isa(analysis, IonChannel.HistPeakAnalysis)
    @test analysis.left_peak_index != analysis.right_peak_index

    # threshold width for small epsilon should be between min and max edges
    thr = IonChannel.get_threshold_width(analysis, 0.1f0)
    @test thr.x₁ >= minimum(analysis.edges) && thr.x₂ <= maximum(analysis.edges)
end
