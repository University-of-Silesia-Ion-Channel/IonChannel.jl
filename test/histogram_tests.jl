using Test
using StatsBase

@testset "histogram_calculator and related" begin
    # simple bimodal data: half around -1, half around +2
    data = vcat(fill(-1.0f0, 50), fill(2.0f0, 50))
    hist = IonChannel.histogram_calculator(data, UInt16(10))

    @test isa(hist, Histogram)
    # edges should be bins+1
    @test length(hist.edges[1]) == 11
    # total weight should be positive and no larger than number of samples
    @test 0 < sum(hist.weights) <= length(data)

    # probability histogram should sum to ~1.0
    prob = IonChannel.calculate_probability_histogram(hist)
    @test abs(sum(prob.weights) - 1.0f0) < 1e-6

    # analyze peaks: expect two peaks roughly around the two clusters
    analysis = IonChannel.analyze_histogram_peaks(prob)
    @test isa(analysis, HistPeakAnalysis)
    @test analysis.pmax1_index != analysis.pmax2_index

    # threshold width for small epsilon should be between min and max edges
    thr = IonChannel.get_threshold_width(analysis, 0.1f0)
    @test thr.x₁ >= minimum(analysis.edges) && thr.x₂ <= maximum(analysis.edges)
end
