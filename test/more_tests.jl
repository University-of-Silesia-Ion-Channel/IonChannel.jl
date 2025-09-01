using Test
using Random

@testset "IO and utility functions" begin
    # create temporary files for read_data
    tf_data = tempname()  # VERY COOL!|
    tf_dwell = tempname() # <---------/
    open(tf_data, "w") do io
        write(io, "1.0\n2.0\n3.0\n")
    end
    open(tf_dwell, "w") do io
        write(io, "0.1\n0.2\n")
    end

    x, y = IonChannel.read_data(tf_data, tf_dwell)
    @test x == Float32[1.0, 2.0, 3.0]
    @test y == Float32[0.1, 0.2]

    # combine_time_with_data
    pairs = IonChannel.combine_time_with_data(Float32[0.0, 1.0, 2.0], 0.1f0)
    @test length(pairs) >= 1
    @test pairs[1][1] == 0.0f0

    # normalize_data
    d = Dict("x" => Float32[1.0, 2.0, 3.0])
    nd = IonChannel.normalize_data(d)
    @test abs(mean(nd)) < 1e-6
    @test abs(std(nd) - 1.0f0) < 1e-5
end

@testset "get_specified_datapoints" begin
    x = Float32[1.0, 2.0, 3.0, 4.0]
    y = Float32[0.1, 0.1, 0.1, 0.1]
    data = IonChannel.get_specified_datapoints(x, y, 0.1f0, UInt32(2))
    @test length(data["x"]) == 2
    @test all(t -> t <= 0.2f0, data["dwell times"])
end

@testset "deviation_from_mean_method behavior" begin
    # construct a simple step signal: 10 samples at 0, then 10 at 1
    step = vcat(fill(0.0f0, 10), fill(1.0f0, 10))
    m = IonChannel.MeanDeviationMethod(0.0f0)
    out = IonChannel.deviation_from_mean_method(step, 0.001f0, m)
    @test isa(out, IonChannel.MeanDeviationMethodOutput)
    # expect at least one dwell time detected (transition)
    @test length(out.dwell_times_approx) >= 1
end

@testset "naive and mika methods on synthetic signal" begin
    # synthetic noisy two-state signal: values near 0 and 1 with longer dwell times
    Random.seed!(1234)
    signal = Float32[]
    for i in 1:4
        v = (i % 2 == 0) ? 1.0f0 : 0.0f0
        # longer segments (20 samples) with small noise to guarantee threshold crossings
        append!(signal, v .+ 0.01f0 * randn(Float32, 20))
    end
    Δt = 0.001f0

    # Naive method
    nm = IonChannel.NaiveMethod(UInt16(50))
    nout = IonChannel.naive_method(signal, Δt, nm)
    @test isa(nout, IonChannel.NaiveMethodOutput)

    # Mika method with small number of bins to keep deterministic
    mm = IonChannel.MikaMethod(UInt16(20))
    mout = IonChannel.mika_method(signal, Δt, mm)
    @test isa(mout, IonChannel.MikaMethodOutput)
    @test length(mout.idealized_data) == length(signal)
end
