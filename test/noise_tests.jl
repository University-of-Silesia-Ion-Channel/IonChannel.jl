using Test

@testset "noise computation" begin
    data = Float32[1.0, 2.0, 3.0, 4.0]
    ideal = Float32[1.0, 2.0, 2.0, 4.0]

    n = IonChannel.noise(data, ideal)
    @test isa(n, IonChannel.Noise)
    # residuals should be [0,0,1,0]
    @test all(n.ξ .== Float32[0.0, 0.0, 1.0, 0.0])
    @test IonChannel.μ(n) ≈ mean(n.ξ)
    @test IonChannel.σ(n) ≈ std(n.ξ)
    @test IonChannel.noise_data(n) === n.ξ
end
