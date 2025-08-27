import Pkg
Pkg.activate(".")

using Test

include("../src/IonChannel.jl")
using .IonChannel

println("Running IonChannel tests...")
Test.@testset "IonChannel tests" begin
    include("histogram_tests.jl")
    include("noise_tests.jl")
    include("more_tests.jl")
    include("deep_channel_tests.jl")
end
