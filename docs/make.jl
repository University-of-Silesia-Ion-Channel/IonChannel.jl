import Pkg
include("../src/IonChannel.jl")
using .IonChannel

Pkg.activate(joinpath(@__DIR__, "."))
Pkg.instantiate()
using Documenter

makedocs(
    modules=[IonChannel],
    sitename = "IonChannel.jl",
    authors = "Piotr Mika"
)
