using Documenter
using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
Pkg.instantiate()

using IonChannel

makedocs(
    modules=[IonChannel],
    sitename = "IonChannel.jl",
    authors = ["IonChannel Contributors"],
    format = Documenter.HTML()
)
