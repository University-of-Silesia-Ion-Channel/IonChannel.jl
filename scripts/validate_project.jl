# validate_project.jl
# Simple project validation script for CI
import Pkg
Pkg.instantiate()
println("Project instantiated successfully.")
# verify Project.toml contains a name and version
using TOML
proj = TOML.parsefile("Project.toml")
if !haskey(proj, "name") || !haskey(proj, "version")
    error("Project.toml missing required keys 'name' or 'version'")
end
println("Project.toml has name=$(proj["name"]) version=$(proj["version"]).")
