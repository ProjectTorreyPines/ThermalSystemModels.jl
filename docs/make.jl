using Documenter
using ThermalSystem_Models

makedocs(
    root = joinpath(dirname(pathof(ThermalSystem_Models)), "..", "docs")
    sitename = "ThermalSystem_Models",
    format = Documenter.HTML(),
    modules = [ThermalSystem_Models]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
