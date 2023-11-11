using Documenter
using ThermalSystem_Models

makedocs(
    root = joinpath(dirname(pathof(ThermalSystem_Models)), "..", "docs"),
    source = "src",
    checkdocs = :all,
    sitename = "ThermalSystemModels",
    format = Documenter.HTML(; prettyurls=false, sidebar_sitename=false,),
    modules = [ThermalSystem_Models, ThermalSystem_Models.Dynamics,ThermalSystem_Models.Dynamics.Gas,ThermalSystem_Models.Dynamics.Liq,ThermalSystem_Models.Dynamics.Steam],
    pages = [
        "Overview" => "index.md",
        "tutorial.md",
        "Tutorial" => [
            "Introduction" => "introduction.md",
            "Examples" => "examples.md",
        ],
        "blocks.md",
        "Blocks" => [
            "Multiphase Library" => "blocks_mf.md",
            "Incompressible Liquid Library" => "blocks_il.md",
            "Gas Library" => "blocks_gas.md",
        ],
        "resources.md",
        "Resources" => [
            "Thermodynamics" => "theory_thermo_basics.md",
            "Gas systems" => "theory_gas.md",
            "Multiphase and Liquid systems" => "theory_multiphase.md",
            "Fast Brayton" => "fast_brayton.md",
        ],
        "apireference.md",
        "gallery.md",
    ]
    )
# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
