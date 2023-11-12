include("docutils.jl")


# Reading Doc Strings
complibs = Dict("multiphase" =>      "src/ODE_Systems/01-MultiPhase.jl",
                "incompressible" =>  "src/ODE_Systems/01-Incompressible.jl",
                "gaseous" =>         "src/ODE_Systems/01-ThermoGas.jl")

for k in collect(keys(complibs))
    generate_libtxt(k,library_docstrings(complibs[k]))
end

# all text files in the pre folder
tomake = [sp[1] for sp in split.(readdir("docs/pre"),".")]
for tm in tomake
    buildhtml(tm)
end
