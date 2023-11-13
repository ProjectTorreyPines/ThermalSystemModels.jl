include("docutils.jl")


# Reading Doc Strings
complibs = Dict("MULTIPHASE" =>      "src/ODE_Systems/01-MultiPhase.jl",
                "INCOMPRESSIBLE" =>  "src/ODE_Systems/01-Incompressible.jl",
                "GASEOUS" =>         "src/ODE_Systems/01-ThermoGas.jl")

for k in collect(keys(complibs))
    generate_libtxt(k,complibs[k])
end
ulibs = Dict("UTILITIES" =>      "src/ODE_Systems/03-MTK_UTILS.jl",
                "OTHER" =>          "src/ODE_Systems/Dynamics.jl")
for k in collect(keys(ulibs))
    generate_utiltxt(k,ulibs[k])
end

# all text files in the pre folder
tomake = [sp[1] for sp in split.(readdir("docs/pre"),".")]
for tm in tomake
    buildhtml(tm)
end
