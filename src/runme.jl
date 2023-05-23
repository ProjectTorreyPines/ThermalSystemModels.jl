# using ThermalSystem_Models
# using Test
using Roots, DataFrames, XSteam, Printf
using Graphs, GraphPlot , GraphRecipes
using Plots, GraphRecipes, LinearAlgebra
using  PlotlyJS

using ModelingToolkit, ThermalSystem_Models, Logging, Revise, CoolProp, Printf, OrdinaryDiffEq, Plots
Logging.disable_logging(Logging.Warn)

TSM = ThermalSystem_Models
TSMD = TSM.Dynamics
TSMT = TSM.DeadTime
MTK = ModelingToolkit
Steam = TSMD.Steam
Gas = TSMD.Gas
Liq = TSMD.Liq
# cycle_vec   
#   1    divertor_circuit     
#   2    blanket_circuit
#   3    breeder_circuit
#   4    extraction_network,
#   5    rankine

# hx_vec  
#   1    divertor_heat_exchanger
#   2    blanket_heat_exchanger
#   3    breeder_heat_exchanger
#   4    boiler_heat_exchanger

function init_test()
    sys = TSM.rankine_with_loops()
    TSM.run_sys(sys; NITER = 5)
    gdawg, Edict, Cdict, name_list = TSM.network2graph2(sys)
    return gdawg
end

# g,ed,cd,nd,f= TSM.network2graph2(sys)
# display(fig2)
# using PlotlyJS
# #                  0             1          2          3                   4               5               6           7
# Energy_Sys = ["Fusion Core", "Divertor" , "blanket" , "breeder" ,"Intermediate_loop", "Thermal Cycle", "Cold Utility", "Electric"];
# Edict = Dict() # Dict([source,targetr]) => value
# FUSION_IDX      = 0
# INTER_IDX       = 4
# CYCLE_IDX       = 5
# ELECTRIC_IDX    = 7
# UTILITY_IDX     = 6

# # imediate loops
# for i = 1:3
#     Edict[(FUSION_IDX,i)]     = df.Qh[i]/10^6
#     Edict[(ELECTRIC_IDX,i)]   = df.Win[i]/10^6
#     Edict[(i,INTER_IDX)]      = -df.Qtx[i]/10^6
#     Edict[(i,UTILITY_IDX)]    = df.Qc[i]/10^6
# end
# Edi
# # primary cycles
# Edict[(INTER_IDX,CYCLE_IDX)]    = df.Qtx[CYCLE_IDX] / 10^6
# Edict[(CYCLE_IDX,ELECTRIC_IDX)] = df.Wout[CYCLE_IDX] / 10^6
# Edict[(CYCLE_IDX,UTILITY_IDX)]  = df.Qc[CYCLE_IDX] / 10^6
# Edict[(CYCLE_IDX,CYCLE_IDX)]    = df.Win[CYCLE_IDX] /10^6
# dict_k = keys(Edict)
# source_t = [keyval[1] for keyval ∈ dict_k]
# target_t = [keyval[2] for keyval ∈ dict_k]
# vals_t    = [v for v in values(Edict)]
# p2 = PlotlyJS.plot(sankey(
#     node = attr(
#       pad = 15,
#       thickness = 20,
#       label = Energy_Sys,
#     ),
#     link = attr(
#       source    = source_t, # indices correspond to labels, eg A1, A2, A1, B1, ...
#       target    = target_t,
#       value     = vals_t,
#       label     = vals_t,
#     )), 
#   Layout(title_text="Basic Sankey Diagram", font_size=10)
# )
# display(p2)
# display(fig2)