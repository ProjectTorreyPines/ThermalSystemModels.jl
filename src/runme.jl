using ThermalSystem_Models
using Test
using Roots, DataFrames, XSteam, Printf
using Graphs, GraphPlot , GraphRecipes
using Plots, GraphRecipes, LinearAlgebra
using  PlotlyJS

TSM = ThermalSystem_Models




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

sys = TSM.rankine_with_loops()
TSM.run_sys(sys; NITER = 5)
df = TSM.eval_sys(sys)
Wnet_cycle = df.Wout[5]-df.Win[5]
η_cycle = Wnet_cycle/df.Qtx[5]
η_total = (df.Wout[6]-df.Win[6])/df.Qh[6]
println("\t η_cyle = $(η_cycle)" )
println(" η_tot = $(η_total)")
display(df)
TSM.component_info(sys.cycles[4])
TSM.show_node_simple(sys.cycles[4])

grp = Vector{Any}(undef,length(sys.cycles))
elabs = Vector{Any}(undef,length(sys.cycles))
nprev = 0;
full_dict = Dict()

gfull,  comp_names, edgelabel_mat, full_dict, comp_dict = TSM.network2graph2(sys.cycles[1]; verbose = false)
offset = ne(gfull)

for idx = 2:5
    gobj,  nms, edgelabel_mat, edgelabel_dict, comp_dict = TSM.network2graph2(sys.cycles[idx]; verbose = false, component_dict = comp_dict)
    source_t = [keyval[1] for keyval ∈ collect(keys(edgelabel_dict))]
    target_t = [keyval[2] for keyval ∈ collect(keys(edgelabel_dict))]
    for i = 1:ne(gobj)
        full_dict[(source_t[i]+offset,target_t[i]+offset)] = edgelabel_dict[(source_t[i],target_t[i])]
    end
    comp_names = vcat(comp_names,nms)
    gfull = blockdiag(gfull,gobj)
    offset = ne(gfull)
end


for hx in sys.heat_exchangers
    hot_element     = hx.hot_stream
    cold_element    = hx.cold_stream

    #initializing src trg
    src = -1
    trg = -1
    Qmag = hot_element.Q
    hot_element.Q > 0 ? dir = 1 : dir = -1

    #first find by name
    hname = hx.hot_stream.name  # strings
    bname = hx.cold_stream.name

    hot_opt = findall(x -> x == hname,comp_names)   # node indexes

    # if more than 1 index, i.e. hot_elements have the same name
    # check the comp_dict to find matches
    if length(hot_opt) > 1
        for opt in hot_opt
            actual_element = comp_dict[opt]
            if TSM.get_inlet_temp(actual_element) == TSM.get_inlet_temp(hot_element)
                # actual element found
                println("Found element $(TSM.get_inlet_temp(actual_element)) and $(TSM.get_inlet_temp(hot_element))")
                src=opt
            end
        end
    else
        src = hot_opt[1]
    end
    cld_opt = findall(x-> x == bname,comp_names)
    if length(cld_opt) > 1
        for opt in cld_opt
            actual_element = comp_dict[opt]
            if TSM.get_inlet_temp(actual_element) == TSM.get_inlet_temp(cld_element)
                # actual element found
                println("Found element $(TSM.get_inlet_temp(actual_element)) and $(TSM.get_inlet_temp(cld_element))")
                src=opt
            end
        end
    else
        trg = cld_opt[1]
    end

    if dir == -1
        tmp = trg
        trg = src
        src = tmp
    end
    add_edge!(gfull,src,trg)
    full_dict[(src,trg)] = "Heat Exchanger Q = $(round(Qmag/10^6))"

end



using NetworkLayout
layout = Spring(Ptype=Float32)
positions = Spring(Ptype = Float64)(adjacency_matrix(gfull))
# typeof(position)
# y = positions[:,2
# display(x)
length(positions)
x = zeros(length(positions),1)
y = zeros(length(positions),1)
for (idx,pos) in enumerate(positions)
    x[idx]=pos[1] * 10
    y[idx]=pos[2] * 10
end
# position[1]


fig2=GraphRecipes.graphplot(gfull, names = comp_names,
    curves=false,
    x=x, y=y,
    nodeshape   =  :rect,
    edgelabel   = full_dict,
    nodecolor = :lightblue,
    nodesize    = 1,
    edge_label_box = false,
    fontsize = 20,
    axis_buffer = 0.0)

    display(fig2)

using PlotlyJS
#                  0             1          2          3                   4               5               6           7
Energy_Sys = ["Fusion Core", "Divertor" , "blanket" , "breeder" ,"Intermediate_loop", "Thermal Cycle", "Cold Utility", "Electric"];

Edict = Dict() # Dict([source,targetr]) => value
FUSION_IDX      = 0
INTER_IDX       = 4
CYCLE_IDX       = 5
ELECTRIC_IDX    = 7
UTILITY_IDX     = 6

# imediate loops
for i = 1:3
    Edict[(FUSION_IDX,i)]     = df.Qh[i]/10^6
    Edict[(ELECTRIC_IDX,i)]   = df.Win[i]/10^6
    Edict[(i,INTER_IDX)]      = -df.Qtx[i]/10^6
    Edict[(i,UTILITY_IDX)]    = df.Qc[i]/10^6
end

Edict[(INTER_IDX,UTILITY_IDX)] = df.Qc[INTER_IDX] / 10^6

# primary cycles
Edict[(INTER_IDX,CYCLE_IDX)]    = df.Qtx[CYCLE_IDX] / 10^6
Edict[(CYCLE_IDX,ELECTRIC_IDX)] = df.Wout[CYCLE_IDX] / 10^6
Edict[(CYCLE_IDX,UTILITY_IDX)]  = df.Qc[CYCLE_IDX] / 10^6
Edict[(CYCLE_IDX,CYCLE_IDX)]    = df.Win[CYCLE_IDX] /10^6
 

dict_k = keys(Edict)
source_t = [keyval[1] for keyval ∈ dict_k]
target_t = [keyval[2] for keyval ∈ dict_k]
vals_t    = [v for v in values(Edict)]

p2 = PlotlyJS.plot(sankey(
    node = attr(
      pad = 15,
      thickness = 20,
      label = Energy_Sys,
    ),
    link = attr(
      source    = source_t, # indices correspond to labels, eg A1, A2, A1, B1, ...
      target    = target_t,
      value     = vals_t,
      label     = vals_t,
    )), 
  Layout(title_text="Basic Sankey Diagram", font_size=10)
)

display(p2)