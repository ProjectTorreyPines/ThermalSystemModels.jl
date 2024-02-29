import MetaGraphs as mg
import MetaGraphsNext as mgn
using Graphs
import ModelingToolkit as MTK
using Printf
import ThermalSystem_Models
import DifferentialEquations as diffeq
tsm = TSM = ThermalSystem_Models
#=============== 
MetaGraphsNext (MGN)
    Constructor - 
    A graph type with custom vertex labels containing vertex-, edge- and graph-level metadata.

    Vertex labels have type Label, 
    vertex (resp. edge, resp. graph) metadata has type VertexData (resp. EdgeData, resp. GraphData). 

    It is recommended not to set Label to an integer type, 
    so as to avoid confusion between vertex labels (which 
    do not change as the graph evolves) and vertex codes (which 
    have type Code<:Integer and can change as the graph evolves).

        MetaGraph{
            Code<:Integer,
            Graph<:AbstractGraph{Code},
            Label,
            VertexData,
            EdgeData,
            GraphData,
            WeightFunction,
            Weight
        } <: AbstractGraph{Code}

        A graph type with custom vertex labels containing vertex-, edge- and graph-level metadata.
        Vertex labels have type Label, while vertex (resp. edge, resp. graph) metadata has type VertexData (resp. EdgeData, resp. GraphData). It is recommended not to set Label to an integer type, so as to avoid confusion between vertex labels (which do not change as the graph evolves) and vertex codes (which have type Code<:Integer and can change as the graph evolves).

        Fields

        graph::Graph: underlying, data-less graph with vertex codes of type Code
        vertex_labels::Dict{Code,Label}: dictionary mapping vertex codes to vertex labels
        vertex_properties::Dict{Label,Tuple{Code,VertexData}}: dictionary mapping vertex labels to vertex codes & metadata
        edge_data::Dict{Tuple{Label,Label},EdgeData}: dictionary mapping edge labels such as (label_u, label_v) to edge metadata
        graph_data::GraphData: metadata for the graph object as a whole
        weight_function::WeightFunction: function computing edge weight from edge metadata, its output must have the same type as default_weight
        default_weight::Weight: default weight used when an edge doesn't exist
=#


#= MGN DEMO =#
begin

    colors = mgn.MetaGraph(
            Graph();  # underlying graph structure
            label_type=Symbol,  # color name
            vertex_data_type=NTuple{3,Int},  # RGB code
            edge_data_type=Symbol,  # result of the addition between two colors
            graph_data="additive colors",  # tag for the whole graph
        )

        # Add Vertices and edges
        @show begin
            # Add nodes
            colors[:red]    = (255, 0, 0);
            colors[:green]  = (0, 255, 0);
            colors[:blue]   = (0, 0, 255);
            colors[Symbol("orange")] = (255, 165, 0)

        
            # Add edges
            colors[:red, :green] = :yellow;
            colors[:red, :blue] = :magenta;
            colors[:green, :blue] = :cyan;
            colors[:green, :blue]
            colors.vertex_properties
        
        end

        collect(mgn.edge_labels(colors))

end

#=
    Creating graph edge and node data structures

Graph
    name
    Plotting
        :reduced_graph  => {7, 8} directed Int64 metagraph with Float64 weights defined by :weight (defau…
        :system_labels  => ["divertor_", "steam_", "breeder_", "inter_loop_", "wall_", "ColdUtility", "Ho…
        :paths          => Dict{SimpleEdge{Int64}, Tuple{Vector{Float64}, Vector{Float64}}}(Edge 9 => 15=…
    
    Technical
        :soln           => soln
        :sys            => ODESystem(0x000000000000054c, Equation[steam_supply₊n₊ṁ(t) + steam_boiler₊p₊ṁ(…
        :utility_vector => [:HotUtility, :ColdUtility, :Electric]
        :util_labels    => ["HotUtility", "ColdUtility", "Electric"]


Vertex
    name
    :nodeType      => :real
    :name          => :divertor_hx
    :sys           => ODESystem(0x00000000000004f2, Equation[0 ~ n₊ṁ(t) + p₊ṁ(t), q₊Q̇(t) ~ p₊Φ(t) + n…
    :parent        => :divertor_

    Plotting
        :pos           => Float32[15.0, 1.77636f-15]
        :height        => 1.0
        :nshape        => :circle
        :normheight    => 0.333333
        :normwidth     => 1
        :displayName   => "hx"
        :markeroutline => :magenta2
        :layer         => 5
        :width         => 1

    Redundant
        :output        => -1
        :input         => -1


Edge
  :evar  => steam_boiler₊p₊Φ(t)
  :etype => :flow
  :mflow => steam_boiler₊p₊ṁ(t)

  :etype      => :transferheat
  :color      => :magenta2
  :isReversed => false
  :evar       => inter_loop_hx2₊Q̇(t)
  :weight     => 1
  :path       => (Float32[15.0, 18.5], Float32[1.77636f-15, 1.77636f-15])
  :sectionId  => 1
  :etype      => :transferheat
  :nSections  => 1

=#
begin
    TSMD  = TSM.Dynamics
    MTK   = ModelingToolkit
    Steam = TSMD.Steam
    Gas   = TSMD.Gas
    Liq   = TSMD.Liq
    TSMfb = TSM.FastBrayton
    #=========================================================#
    #                Rankine feedwater
    #=========================================================#

    ##########################
    ModelingToolkit.@variables t

    energy_sys, sts, edict = TSMD.default_energy_sys();
    wall_sys, wall_connections, wparams, wdict =
        TSMD.wall_circuit(; load = 100e6, Tmin = 250 + 273.15, Tmax = 450 + 273.15);
    divertor_sys, divertor_connections, dparams, ddict =
        TSMD.divertor_circuit(; load = 150e6, Tmin = 300 + 273.15, Tmax = 650 + 273.15);
    breeder_sys, breeder_connections, bparams, bdict =
        TSMD.breeder_circuit(; load = 250e6, Tmin = 500 + 273.15, Tmax = 800 + 273.15);
    inter_loop_sys, inter_loop_connections, iparams, idict =
        TSMD.intermediate_loop(; Nhx = 4, flowrate = 200, Tmin = 220 + 273.15);
    steam_systems, steam_connections, sparams, sdict = TSMD.feedwater_rankine(; flowrate = 550)
    η_cycle, η_bop = sts;

    energy_con = vcat(
        TSMD.work_connect(
            edict[:Electric],
            wdict[:wall_circulator].w,
            ddict[:divertor_circulator].w,
            bdict[:breeder_circulator].w,
            idict[:inter_loop_circulator].w,
            sdict[:steam_hp_pump].w,
            sdict[:steam_lp_pump].w,
            sdict[:steam_turbine].hp.w,
            sdict[:steam_turbine].lp.w,
        ),
        TSMD.heat_connect(
            edict[:HotUtility],
            wdict[:wall_heat].q,
            ddict[:divertor_heat].q,
            bdict[:breeder_heat].q,
        ),
        TSMD.heat_connect(
            edict[:ColdUtility],
            wdict[:wall_relief].q,
            ddict[:divertor_relief].q,
            bdict[:breeder_relief].q,
            idict[:inter_loop_relief].q,
            sdict[:steam_condensor].q,
        ),
        η_cycle ~ 1 - abs(sdict[:steam_condensor].q.Q̇ / sdict[:steam_boiler].q.Q̇),
        η_bop ~ 1 - abs(edict[:ColdUtility].Q̇ / edict[:HotUtility].Q̇),
    );

    params = vcat(wparams, dparams, bparams, iparams, sparams);
    connections = vcat(
        steam_connections,
        inter_loop_connections,
        wall_connections,
        divertor_connections,
        breeder_connections,
        energy_con,
    )
    plant_systems =
        vcat(steam_systems, inter_loop_sys, wall_sys, divertor_sys, breeder_sys, energy_sys)

    mtk.@named hx1 = TSMD.Gen_HeatExchanger(
        B = idict[:inter_loop_hx1],
        A = wdict[:wall_hx],
        returnmode = :eq,
    )
    mtk.@named hx2 = TSMD.Gen_HeatExchanger(
        B = idict[:inter_loop_hx2],
        A = ddict[:divertor_hx],
        returnmode = :eq,
    )
    mtk.@named hx3 = TSMD.Gen_HeatExchanger(
        B = idict[:inter_loop_hx3],
        A = bdict[:breeder_hx],
        returnmode = :eq,
    )
    mtk.@named boilhx = TSMD.Gen_HeatExchanger(
        A = sdict[:steam_boiler],
        B = idict[:inter_loop_hx4],
        returnmode = :eq,
    )

    push!(connections, hx1...)
    push!(connections, hx2...)
    push!(connections, hx3...)
    push!(connections, boilhx...)

    mtk.@named sys = mtk.ODESystem(connections, t, sts, params; systems = plant_systems)
    utility_vector = [:HotUtility, :ColdUtility, :Electric]
    TSMD.system_details(sys)
    simple_sys = mtk.structural_simplify(sys)
    tspan = (0.0, 1.0)
    prob = mtk.ODEProblem(simple_sys, [], tspan)
    sol = diffeq.solve(prob, diffeq.ImplicitEuler())
    soln(v) = sol[v][end]

    
    utility_vector = [:HotUtility, :ColdUtility, :Electric];
    G = TSMD.system2metagraph(sys, utility_vector; soln = soln, verbose = false);
    
    para_vars = MTK.parameters(simple_sys);
    para_syms = TSMD.variable2symbol(MTK.parameters(simple_sys))
    para_vals = prob.p;
    gcopy = TSMD.create_plot_graph(G; toignore = [:steam_condensor], verbose = false);
    xLayReqs, vSortReqs, xs, ys, paths, lay2node = TSMD.layers_to_force!(gcopy);
    TSMD.initialize_plot_props!(gcopy, lay2node,xs,ys,paths);
    TSMD.add_plot_elments!(gcopy);
    TSMD.set_default_node_prop!(gcopy, :height, 1.0);
    xLayReqs, vSortReqs, xs, ys, paths, lay2node = TSMD.layers_to_force!(gcopy);
    x, y = TSMD.setVerticalSpacing!(gcopy; vspan = 40.0);
    TSMD.setLayerWidth!(gcopy; pad = 2.5, verbose = false);
    xLayReqs, vSortReqs, xs, ys, paths, lay2node = TSMD.layers_to_force!(gcopy);
    TSMD.edgeroute_nodes(gcopy; voff = 0.1);
    TSMD.set_plot_props!(gcopy);
    gcopy;
end