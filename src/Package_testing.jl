using ModelingToolkit, ThermalSystem_Models, Logging, Revise, CoolProp, Printf, OrdinaryDiffEq
using GraphMakie

Revise.retry()
Logging.disable_logging(Logging.Warn)
#=============================================================#
#                   Preliminaries     
#=============================================================#
begin
    TSM = ThermalSystem_Models
    TSMD = TSM.Dynamics
    TSMT = TSM.DeadTime
    MTK = ModelingToolkit
    Steam = TSMD.Steam
    Gas = TSMD.Gas
    Liq = TSMD.Liq
end

#=============================================================#
#                       TESTING METADATA
# Notes: Didnt work very well
#=============================================================#
    # @named tstpart = Steam.AdiabaticTurbine()
    # @named tstpart2 = Steam.IdealBoiler()
    # @variables t
    # @named copysys = ODESystem(equations(tstpart),t; metadata = Dict(:description => "SteamTurbine"))
    # @named copysys2 = ODESystem(equations(tstpart2),t; metadata = Dict(:description => "boiler"))
    # @named cyss = extend(tstpart, ODESystem(Equation[],t; name = :not, metadata = Dict(:description => "SteamTurbine")))
    # @named dyss = extend(tstpart2, ODESystem(Equation[],t; name = :not, metadata = Dict(:description => "SteamBoiler")))
    # eq = Steam.hydro_connect(cyss.n,dyss.p)
    # @named newsys = ODESystem(eq,t; systems = [cyss,dyss])
    # # substitute(equations(copysys), Steam.hydro_prop_dict)

#=============================================================#
#        OPEN FEEDWATER AND CLOSED FEEDWATER RANKINE LOOP  
#=============================================================#

# Creating and solving system
begin
    Pmax = 150
    Pmid2 = 40
    Pmid1 = 5
    Pmin = 0.1
    Tmax = 600+273.15
    @variables t
    @named ElectricUtil     = Steam.WorkPin()
    @named ElectricGen      = Steam.WorkPin()
    @named ColdUtil         = Steam.HeatTransferPin()
    @named HotUtil          = Steam.HeatTransferPin()

    @named reservoir    = Steam.ContinuityReservoir()
    @named valve        = Steam.SteamFlowSource(ṁ = 1.0)
    @named openfw       = Steam.OpenFeedwaterHeater()
    @named closedfw     = Steam.ClosedFeedwaterHeater()
    @named mixer        = Steam.MixingChamber()

    @named pump1        = Steam.AdiabaticPump(Pout = Pmid1, setpressure = false, η=1.0)
    @named pump2        = Steam.AdiabaticPump(Pout = Pmax, setpressure = false, η=1.0)
    @named pump3        = Steam.AdiabaticPump(Pout = Pmax, setpressure = true, η=1.0, controlinlet = true)

    @named turbine1    = Steam.AdiabaticTurbine(setpressure = true, Pout = Pmid2)
    @named turbine2    = Steam.SIMOAdiabaticTurbine(setpressure = true, Pyin = Pmid1, Pzin = Pmin, ηin = 1.0) 

    @named boiler       = Steam.IdealBoiler(Tout = Tmax)
    @named reheat       = Steam.IdealBoiler(Tout = Tmax)
    @named condensor    = Steam.IdealCondensor()


    connections = vcat(Steam.hydro_connect(reservoir.n,valve.p),
                        Steam.hydro_connect(valve.n, boiler.p),
                        Steam.hydro_connect(boiler.n, turbine1.p),
                        Steam.hydro_connect(turbine1.n, reheat.p, closedfw.p1),
                        Steam.hydro_connect(reheat.n, turbine2.p),
                        Steam.hydro_connect(turbine2.hp.n, openfw.p1),
                        Steam.hydro_connect(turbine2.lp.n, condensor.p),
                        Steam.hydro_connect(condensor.n, pump1.p),
                        Steam.hydro_connect(pump1.n,openfw.p2),
                        Steam.hydro_connect(openfw.n,pump2.p),
                        Steam.hydro_connect(pump2.n,closedfw.p2),
                        Steam.hydro_connect(closedfw.n1,pump3.p),
                        Steam.hydro_connect(pump3.n, mixer.p1),
                        Steam.hydro_connect(closedfw.n2, mixer.p2),
                        Steam.hydro_connect(mixer.n, reservoir.p),
                        TSMD.work_connect(ElectricGen, turbine1.w, turbine2.hp.w, turbine2.lp.w),
                        TSMD.work_connect(ElectricUtil,pump1.w, pump2.w, pump3.w),
                        TSMD.heat_connect(ColdUtil, condensor.q),
                        TSMD.heat_connect(HotUtil, boiler.q, reheat.q));
            

    systems = [reservoir, valve, openfw, closedfw, mixer, pump1, pump2,pump3,turbine1,turbine2, boiler,reheat,condensor, ElectricGen, ElectricUtil,HotUtil,ColdUtil]
    @named odesystem = ODESystem(connections, t; systems = systems)
    TSMD.system_details(odesystem)

    defaults = [valve.Ṁ => 1.0,
        pump1.η => 1.0,
        pump1.P => Pmid2,
        pump2.η => 1.0,
        pump2.P => Pmax,
        pump3.η => 1.0,
        pump3.P => Pmax,
        turbine1.P => Pmid2,
        turbine1.η => 1.0,
        turbine2.Py => Pmid1,
        turbine2.Pz => Pmin,
        turbine2.η => 1.0,
        turbine2.hp.P => Pmid1,
        turbine2.hp.η => 1.0,
        turbine2.lp.P => Pmin,
        turbine2.lp.η => 1.0,
        boiler.T => Tmax,
        reheat.T => Tmax]
        
        simple_sys = structural_simplify(odesystem)
        tspan = (0.0,1.0)
        ode_prob = ODEProblem(simple_sys,[],tspan,defaults)
        sol = solve(ode_prob, Rodas4())
        # Printing solution
        TSMD.showsol(systems,sol)
end

#=============================================================================#
#               Visualization / Plotting
# Notes
#   Cdict = dict(:variable_symbol -> index)
#   vdict = dict(index -> :variable_symbol)
#   epack = dict[ (source_index, destination_index) ] => vairable_symbol
#=============================================================================# 
using Graphs, NetworkLayout, GraphMakie
using GLMakie
using GLMakie.Colors

# Creating Graph Object
begin
    G,cdict,vdict, epack = TSMD.system2graph(odesystem, verbose = false);
    TSMD.add_graph_connection!(odesystem,G,cdict,vdict, ElectricGen.Ẇ ,epack);
    TSMD.add_graph_connection!(odesystem,G,cdict,vdict, ElectricUtil.Ẇ ,epack);
    TSMD.add_graph_connection!(odesystem,G,cdict,vdict, HotUtil.Q̇, epack);
    TSMD.add_graph_connection!(odesystem,G,cdict,vdict, ColdUtil.Q̇, epack);
end

# Creating Graph Object
begin
    edgw_wid_func(s,d,w) = max(w*abs(sol[epack[(s,d)]][end])/10^6,1)
    edgw_wid_func(s,d,w) = 1
    #   GRAPH PLOTTING
    topin = [:HotUtil, :ColdUtil, :ElectricGen, :ElectricUtil]
    locs = [(0.0,5.0), (10.0,5.0), (5.0,10.0), (5.0,0.0)]
    pin_dict = Dict(cdict[topin[i]] => locs[i] for i = 1:4)

    nshape = Dict(i => :rect for i = 1:nv(G))
    ncolo = Dict(i => :lightblue for i =1:nv(G))
    ecolo = Dict(i => :black for i = 1:ne(G))
    nsize = Dict(i => abs.(randn(1))[1]*100 for i = 1:ne(G))

    nshape = Dict(i => :circle for i = 1:nv(G))

    for i in topin--+-
        nshape[cdict[i]] = :circle
        ncolo[cdict[i]] = :yellow
    end

    node_weights = Dict()
    edge_nodes = collect(edges(G))

    for (i,ed) in enumerate(edge_nodes)
        str_end = TSMD.split_variable_sym(epack[(ed.src,ed.dst)])
        if str_end[end] == "Q̇"
            ecolo[i] = :red
        end
        if TSMD.check_for_Ẇ([epack[(ed.src,ed.dst)],epack[(ed.src,ed.dst)]]) == true
            ecolo[i] = :blue
        end
    end
end

# Plot Bezier
begin
    turbine_path = BezierPath([
        MoveTo(Point(0, 0)),
        MoveTo(Point(-0.5,-0.5)),
        LineTo(Point(-0.5, 0.5)),
        LineTo(Point(0.5, 0.75)),
        LineTo(Point(0.5, -0.75)),
        ClosePath()
    ])

    compressor_path = BezierPath([
        MoveTo(Point(0, 0)),
        MoveTo(Point(-0.5,-0.75)),
        LineTo(Point(-0.5, 0.75)),
        LineTo(Point(0.5, 0.5)),
        LineTo(Point(0.5, -0.5)),
        ClosePath()
    ])
end
Plant,sts = TSMD.para_full_system_rankine();
simple_sys = structural_simplify(Plant)
defs =[simple_sys.div_mass_fcn => 50
    simple_sys.breeder_mass_fcn => 50
    simple_sys.cycle_mass_fcn => 30
    simple_sys.loop_mass_fcn => 40
    simple_sys.wall_mass_fcn => 50
    simple_sys.Qdiv => 150e6
    simple_sys.Qbrd => 700e6
    simple_sys.Qwall => 100e6]

tspan = (0.0,10.0)
prob = ODEProblem(simple_sys,[],tspan,defs)
kwargs = (abstol=1e-5, reltol=1e-1)
# sol   = solve(prob, Rodas3();kwargs...);
# copyplant,Plant,sol = TSMD.test_full_system_rankine();

TSMD.showsol()
begin
    allsys = TSMD.systems(copyplant);

    gobj  = Vector{DiGraph}(undef,5)
    cdict = Vector{Dict}(undef,5)
    vdict = Vector{Dict}(undef,5)
    epack = Vector{Dict}(undef,5)

    for i=1:length(allsys)
        gobj[i],cdict[i],vdict[i], epack[i] = TSMD.system2graph(allsys[i], verbose = true);
    end
end

begin
    G,cdict,vdict, epack = TSMD.system2graph(copyplant, verbose = false);
    TSMD.add_graph_connection!(odesystem,G,cdict,vdict, ElectricGen.Ẇ ,epack);
    TSMD.add_graph_connection!(odesystem,G,cdict,vdict, ElectricUtil.Ẇ ,epack);
    TSMD.add_graph_connection!(odesystem,G,cdict,vdict, HotUtil.Q̇, epack);
    TSMD.add_graph_connection!(odesystem,G,cdict,vdict, ColdUtil.Q̇, epack);
end

# Shell(; nlist = [[5,8],]), 
#SFDP(Ptype=Float32, tol=0.01, C=1.0, K=1.5, pin = pin_dict), 
#Spring(;C = 1.2, pin = pin_dict),  #( seed = 2,  C =0.5, Ptype = Float64),
#Stress()

GLMakie.activate!(renderloop = GLMakie.renderloop)

inipos = squaregrid(adjacency_matrix(G); Ptype = Float64, cols = :auto )
len = length(inipos)
order_vec = randn(len)
new_order = sortperm(order_vec)
inipos2 = [inipos[i] for i in new_order]

f,ax,p = graphplot(G; ilabels = [vdict[i] for i = 1:nv(G)], 
                    layout = Stress(), #Stress(; initialpos = inipos2, pin = pin_dict, iterations = 400),
                    arrow_shift=:end,
                    node_marker = :circle,
                    node_attr = (; markersize = 80),
                    curve_distance=.1, curve_distance_usage=true,
                    axis_buffer = 0.1,
                    node_color = [ncolo[i] for i =1:nv(G)],
                    edge_color = [ecolo[i] for i = 1:ne(G)])

f, ax, p = GraphMakie.graphplot(G, edge_width=4)
# gp.node_size = 100

g = wheel_graph(10)
f, ax, p = graphplot(G,
                     edge_width = [2.0 for i in 1:ne(G)],
                     edge_color =  [ecolo[i] for i = 1:ne(G)],
                     arrow_shift=:end,
                     node_size = [100 for i in 1:nv(G)],
                     node_strokewidth = 1,
                     nlabels_align = (:center, :center),
                     node_color = [ncolo[i] for i =1:nv(G)],
                     nlabels = [String(vdict[i]) for i = 1:nv(G)])
begin
    hidedecorations!(ax); hidespines!(ax)
    ax.aspect = DataAspect()


    deregister_interaction!(ax, :rectanglezoom)

    function node_hover_action5(state, idx, event, axis)
        p.node_size[][idx] = state ? 110 : 100
        p.node_size[] = p.node_size[] # trigger observable
    end
    NHOVER = NodeHoverHandler(node_hover_action5)
    register_interaction!(ax, :nhover, NHOVER)



    function node_drag_action5(state, idx, event, axis)
        p[:node_pos][][idx] = event.data
        p[:node_pos][] = p[:node_pos][]
    end

    ndrag = NodeDragHandler(node_drag_action5)
    register_interaction!(ax, :ndrag, ndrag)
end
#   INTERACTIVITY
begin
    deregister_interaction!(ax, :rectanglezoom)
    function node_hover_action(state, idx, event, axis)
        p.node_size[][idx] = state ? 20 : 10
        p.node_size[] = p.node_size[] # trigger observable
    end
    nhover = NodeHoverHandler(node_hover_action)
    function node_drag_action(state, idx, event, axis)
        p[:node_pos][][idx] = event.data
        p[:node_pos][] = p[:node_pos][]
    end
    ndrag = NodeDragHandler(node_drag_action)
    register_interaction!(ax, :nhover, nhover)
    register_interaction!(ax, :ndrag, ndrag)
end
                    # gp = GraphMakie.graphplot(G; ilabels = [vdict[i] for i = 1:nv(G)], 
                    # edgewidth = edgw_wid_func,
                    # nodesize = 0.3,
                    # curves = false,
                    # axis_buffer = 0.1,
                    # nodeshape = [nshape[i] for i =1:nv(G)],
                    # nodecolor = [colo[i] for i =1:nv(G)],
                    # ticks = :auto)sing Makie.Colors
g = wheel_digraph(10)
f, ax, p = graphplot(g, edge_width=4, edge_color=[colorant"black" for i in 1:ne(g)])

g = wheel_graph(10)
f, ax, p = graphplot(g,
                     edge_width = [2.0 for i in 1:ne(g)],
                     edge_color = [colorant"gray" for i in 1:ne(g)],
                     node_size = [10 for i in 1:nv(g)],
                     node_color = [colorant"red" for i in 1:nv(g)])
hidedecorations!(ax); hidespines!(ax)
ax.aspect = DataAspect()

deregister_interaction!(ax, :rectanglezoom)
function node_hover_action(state, idx, event, axis)
    p.node_size[][idx] = state ? 20 : 10
    p.node_size[] = p.node_size[] # trigger observable
end
nhover = NodeHoverHandler(node_hover_action)
register_interaction!(ax, :nhover, nhover)