using ThermalSystem_Models
using ModelingToolkit, DifferentialEquations
using Test
using MetaGraphs
using Plots
TSM = ThermalSystem_Models
TSMD = TSM.Dynamics
MTK = ModelingToolkit
Steam = TSMD.Steam
Gas = TSMD.Gas
Liq = TSMD.Liq


@testset "ThermalSystem_Models.jl" begin
    TSM = ThermalSystem_Models
    TSMD = TSM.Dynamics
    MTK = ModelingToolkit
    Steam = TSMD.Steam
    Gas = TSMD.Gas
    Liq = TSMD.Liq
    #=========================================================#
    #                Rankine feedwater
    #=========================================================#

    ##########################
    ModelingToolkit.@variables t

    energy_sys, sts, edict = TSMD.default_energy_sys()
    wall_sys, wall_connections, wparams, wdict =
        TSMD.wall_circuit(; load = 100e6, Tmin = 250 + 273.15, Tmax = 450 + 273.15)
    divertor_sys, divertor_connections, dparams, ddict =
        TSMD.divertor_circuit(; load = 150e6, Tmin = 300 + 273.15, Tmax = 650 + 273.15)
    breeder_sys, breeder_connections, bparams, bdict =
        TSMD.breeder_circuit(; load = 250e6, Tmin = 500 + 273.15, Tmax = 800 + 273.15)
    inter_loop_sys, inter_loop_connections, iparams, idict =
        TSMD.intermediate_loop(; Nhx = 4, flowrate = 200, Tmin = 220 + 273.15)
    steam_systems, steam_connections, sparams, sdict = TSMD.feedwater_rankine(; flowrate = 550)
    η_cycle, η_bop = sts

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
    )

    params = vcat(wparams, dparams, bparams, iparams, sparams)
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

    @named hx1 = TSMD.Gen_HeatExchanger(
        B = idict[:inter_loop_hx1],
        A = wdict[:wall_hx],
        returnmode = :eq,
    )
    @named hx2 = TSMD.Gen_HeatExchanger(
        B = idict[:inter_loop_hx2],
        A = ddict[:divertor_hx],
        returnmode = :eq,
    )
    @named hx3 = TSMD.Gen_HeatExchanger(
        B = idict[:inter_loop_hx3],
        A = bdict[:breeder_hx],
        returnmode = :eq,
    )
    @named boilhx = TSMD.Gen_HeatExchanger(
        A = sdict[:steam_boiler],
        B = idict[:inter_loop_hx4],
        returnmode = :eq,
    )

    push!(connections, hx1...)
    push!(connections, hx2...)
    push!(connections, hx3...)
    push!(connections, boilhx...)

    @named sys = ODESystem(connections, t, sts, params; systems = plant_systems)
    utility_vector = [:HotUtility, :ColdUtility, :Electric]
    TSMD.system_details(sys)
    simple_sys = structural_simplify(sys)
    tspan = (0.0, 1.0)
    prob = ODEProblem(simple_sys, [], tspan)
    sol = DifferentialEquations.solve(prob)
    soln(v) = sol[v][end]
    GG = TSMD.system2metagraph(sys, utility_vector; soln = soln, verbose = false)
    # TSMD.showsol(vcat(steam_systems,[sdict[:steam_turbine].lp,sdict[:steam_turbine].hp]),sol)
    println(sol[bdict[:breeder_relief].q.Q̇][end])
    println(sol[wdict[:wall_relief].q.Q̇][end])
    println(sol[ddict[:divertor_relief].q.Q̇][end])
    println(sol[idict[:inter_loop_relief].q.Q̇][end])
    println(sol[idict[:inter_loop_circulator].p.T][end])
    println(sol[idict[:inter_loop_circulator].n.T][end])
    println(sol[idict[:inter_loop_hx1].p.T][end])
    println(sol[idict[:inter_loop_hx2].p.T][end])
    println(sol[idict[:inter_loop_hx3].p.T][end])
    println(sol[idict[:inter_loop_hx3].n.T][end])
    println("cycle eff  $(sol[η_cycle][end])")
    println("plant eff $(sol[η_bop][end])")
    TSMD.reverse_edge!(GG, 12, 21)
    
    gcopy = TSMD.create_plot_graph(GG)
    xLayReqs, vSortReqs, xs, ys, paths, lay2node = TSMD.layers_to_force!(gcopy)
    TSMD.initialize_plot_props!(gcopy, lay2node,xs,ys,paths)
    gc = TSMD.add_plot_elments(gcopy; verbose = false)
    TSMD.set_default_node_prop!(gc, :height, 1.0)

    xLayReqs, vSortReqs, xs, ys, paths, lay2node = TSMD.layers_to_force!(gc)
    x, y = TSMD.setVerticalSpacing!(gc; vspan = 40.0)
    TSMD.setLayerWidth!(gc; pad = 2.5, verbose = true)
    xLayReqs, vSortReqs, xs, ys, paths, lay2node = TSMD.layers_to_force!(gc)
    TSMD.edgeroute_nodes(gc; voff = 0.1)
    TSMD.set_plot_props!(gc)
    gr()
    syslabs = get_prop(gc, :system_labels)
    sysnamedict = Dict([
        "cycle_" => "Brayton Helium",
        "steam_" => "Feedwater Rankine",
        "divertor_" => "Divertor Helium",
        "breeder_" => "Breeder PbLi",
        "inter_loop_" => "Inter. loop Helium",
        "wall_" => "first wall helium",
        "ColdUtility" => "Sink",
        "HotUtility" => "Fusion Core",
    ])

    p = TSMD.plotplant(
        gc;
        numbering = false,
        mode = :path,
        nsize = 2.0,
        compnamesubs = (
            "enfw" => "en\nfw",
            "_" => "\n",
            "circulator" => "pump",
            "hotutility" => "Fusion\nCore",
            "coldutility" => "Sink",
            "relief" => "Trim\ncooler",
            "condensor" => "Cool\nHX",
        ),
        compnameattr = (:right,  4),
        compnamerot = 90,
        sysnamedict = sysnamedict,
        legpad = 0.5,
        legwid =13,
        legheight = 1.5,
        legoffset = 2.0,
        pathattr = (
            linewidth = 1,
            marker = false,
            markersize = 0.0,
            markercolor = :red,
            alpha = 0.7,
            legend = false,
        ),
        figattr = (
            grid = false,
            aspect_ratio = :equal,
            showaxis = false,
            xlim = [-15, 75],
            ylim = [-21, 21],
            xticks = [0, 1, 2, 3, 4, 5, 6, 7],
            plot_title = "Test Feedwater Rankine",
            plot_titlefonthalign = :hcenter,
            plot_titlefontvalign = :bottom,
            dpi = 200,
            plot_titlevspan = .0001,
        ), #aspect_ratio = :equal
    )
    display(p)
        # end?
        # GG = TSMD.system2metagraph(sys,utility_vector; soln = soln, verbose = false);
        # TSMD.showsol(plant_systems,prob)
end
