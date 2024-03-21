using ThermalSystemModels
using ModelingToolkit
using Test
using MetaGraphs
using Plots, DifferentialEquations

display("Running Tests")

@testset "ThermalSystemModels.jl" begin
    TSM   = ThermalSystemModels
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
    sol = DifferentialEquations.solve(prob, DifferentialEquations.ImplicitEuler())

    println("Solved Rankine")
    
end
