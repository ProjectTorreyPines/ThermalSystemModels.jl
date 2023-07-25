using ThermalSystem_Models
using ModelingToolkit
using Test

TSM = ThermalSystem_Models
TSMD = TSM.Dynamics
MTK = ModelingToolkit
Steam = TSMD.Steam
Gas = TSMD.Gas
Liq = TSMD.Liq


@testset "ThermalSystem_Models.jl" begin
#=========================================================#
#                Rankine feedwater
#=========================================================#
    
    ##########################
    function node_prop(G,prop)
        return [get_prop(G,i,prop) for i =1:nv(G)]
    end
    # using HiGHS
    # Logging.disable_logging(Logging.Warn)
    ModelingToolkit.@variables t
    # function default_energy_sys()
    #     sts = ModelingToolkit.@variables η_cycle(t)=0.5 η_bop(t)=0.5
    #     @named Electric      = Gas.WorkPin()                # Connects to all pumps, turbines, compressors
    #     @named HotUtility    = Gas.HeatTransferPin()
    #     @named ColdUtility   = Gas.HeatTransferPin()
    #     energy_sys = [Electric,HotUtility,ColdUtility]

    #     sysdict = TSMD.sys2dict(energy_sys)

    #     return energy_sys,sts, sysdict
    # end

    # function wall_circuit(; max_pressure = 80, pressrue_drop = 5, Tmin = 450 +273.15, Tmax = 550+273.15)
    #     params = @parameters Qwall = 100e6
    #     pressure_max_wall   = max_pressure;  # bar
    #     pressure_drop_wall  = pressrue_drop;
    #     pressure_min_wall = pressure_max_wall-pressure_drop_wall;
    #     Tmin_wall = Tmin;
    #     Tmax_wall = Tmax;
        
    #     @named wall_supply          = Gas.SinglePortReservoir(P = pressure_min_wall, T = Tmin_wall);
    #     @named wall_circulator      = Gas.PassiveThermoCompressor2(η = 0.9);
    #     @named wall_const_pressure  = Gas.SetPressure(P=pressure_max_wall);
    #     @named wall_heat            = Gas.FlowControlThermoHeatTransfer(ΔP = pressure_drop_wall,Tout = Tmax_wall);
    #     @named wall_hx              = Gas.ThermoHeatTransfer() ; 
    #     @named wall_relief          = Gas.ReliefElement();
        
    #     wall_sys = [wall_supply,wall_circulator,wall_const_pressure,wall_heat,wall_hx,wall_relief];

    #     sysdict = TSMD.sys2dict(wall_sys)

    #     wall_connections = vcat(Gas.gas_connect(wall_supply.n,wall_circulator.p,wall_relief.n),
    #                             Gas.gas_connect(wall_circulator.n,wall_const_pressure.p),
    #                             Gas.gas_connect(wall_const_pressure.n,wall_heat.p),
    #                             Gas.gas_connect(wall_heat.n,wall_hx.p),
    #                             Gas.gas_connect(wall_hx.n,wall_relief.p),
    #                             wall_heat.Q̇ ~ Qwall);
    #     return wall_sys, wall_connections, params, sysdict
    # end

    # function divertor_circuit(; max_pressure = 80, pressrue_drop = 5, Tmin = 450 +273.15, Tmax = 550+273.15)
    #     params = @parameters Qdivertor = 100e6
    #     pressure_max_divertor = max_pressure;  # bar
    #     pressure_drop_divertor = pressrue_drop;
    #     pressure_min_divertor = pressure_max_divertor-pressure_drop_divertor;
    #     Tmin_divertor = Tmin;
    #     Tmax_divertor = Tmax;

    #     @named divertor_supply          = Gas.SinglePortReservoir(P = pressure_min_divertor, T = Tmin_divertor)
    #     @named divertor_circulator      = Gas.PassiveThermoCompressor2(η = 0.9)
    #     @named divertor_const_pressure  = Gas.SetPressure(P=pressure_max_divertor)
    #     @named divertor_heat            = Gas.FlowControlThermoHeatTransfer(ΔP = pressure_drop_divertor,Tout = Tmax_divertor)
    #     @named divertor_hx              = Gas.ThermoHeatTransfer()  
    #     @named divertor_relief          = Gas.ReliefElement()
    #     divertor_sys = [divertor_supply,divertor_circulator,divertor_const_pressure,divertor_heat,divertor_hx,divertor_relief]
    #     sysdict = TSMD.sys2dict(divertor_sys)
    #     divertor_connections = vcat(Gas.gas_connect(divertor_supply.n,divertor_circulator.p,divertor_relief.n),
    #                             Gas.gas_connect(divertor_circulator.n,divertor_const_pressure.p),
    #                             Gas.gas_connect(divertor_const_pressure.n,divertor_heat.p),
    #                             Gas.gas_connect(divertor_heat.n,divertor_hx.p),
    #                             Gas.gas_connect(divertor_hx.n,divertor_relief.p),
    #                             divertor_heat.Q̇ ~ Qdivertor)

    #     return divertor_sys, divertor_connections, params, sysdict
    # end

    # function breeder_circuit(; max_pressure = 40, pressrue_drop = 8, Tmin = 750 +273.15, Tmax = 900+273.15)
    #     params = @parameters Qbreeder = 100e6
    #     pressure_max_breeder = max_pressure;  # bar
    #     pressure_drop_breeder = pressrue_drop;
    #     pressure_min_breeder = pressure_max_breeder-pressure_drop_breeder;
    #     Tmin_breeder = Tmin;
    #     Tmax_breeder = Tmax;
        
    #     @named breeder_supply           = Liq.SinglePortReservoir(P = pressure_min_breeder,T = Tmin_breeder)
    #     @named breeder_circulator       = Liq.PassiveIncompressiblePump2(η = 0.9)
    #     @named breeder_const_pressure   = Liq.SetPressure2(P=pressure_max_breeder)
    #     @named breeder_heat       		= Liq.FlowControlIncompressibleHeatTransfer(ΔP = pressure_drop_breeder, Tout = Tmax_breeder)
    #     @named breeder_hx               = Liq.IncompressibleHeatTransfer()
    #     @named breeder_relief           = Liq.ReliefElement()
        
    #     breeder_sys = [breeder_supply,breeder_circulator,breeder_const_pressure,breeder_heat,breeder_hx,breeder_relief]
    #     sysdict = TSMD.sys2dict(breeder_sys)
    #     breeder_connections = vcat(Liq.incompressible_connect(breeder_supply.n,breeder_circulator.p,breeder_relief.n),
    #                             Liq.incompressible_connect(breeder_circulator.n,breeder_const_pressure.p),
    #                             Liq.incompressible_connect(breeder_const_pressure.n,breeder_heat.p),
    #                             Liq.incompressible_connect(breeder_heat.n,breeder_hx.p),
    #                             Liq.incompressible_connect(breeder_hx.n,breeder_relief.p),
    #                             breeder_heat.Q̇ ~ Qbreeder)


    #     return breeder_sys, breeder_connections, params, sysdict
    # end

    # function feedwater_rankine(;max_pressure = 150,mid_pressure = 12, min_pressure = 0.1, flowrate = 50)
    #     params = @parameters steam_ṁ = flowrate
    #     steam_pmid = mid_pressure
    #     steam_pmin = min_pressure
    #     steam_pmax = max_pressure
    #     @named steam_supply         = Steam.ContinuityReservoir()
    #     @named steam_hp_pump        = Steam.AdiabaticPump(Pout = steam_pmax,setpressure = true, η=1.0)
    #     @named steam_boiler         = Steam.SteamHeatTransfer()
    #     @named steam_turbine        = Steam.SIMOAdiabaticTurbine(setpressure = true, Pyin = steam_pmid, Pzin = steam_pmin,ηin = 1.0) 
    #     @named steam_lp_pump        = Steam.AdiabaticPump(Pout = steam_pmid, setpressure = false, η=1.0)
    #     @named steam_condensor      = Steam.IdealCondensor()
    #     @named steam_openfw         = Steam.OpenFeedwaterHeater()

    #     steam_connections = vcat(Steam.hydro_connect(steam_supply.n, steam_boiler.p),
    #                         Steam.hydro_connect(steam_boiler.n, steam_turbine.p),
    #                         Steam.hydro_connect(steam_turbine.hp.n, steam_openfw.p1),
    #                         Steam.hydro_connect(steam_turbine.lp.n, steam_condensor.p),
    #                         Steam.hydro_connect(steam_condensor.n, steam_lp_pump.p),
    #                         Steam.hydro_connect(steam_lp_pump.n,steam_openfw.p2),
    #                         Steam.hydro_connect(steam_openfw.n,steam_hp_pump.p),
    #                         Steam.hydro_connect(steam_hp_pump.n,steam_supply.p),
    #                         steam_boiler.n.ṁ ~ steam_ṁ)

    #     steam_systems = [steam_boiler,steam_turbine,steam_lp_pump,steam_condensor,steam_openfw,steam_hp_pump, steam_supply];
    #     sysdict = TSMD.sys2dict(steam_systems)
    #     return steam_systems, steam_connections, params,  sysdict
    # end

    # function intermediate_loop(;Pmax = 40 ,Pmin = 32, Nhx = 3, Tmin = 350 + 273.15, )
    #     params = @parameters inter_loop_ṁ = 100
    #     pressure_max_loop = 40;  # bar
    #     pressure_drop_loop = 8;
    #     pressure_min_loop = pressure_max_loop-pressure_drop_loop;
    #     Tmin_loop = 350+273.15;

    #     pdrop_per = pressure_drop_loop / Nhx

    #     @named inter_loop_supply            = Gas.SinglePortReservoir(P = pressure_min_loop, T = Tmin_loop);
    #     @named inter_loop_circulator        = Gas.PassiveThermoCompressor2(η = 0.9);
    #     @named inter_loop_const_pressure    = Gas.SetPressure(P=pressure_max_loop);
    #     @named inter_loop_hx1		    	= Gas.ThermoHeatTransfer(ΔP = pdrop_per) ;
    #     @named inter_loop_relief          	= Gas.ReliefElement();

    #     inter_loop_sys = [inter_loop_supply,inter_loop_relief,inter_loop_circulator,inter_loop_const_pressure,inter_loop_hx1]
        
    #     inter_loop_connections = vcat(Gas.gas_connect(inter_loop_supply.n, inter_loop_circulator.p,inter_loop_relief.n),
    #                                     Gas.gas_connect(inter_loop_circulator.n,inter_loop_const_pressure.p),
    #                                     Gas.gas_connect(inter_loop_const_pressure.n,inter_loop_hx1.p),
    #                                     inter_loop_circulator.p.ṁ ~ inter_loop_ṁ);
    #     for i = 2:Nhx
    #         push!(inter_loop_sys,Gas.ThermoHeatTransfer(name = Symbol("inter_loop_hx" * "$(i)") ,ΔP = pdrop_per));
    #         inter_loop_connections = vcat(inter_loop_connections, Gas.gas_connect(inter_loop_sys[end-1].n, inter_loop_sys[end].p))
    #     end
    #     inter_loop_connections = vcat(inter_loop_connections, 
    #                                 Gas.gas_connect(inter_loop_sys[end].n, inter_loop_relief.p))

    #     sysdict = TSMD.sys2dict(inter_loop_sys)
    #     return inter_loop_sys, inter_loop_connections, params, sysdict
    # end

    # function rem_all_edge!(G,nds; ignore = [],verbose = false)
    #     remd = 0
    #     for na in nds
    #         e_na = all_neighbors(G,na)
    #         for ena in e_na
    #             if ena ∈ ignore
    #                 continue
    #             end
    #             verbose ? println("$(e_na) is not in $(ignore)") : nothing

    #             if has_edge(G,na,ena) 
    #                 rem_edge!(G,na,ena)
    #                 remd += 1
    #             end
    #             if has_edge(G,ena,na) 
    #                 rem_edge!(G,ena,na)
    #                 remd += 1
    #             end
    #         end
    #     end
    #     verbose ? println("removed $(remd) edges") : nothing
    #     return G
    # end


    energy_sys,sts, edict = TSMD.default_energy_sys();
    wall_sys, wall_connections, wparams, wdict = TSMD.wall_circuit();
    divertor_sys, divertor_connections, dparams, ddict = TSMD.divertor_circuit();
    breeder_sys, breeder_connections, bparams, bdict = TSMD.breeder_circuit();
    inter_loop_sys, inter_loop_connections, iparams, idict = TSMD.intermediate_loop(;Nhx = 4);
    steam_systems, steam_connections, sparams,  sdict = TSMD.feedwater_rankine();
    η_cycle, η_bop = sts

    energy_con = vcat(TSMD.work_connect(edict[:Electric], wdict[:wall_circulator].w, ddict[:divertor_circulator].w, 
                                        bdict[:breeder_circulator].w, idict[:inter_loop_circulator].w, 
                                        sdict[:steam_hp_pump].w,sdict[:steam_lp_pump].w, sdict[:steam_turbine].hp.w, 
                                        sdict[:steam_turbine].lp.w),
                    TSMD.heat_connect(edict[:HotUtility], wdict[:wall_heat].q,ddict[:divertor_heat].q, bdict[:breeder_heat].q),
                    TSMD.heat_connect(edict[:ColdUtility], wdict[:wall_relief].q, ddict[:divertor_relief].q, bdict[:breeder_relief].q, idict[:inter_loop_relief].q, sdict[:steam_condensor].q),
                    η_cycle ~ 1 - abs(sdict[:steam_condensor].q.Q̇/sdict[:steam_boiler].q.Q̇),
                    η_bop ~ 1 - abs(edict[:ColdUtility].Q̇/edict[:HotUtility].Q̇));
    # 

    params = vcat(wparams,dparams,bparams,iparams,sparams)

    connections = vcat(steam_connections,inter_loop_connections,wall_connections, divertor_connections, breeder_connections,energy_con);

    plant_systems = vcat(steam_systems,inter_loop_sys,wall_sys, divertor_sys, breeder_sys, energy_sys);

    @named hx1      = TSMD.Gen_HeatExchanger(B = idict[:inter_loop_hx1], A = wdict[:wall_hx], returnmode =:eq);
    @named hx2      = TSMD.Gen_HeatExchanger(B = idict[:inter_loop_hx2], A = ddict[:divertor_hx], returnmode =:eq);
    @named hx3      = TSMD.Gen_HeatExchanger(B = idict[:inter_loop_hx3], A = bdict[:breeder_hx], returnmode =:eq);
    @named boilhx   = TSMD.S2G_HeatExchanger(A = sdict[:steam_boiler], B = idict[:inter_loop_hx4], returnmode =:eq);

    push!(connections,hx1...);
    push!(connections,hx2...);
    push!(connections,hx3...);
    push!(connections,boilhx...);

    @named sys = ODESystem(connections, t,sts,params; systems = plant_systems);
    utility_vector = [:HotUtility, :ColdUtility, :Electric];
    TSMD.system_details(sys)
    simple_sys = structural_simplify(sys)
    tspan = (0.0,1.0)
    prob = ODEProblem(simple_sys,[],tspan)
    soln(v) = prob.f.observed(v, prob.u0, prob.p,0.0)
    GG = TSMD.system2metagraph(sys,utility_vector; soln = soln, verbose = false);
end
