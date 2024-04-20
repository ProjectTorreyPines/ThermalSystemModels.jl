module Dynamics
using ModelingToolkit, Plots, Revise, Unitful, Logging,  Printf
using OrdinaryDiffEq, NonlinearSolve, DifferentialEquations
using Symbolics
using LayeredLayouts, MetaGraphs, Graphs, Plots, Random
using Statistics, GeometryBasics
ModelingToolkit.@variables t
# Logging.disable_logging(Logging.Warn) NonlinearSolve, DifferentialEquations, 

include("03-MTK_UTILS.jl")
module Gas
    include("01-ThermoGas.jl")
end

module Steam
    include("01-MultiPhase.jl")
end

module Liq
    include("01-Incompressible.jl")
end

"""
    L2G_HeatExchanger(; name, ϵ = 0.95, A = Liq.IncompressibleHeatTransfer(), B = Gas.ThermoHeatTransfer())

DOCSTRING
Liquid to Gas heat exchanger

OPTIONAL INPUTS
- `name`: DESCRIPTION
- `ϵ = 0.95`: DESCRIPTION
- `A = Liq.IncompressibleHeatTransfer()`: DESCRIPTION
- `B = Gas.ThermoHeatTransfer()`: DESCRIPTION
"""
function L2G_HeatExchanger(;
    name,
    ϵ = 0.95,
    A = Liq.IncompressibleHeatTransfer(),
    B = Gas.ThermoHeatTransfer(),
)

    ps = @parameters ϵ = ϵ

    @variables Q̇(t) = 0.0 Cmin(t) = 0.0

    eqs = [
        Q̇ ~ ϵ * (A.p.T - B.p.T) * ((A.C < B.C) * A.C + (A.C >= B.C) * B.C)      #heat transfer out of A -> B , if A/T > B/T 
        0 ~ A.Q̇ + B.Q̇
        A.Q̇ ~ -Q̇
    ]
    ODESystem(eqs, t, [Q̇], ps; name = name, systems = [A, B], defaults = [ϵ => 0.95])
end


"""
    S2G_HeatExchanger(; name, ϵ = 0.95, A = Steam.SteamHeatTransfer(), B = Gas.ThermoHeatTransfer(), returnmode = :eqs)

DOCSTRING
Steam to Gas heat exchanger

OPTIONAL INPUTS
- `name`: DESCRIPTION
- `ϵ = 0.95`: DESCRIPTION
- `A = Steam.SteamHeatTransfer()`: DESCRIPTION
- `B = Gas.ThermoHeatTransfer()`: DESCRIPTION
- `returnmode = :eqs`: DESCRIPTION
"""
function S2G_HeatExchanger(;
    name,
    ϵ = 0.95,
    A = Steam.SteamHeatTransfer(),
    B = Gas.ThermoHeatTransfer(),
    returnmode = :eqs,
)
    if returnmode == :ode

        ps = @parameters ϵ = ϵ

        eqs = [                                                 # max temperature change by gas, + if heat from gas to steam
            A.Q̇ ~
                (
                    (
                        A.p.ṁ * (
                            Steam.stm_hptfunc(A.p.P, A.p.T + ((B.p.T - A.p.T) * ϵ)) - A.p.h
                        )
                    ) <= (B.C * ((B.p.T - A.p.T) * ϵ))
                ) * (
                    A.p.ṁ *
                    (Steam.stm_hptfunc(A.p.P, A.p.T + ((B.p.T - A.p.T) * ϵ)) - A.p.h)
                ) +
                (
                    (B.C * ((B.p.T - A.p.T) * ϵ)) <= (
                        A.p.ṁ * (
                            Steam.stm_hptfunc(A.p.P, A.p.T + ((B.p.T - A.p.T) * ϵ)) - A.p.h
                        )
                    )
                ) * (B.C * ((B.p.T - A.p.T) * ϵ))     # heat transfer out of A -> B , if A/T > B/T 
            0 ~ A.Q̇ + B.Q̇
        ]
        ODESystem(eqs, t, [], ps; name = name, systems = [A, B])
    else
        eqs = [                                                 # max temperature change by gas, + if heat from gas to steam
            A.Q̇ ~
                (
                    (
                        A.p.ṁ * (
                            Steam.stm_hptfunc(A.p.P, A.p.T + ((B.p.T - A.p.T) * ϵ)) - A.p.h
                        )
                    ) <= (B.C * ((B.p.T - A.p.T) * ϵ))
                ) * (
                    A.p.ṁ *
                    (Steam.stm_hptfunc(A.p.P, A.p.T + ((B.p.T - A.p.T) * ϵ)) - A.p.h)
                ) +
                (
                    (B.C * ((B.p.T - A.p.T) * ϵ)) <= (
                        A.p.ṁ * (
                            Steam.stm_hptfunc(A.p.P, A.p.T + ((B.p.T - A.p.T) * ϵ)) - A.p.h
                        )
                    )
                ) * (B.C * ((B.p.T - A.p.T) * ϵ))     # heat transfer out of A -> B , if A/T > B/T 
            0 ~ A.Q̇ + B.Q̇
        ]
        return eqs
    end
end


"""
    L2S_HeatExchanger(; name, ϵ = 0.95, A = Steam.SteamHeatTransfer(), B = Liq.IncompressibleHeatTransfer(), returnmode = :eqs)

DOCSTRING
Liquid to steam heat exchanger

OPTIONAL INPUTS
- `name`: DESCRIPTION
- `ϵ = 0.95`: DESCRIPTION
- `A = Steam.SteamHeatTransfer()`: DESCRIPTION
- `B = Liq.IncompressibleHeatTransfer()`: DESCRIPTION
- `returnmode = :eqs`: DESCRIPTION
"""
function L2S_HeatExchanger(;
    name,
    ϵ = 0.95,
    A = Steam.SteamHeatTransfer(),
    B = Liq.IncompressibleHeatTransfer(),
    returnmode = :eqs,
)

    if returnmode == :ode

        ps = @parameters ϵ = ϵ

        eqs = [                                                 # max temperature change by gas, + if heat from gas to steam
            A.Q̇ ~
                (
                    (
                        A.p.ṁ * (
                            Steam.stm_hptfunc(A.p.P, A.p.T + ((B.p.T - A.p.T) * ϵ)) - A.p.h
                        )
                    ) <= (B.C * ((B.p.T - A.p.T) * ϵ))
                ) * (
                    A.p.ṁ *
                    (Steam.stm_hptfunc(A.p.P, A.p.T + ((B.p.T - A.p.T) * ϵ)) - A.p.h)
                ) +
                (
                    (B.C * ((B.p.T - A.p.T) * ϵ)) <= (
                        A.p.ṁ * (
                            Steam.stm_hptfunc(A.p.P, A.p.T + ((B.p.T - A.p.T) * ϵ)) - A.p.h
                        )
                    )
                ) * (B.C * ((B.p.T - A.p.T) * ϵ))     # heat transfer out of A -> B , if A/T > B/T 
            0 ~ A.Q̇ + B.Q̇
        ]
        ODESystem(eqs, t, [], ps; name = name, systems = [A, B])
    else
        eqs = [                                                 # max temperature change by gas, + if heat from gas to steam
            A.Q̇ ~
                (
                    (
                        A.p.ṁ * (
                            Steam.stm_hptfunc(A.p.P, A.p.T + ((B.p.T - A.p.T) * ϵ)) - A.p.h
                        )
                    ) <= (B.C * ((B.p.T - A.p.T) * ϵ))
                ) * (
                    A.p.ṁ *
                    (Steam.stm_hptfunc(A.p.P, A.p.T + ((B.p.T - A.p.T) * ϵ)) - A.p.h)
                ) +
                (
                    (B.C * ((B.p.T - A.p.T) * ϵ)) <= (
                        A.p.ṁ * (
                            Steam.stm_hptfunc(A.p.P, A.p.T + ((B.p.T - A.p.T) * ϵ)) - A.p.h
                        )
                    )
                ) * (B.C * ((B.p.T - A.p.T) * ϵ))     # heat transfer out of A -> B , if A/T > B/T 
            0 ~ A.Q̇ + B.Q̇
        ]
        return eqs
    end
end


"""
    Gen_HeatExchanger(; name, ϵ = 0.95, A, B, returnmode = :ode)

DOCSTRING
General model for a heat exchanger which only accounts for temperature and specific heat.

OPTIONAL INPUTS
- name: DESCRIPTION
- ϵ = 0.95: DESCRIPTION
- A: DESCRIPTION
- B: DESCRIPTION
- returnmode = :ode: DESCRIPTION
"""
function Gen_HeatExchanger(; name, ϵ = 0.95, A, B, returnmode = :ode)

    eff = ϵ
    if returnmode == :ode
        ps = @parameters ϵ = ϵ

        eqs = [
            -A.Q̇ ~ ϵ * (A.p.T - B.p.T) * ((A.C < B.C) * A.C + (A.C >= B.C) * B.C)      #   heat transfer out of A -> B , if A/T > B/T 
            0 ~ A.Q̇ + B.Q̇
        ]
        ODESystem(eqs, t, [], ps; name = name, systems = [A, B], defaults = [ϵ => 0.95])
    else
        eqs = [
            -A.Q̇ ~ eff * (A.p.T - B.p.T) * ((A.C < B.C) * A.C + (A.C >= B.C) * B.C)      #   heat transfer out of A -> B , if A/T > B/T 
            0 ~ A.Q̇ + B.Q̇
        ]
        return eqs
    end
end

ModelingToolkit.@variables t

"""
    default_energy_sys()

DOCSTRING

"""
function default_energy_sys()
    sts = ModelingToolkit.@variables η_cycle(t) = 0.5 η_bop(t) = 0.5
    @named Electric   = Gas.WorkPin()                # Connects to all pumps, turbines, compressors
    @named HotUtility = Gas.HeatTransferPin()
    @named ColdUtility = Gas.HeatTransferPin()
    energy_sys = [Electric, HotUtility, ColdUtility]

    sysdict = sys2dict(energy_sys)

    return energy_sys, sts, sysdict
end


"""
    wall_circuit(; max_pressure = 80, pressrue_drop = 5, Tmin = 450 + 273.15, Tmax = 550 + 273.15, load = 1.0e8)

DOCSTRING
Returns a simple fluid circuit for a helium heat extraction circuit. Currently only allows for helium as the fluid.

OPTIONAL INPUTS
- `max_pressure = 80`: DESCRIPTION
- `pressrue_drop = 5`: DESCRIPTION
- `Tmin = 450 + 273.15`: DESCRIPTION
- `Tmax = 550 + 273.15`: DESCRIPTION
- `load = 1.0e8`: DESCRIPTION
"""
function wall_circuit(; max_pressure = 80, pressrue_drop = 5, Tmin = 450 +273.15, Tmax = 550+273.15, load = 100e6)
    params = @parameters Qwall = load
    pressure_max_wall   = max_pressure;  # bar
    pressure_drop_wall  = pressrue_drop;
    pressure_min_wall   = pressure_max_wall-pressure_drop_wall;
    Tmin_wall = Tmin;
    Tmax_wall = Tmax;
    
    @named wall_supply          = Gas.SinglePortReservoir(P = pressure_min_wall, T = Tmin_wall);
    @named wall_circulator      = Gas.PassiveThermoCompressor(η = 0.9);
    @named wall_const_pressure  = Gas.SetPressure(P=pressure_max_wall);
    @named wall_heat            = Gas.FlowControlThermoHeatTransfer(ΔP = pressure_drop_wall,Tout = Tmax_wall);
    @named wall_hx              = Gas.ThermoHeatTransfer() ; 
    @named wall_relief          = Gas.ReliefElement();
    
    wall_sys = [wall_supply,wall_circulator,wall_const_pressure,wall_heat,wall_hx,wall_relief];

    sysdict = sys2dict(wall_sys)

    wall_connections = vcat(Gas.gas_connect(wall_supply.n,wall_circulator.p,wall_relief.n),
                            Gas.gas_connect(wall_circulator.n,wall_const_pressure.p),
                            Gas.gas_connect(wall_const_pressure.n,wall_heat.p),
                            Gas.gas_connect(wall_heat.n,wall_hx.p),
                            Gas.gas_connect(wall_hx.n,wall_relief.p),
                            wall_heat.Q̇ ~ Qwall);
    return wall_sys, wall_connections, params, sysdict
end

"""
    divertor_circuit(; max_pressure = 80, pressrue_drop = 5, Tmin = 450 + 273.15, Tmax = 550 + 273.15, load = 1.0e8)

DOCSTRING
Returns a simple fluid circuit for a helium heat extraction circuit. Currently only allows for helium as the fluid.

OPTIONAL INPUTS
- `max_pressure = 80`: DESCRIPTION
- `pressrue_drop = 5`: DESCRIPTION
- `Tmin = 450 + 273.15`: DESCRIPTION
- `Tmax = 550 + 273.15`: DESCRIPTION
- `load = 1.0e8`: DESCRIPTION
"""
function divertor_circuit(; max_pressure = 80, pressrue_drop = 5, Tmin = 450 +273.15, Tmax = 550+273.15, load = 100e6)
    params = @parameters Qdivertor = load
    pressure_max_divertor = max_pressure;  # bar
    pressure_drop_divertor = pressrue_drop;
    pressure_min_divertor = pressure_max_divertor-pressure_drop_divertor;
    Tmin_divertor = Tmin;
    Tmax_divertor = Tmax;

    @named divertor_supply          = Gas.SinglePortReservoir(P = pressure_min_divertor, T = Tmin_divertor)
    @named divertor_circulator      = Gas.PassiveThermoCompressor(η = 0.9)
    @named divertor_const_pressure  = Gas.SetPressure(P=pressure_max_divertor)
    @named divertor_heat            = Gas.FlowControlThermoHeatTransfer(ΔP = pressure_drop_divertor,Tout = Tmax_divertor)
    @named divertor_hx              = Gas.ThermoHeatTransfer()  
    @named divertor_relief          = Gas.ReliefElement()
    divertor_sys = [divertor_supply,divertor_circulator,divertor_const_pressure,divertor_heat,divertor_hx,divertor_relief]
    sysdict = sys2dict(divertor_sys)
    divertor_connections = vcat(Gas.gas_connect(divertor_supply.n,divertor_circulator.p,divertor_relief.n),
                            Gas.gas_connect(divertor_circulator.n,divertor_const_pressure.p),
                            Gas.gas_connect(divertor_const_pressure.n,divertor_heat.p),
                            Gas.gas_connect(divertor_heat.n,divertor_hx.p),
                            Gas.gas_connect(divertor_hx.n,divertor_relief.p),
                            divertor_heat.Q̇ ~ Qdivertor)

    return divertor_sys, divertor_connections, params, sysdict
end

"""
    breeder_circuit(; max_pressure = 40, pressrue_drop = 8, Tmin = 750 + 273.15, Tmax = 900 + 273.15, load = 1.0e8)

DOCSTRING
Returns a simple fluid circuit for a helium heat extraction circuit. Currently only allows for PbLi as the fluid.

OPTIONAL INPUTS
- `max_pressure = 40`: DESCRIPTION
- `pressrue_drop = 8`: DESCRIPTION
- `Tmin = 750 + 273.15`: DESCRIPTION
- `Tmax = 900 + 273.15`: DESCRIPTION
- `load = 1.0e8`: DESCRIPTION
"""
function breeder_circuit(; max_pressure = 40, pressrue_drop = 8, Tmin = 750 +273.15, Tmax = 900+273.15, load = 100e6)
    params = @parameters Qbreeder = load
    pressure_max_breeder = max_pressure;  # bar
    pressure_drop_breeder = pressrue_drop;
    pressure_min_breeder = pressure_max_breeder-pressure_drop_breeder;
    Tmin_breeder = Tmin;
    Tmax_breeder = Tmax;
    
    @named breeder_supply           = Liq.SinglePortReservoir(P = pressure_min_breeder,T = Tmin_breeder)
    @named breeder_circulator       = Liq.PassiveIncompressiblePump2(η = 0.9)
    @named breeder_const_pressure   = Liq.SetPressure2(P=pressure_max_breeder)
    @named breeder_heat       		= Liq.FlowControlIncompressibleHeatTransfer(ΔP = pressure_drop_breeder, Tout = Tmax_breeder)
    @named breeder_hx               = Liq.IncompressibleHeatTransfer()
    @named breeder_relief           = Liq.ReliefElement()
    
    breeder_sys = [breeder_supply,breeder_circulator,breeder_const_pressure,breeder_heat,breeder_hx,breeder_relief]
    sysdict = sys2dict(breeder_sys)
    breeder_connections = vcat(Liq.incompressible_connect(breeder_supply.n,breeder_circulator.p,breeder_relief.n),
                            Liq.incompressible_connect(breeder_circulator.n,breeder_const_pressure.p),
                            Liq.incompressible_connect(breeder_const_pressure.n,breeder_heat.p),
                            Liq.incompressible_connect(breeder_heat.n,breeder_hx.p),
                            Liq.incompressible_connect(breeder_hx.n,breeder_relief.p),
                            breeder_heat.Q̇ ~ Qbreeder)


    return breeder_sys, breeder_connections, params, sysdict
end

"""
    feedwater_rankine(; max_pressure = 150, mid_pressure = 25, min_pressure = 0.1, ηpump = 1.0, ηturbine = 1.0, flowrate = 50)

DOCSTRING


OPTIONAL INPUTS
- `max_pressure = 150`: DESCRIPTION
- `mid_pressure = 25`: DESCRIPTION
- `min_pressure = 0.1`: DESCRIPTION
- `ηpump = 1.0`: DESCRIPTION
- `ηturbine = 1.0`: DESCRIPTION
- `flowrate = 50`: DESCRIPTION
"""
function feedwater_rankine(;
    max_pressure = 150,
    mid_pressure = 25,
    min_pressure = 0.1,
    ηpump = 1.0,
    ηturbine = 1.0,
    flowrate = 50,
)
    params = @parameters steam_ṁ = flowrate
    steam_pmid = mid_pressure
    steam_pmin = min_pressure
    steam_pmax = max_pressure
    @named steam_supply = Steam.ContinuityReservoir()
    @named steam_hp_pump =
        Steam.AdiabaticPump(Pout = steam_pmax, setpressure = true, η = ηpump)
    @named steam_boiler = Steam.SteamHeatTransfer()
    @named steam_turbine = Steam.SIMOAdiabaticTurbine(
        setpressure = true,
        Pyin = steam_pmid,
        Pzin = steam_pmin,
        ηin = ηturbine,
    )
    @named steam_lp_pump =
        Steam.AdiabaticPump(Pout = steam_pmid, setpressure = false, η = ηpump)
    @named steam_condensor = Steam.IdealCondensor()
    @named steam_openfw = Steam.OpenFeedwaterHeater()

    steam_connections = vcat(
        Steam.hydro_connect(steam_supply.n, steam_boiler.p),
        Steam.hydro_connect(steam_boiler.n, steam_turbine.p),
        Steam.hydro_connect(steam_turbine.hp.n, steam_openfw.p1),
        Steam.hydro_connect(steam_turbine.lp.n, steam_condensor.p),
        Steam.hydro_connect(steam_condensor.n, steam_lp_pump.p),
        Steam.hydro_connect(steam_lp_pump.n, steam_openfw.p2),
        Steam.hydro_connect(steam_openfw.n, steam_hp_pump.p),
        Steam.hydro_connect(steam_hp_pump.n, steam_supply.p),
        steam_boiler.p.ṁ ~ steam_ṁ,
    )

    steam_systems = [
        steam_boiler,
        steam_turbine,
        steam_lp_pump,
        steam_condensor,
        steam_openfw,
        steam_hp_pump,
        steam_supply,
    ]
    sysdict = sys2dict(steam_systems)
    return steam_systems, steam_connections, params, sysdict
end

"""
    intermediate_loop(; Pmax = 40, Pmin = 32, Nhx = 3, Tmin = 350 + 273.15, flowrate = 100)

DOCSTRING


OPTIONAL INPUTS
- `Pmax = 40`: DESCRIPTION
- `Pmin = 32`: DESCRIPTION
- `Nhx = 3`: DESCRIPTION
- `Tmin = 350 + 273.15`: DESCRIPTION
- `flowrate = 100`: DESCRIPTION
"""
function intermediate_loop(;Pmax = 40 ,Pmin = 32, Nhx = 3, Tmin = 350 + 273.15, flowrate = 100)
    params = @parameters inter_loop_ṁ = flowrate
    pressure_max_loop = Pmax;  # bar
    pressure_drop_loop = Pmax - Pmin;
    pressure_min_loop = pressure_max_loop-pressure_drop_loop;
    Tmin_loop = Tmin

    pdrop_per = pressure_drop_loop / Nhx

    @named inter_loop_supply            = Gas.SinglePortReservoir(P = pressure_min_loop, T = Tmin_loop);
    @named inter_loop_circulator        = Gas.PassiveThermoCompressor(η = 0.9);
    @named inter_loop_const_pressure    = Gas.SetPressure(P=pressure_max_loop);
    @named inter_loop_hx1		    	= Gas.ThermoHeatTransfer(ΔP = pdrop_per) ;
    @named inter_loop_relief          	= Gas.ReliefElement();

    inter_loop_sys = [inter_loop_supply,inter_loop_relief,inter_loop_circulator,inter_loop_const_pressure,inter_loop_hx1]
    
    inter_loop_connections = vcat(Gas.gas_connect(inter_loop_supply.n, inter_loop_circulator.p,inter_loop_relief.n),
                                    Gas.gas_connect(inter_loop_circulator.n,inter_loop_const_pressure.p),
                                    Gas.gas_connect(inter_loop_const_pressure.n,inter_loop_hx1.p),
                                    inter_loop_circulator.p.ṁ ~ inter_loop_ṁ);
    for i = 2:Nhx
        push!(inter_loop_sys,Gas.ThermoHeatTransfer(name = Symbol("inter_loop_hx" * "$(i)") ,ΔP = pdrop_per));
        inter_loop_connections = vcat(inter_loop_connections, Gas.gas_connect(inter_loop_sys[end-1].n, inter_loop_sys[end].p))
    end
    inter_loop_connections = vcat(inter_loop_connections, 
                                Gas.gas_connect(inter_loop_sys[end].n, inter_loop_relief.p))

    sysdict = sys2dict(inter_loop_sys)
    return inter_loop_sys, inter_loop_connections, params, sysdict
end


"""
    rankine(; steam_pmax = 30, steam_pmin = 0.75, flowrate = 50)

DOCSTRING


OPTIONAL INPUTS
- `Pmax = 40`: DESCRIPTION
- `Pmin = 32`: DESCRIPTION
- `Nhx = 3`: DESCRIPTION
- `Tmin = 350 + 273.15`: DESCRIPTION
- `flowrate = 100`: DESCRIPTION
"""

function rankine(; steam_pmax = 30, steam_pmin = 0.75, flowrate = 50)
    #presssure units in bar
    params = @parameters steam_ṁ = flowrate
    @named steam_supply = Steam.Reservoir(P = steam_pmin)
    @named steam_valve = Steam.SteamFlowSource(ṁ = flowrate)
    @named steam_pump = Steam.AdiabaticPump(Pout = steam_pmax, setpressure = true)
    @named steam_boiler = Steam.SteamHeatTransfer()
    @named steam_turbine = Steam.AdiabaticTurbine(setpressure = true, Pout = steam_pmin)
    @named steam_condensor = Steam.ReliefElement()
    steam_sys = [
        steam_supply,
        steam_valve,
        steam_pump,
        steam_boiler,
        steam_turbine,
        steam_condensor,
    ]
    steam_connections = vcat(
        Steam.hydro_connect(steam_supply.n, steam_valve.p, steam_condensor.n),
        Steam.hydro_connect(steam_valve.n, steam_pump.p),
        Steam.hydro_connect(steam_pump.n, steam_boiler.p),
        Steam.hydro_connect(steam_boiler.n, steam_turbine.p),
        Steam.hydro_connect(steam_turbine.n, steam_condensor.p),
    )

    sysdict = sys2dict(steam_sys)
    return steam_sys, steam_connections, params, sysdict
end

function brayton_cycle(; flowrate = 50, TminCycle = 300, PminCycle = 15)
    params = @parameters cycle_ṁ = flowrate
    @named cycle_supply        = Gas.TwoPortReservoir(P = PminCycle, T = TminCycle)
    @named cycle_compressor_lp = Gas.ActiveThermoCompressor(rp = 1.7, η = 0.9)
    @named cycle_intercooler_1 = Gas.Intercooler(Tout = TminCycle)
    @named cycle_compressor_mp = Gas.ActiveThermoCompressor(rp = 1.5, η = 0.95)
    @named cycle_intercooler_2 = Gas.Intercooler(Tout = TminCycle)
    @named cycle_compressor_hp = Gas.ActiveThermoCompressor(rp = 1.5, η = 0.95)
    @named cycle_heat          = Gas.ThermoHeatTransfer()
    @named cycle_turbine       = Gas.PassiveThermoTurbine()
    @named cycle_cooler        = Gas.IdealCooler()


    connections = vcat(
        Gas.gas_connect(cycle_supply.n, cycle_compressor_lp.p),
        Gas.gas_connect(cycle_compressor_lp.n, cycle_intercooler_1.p),
        Gas.gas_connect(cycle_intercooler_1.n, cycle_compressor_mp.p),
        Gas.gas_connect(cycle_compressor_mp.n, cycle_intercooler_2.p),
        Gas.gas_connect(cycle_intercooler_2.n, cycle_compressor_hp.p),
        Gas.gas_connect(cycle_compressor_hp.n, cycle_heat.p),
        Gas.gas_connect(cycle_heat.n, cycle_turbine.p),
        Gas.gas_connect(cycle_turbine.n, cycle_cooler.p),
        Gas.gas_connect(cycle_cooler.n, cycle_supply.p),
        cycle_compressor_lp.p.ṁ ~ cycle_ṁ,
    )


    cyclesys = [
        cycle_supply,
        cycle_compressor_lp,
        cycle_intercooler_1,
        cycle_compressor_mp,
        cycle_intercooler_2,
        cycle_compressor_hp,
        cycle_heat,
        cycle_turbine,
        cycle_cooler,
    ]
    # ODESystem(connections,t; name = name, systems = systemNames)
    sysdict = sys2dict(cyclesys)
    return cyclesys, connections, params, sysdict
end




# function FeedwaterRankine2(; name, Pmin = 0.1, Pmid = 10, Pmax = 150)
#     # Control elements

        #     @named gnd = Steam.ContinuityReservoir()
        #     @named valve = Steam.SteamFlowValve()
        #     @named pumpA = Steam.AdiabaticPump(Pout = Pmid, setpressure = false)
        #     @named pumpB = Steam.AdiabaticPump(Pout = Pmax, setpressure = true)

        #     #Boiler
        #     # @named boil         = IdealBoiler(Tout = 600+273)
        #     @named boil = Steam.SteamHeatTransfer()
        #     @named turbn =
        #         Steam.SIMOAdiabaticTurbine(setpressure = true, Pyin = Pmid, Pzin = Pmin, ηin = 1.0)

        #     @named cndnsr = Steam.IdealCondensor()
        #     @named openfw = Steam.OpenFeedwaterHeater()

        #     @named WorkRes = Steam.WorkPin()
        #     @named ColdUtil = Steam.HeatTransferPin()
        #     @named HotUtil = Steam.HeatTransferPin()

        #     connections = vcat(
        #         Steam.hydro_connect(openfw.n, valve.p),
        #         Steam.hydro_connect(valve.n, pumpB.p),          # pump -> boilder
        #         Steam.hydro_connect(pumpB.n, gnd.p),
        #         Steam.hydro_connect(gnd.n, boil.p),
        #         Steam.hydro_connect(boil.n, turbn.p),             # boiler -> turbine
        #         Steam.hydro_connect(turbn.hp.n, openfw.p1),        # turbine -> openfw y
        #         Steam.hydro_connect(turbn.lp.n, cndnsr.p),        #   '------> condensor
        #         Steam.hydro_connect(cndnsr.n, pumpA.p),
        #         Steam.hydro_connect(pumpA.n, openfw.p2),
        #         work_connect(WorkRes, turbn.lp.w, turbn.hp.w, pumpA.w, pumpB.w),
        #         heat_connect(ColdUtil, cndnsr.q),
        #         heat_connect(HotUtil, boil.q),
        #     )

        #     systems =
        #         [valve, pumpA, pumpB, boil, turbn, cndnsr, openfw, gnd, WorkRes, ColdUtil, HotUtil]

        #     ODESystem(connections, t; name = name, systems = systems)
        # end

        # function FwRankine3(; name)
        #     Pmax = 150
        #     Pmid2 = 40
        #     Pmid1 = 5
        #     Pmin = 0.1
        #     Tmax = 600 + 273.15
        #     @variables t
        #     @named ElectricUtil = Steam.WorkPin()
        #     @named ElectricGen = Steam.WorkPin()
        #     @named ColdUtil = Steam.HeatTransferPin()
        #     @named HotUtil = Steam.HeatTransferPin()

        #     @named reservoir = Steam.ContinuityReservoir()
        #     @named valve = Steam.SteamFlowValve()
        #     @named openfw = Steam.OpenFeedwaterHeater()
        #     @named closedfw = Steam.ClosedFeedwaterHeater()
        #     @named mixer = Steam.MixingChamber()

        #     @named pump1 = Steam.AdiabaticPump(Pout = Pmid1, setpressure = false, η = 1.0)
        #     @named pump2 = Steam.AdiabaticPump(Pout = Pmax, setpressure = false, η = 1.0)
        #     @named pump3 =
        #         Steam.AdiabaticPump(Pout = Pmax, setpressure = true, η = 1.0, controlinlet = true)

        #     @named turbine1 = Steam.AdiabaticTurbine(setpressure = true, Pout = Pmid2)
        #     @named turbine2 =
        #         Steam.SIMOAdiabaticTurbine(setpressure = true, Pyin = Pmid1, Pzin = Pmin, ηin = 1.0)

        #     @named boiler = Steam.SteamHeatTransfer()
        #     @named reheat = Steam.SteamHeatTransfer()
        #     @named condensor = Steam.IdealCondensor()


        #     connections = vcat(
        #         Steam.hydro_connect(reservoir.n, valve.p),
        #         Steam.hydro_connect(valve.n, boiler.p),
        #         Steam.hydro_connect(boiler.n, turbine1.p),
        #         Steam.hydro_connect(turbine1.n, reheat.p, closedfw.p1),
        #         Steam.hydro_connect(reheat.n, turbine2.p),
        #         Steam.hydro_connect(turbine2.hp.n, openfw.p1),
        #         Steam.hydro_connect(turbine2.lp.n, condensor.p),
        #         Steam.hydro_connect(condensor.n, pump1.p),
        #         Steam.hydro_connect(pump1.n, openfw.p2),
        #         Steam.hydro_connect(openfw.n, pump2.p),
        #         Steam.hydro_connect(pump2.n, closedfw.p2),
        #         Steam.hydro_connect(closedfw.n1, pump3.p),
        #         Steam.hydro_connect(pump3.n, mixer.p1),
        #         Steam.hydro_connect(closedfw.n2, mixer.p2),
        #         Steam.hydro_connect(mixer.n, reservoir.p),
        #         work_connect(ElectricGen, turbine1.w, turbine2.hp.w, turbine2.lp.w),
        #         work_connect(ElectricUtil, pump1.w, pump2.w, pump3.w),
        #         heat_connect(ColdUtil, condensor.q),
        #     )


        #     systems = [
        #         reservoir,
        #         valve,
        #         openfw,
        #         closedfw,
        #         mixer,
        #         pump1,
        #         pump2,
        #         pump3,
        #         turbine1,
        #         turbine2,
        #         boiler,
        #         reheat,
        #         condensor,
        #         ElectricGen,
        #         ElectricUtil,
        #         ColdUtil,
        #     ]
        #     ODESystem(connections, t; name = name, systems = systems)
        # end

        # # function FeedwaterRankine(; name, Pmin = 0.1, Pmid = 10, Pmax = 150)
        # #     @named iores = Steam.ioReservoir(P = Pmin, fixboth = false)
        # #     @named valve = Steam.SteamFlowValve()
        # #     @named pumpB = Steam.AdiabaticPump(Pout = Pmax, setpressure = true)
        # #     @named boil = Steam.SteamHeatTransfer()
        # #     @named turbine =
        # #         Steam.SIMOAdiabaticTurbine(setpressure = true, Pyin = Pmid, Pzin = Pmin, ηin = 1.0)
        # #     @named pumpA = Steam.AdiabaticPump(Pout = 10, setpressure = true)
        # #     @named condensor = Steam.IdealCondensor()
        # #     @named openfw = Steam.OpenFeedwaterHeater()

        # #     @named WorkRes = Steam.WorkPin()
        # #     @named ColdUtil = Steam.HeatTransferPin()
        # #     @named HotUtil = Steam.HeatTransferPin()

        # #     connections = vcat(
        # #         Steam.hydro_connect(pumpB.n, valve.p),
        # #         Steam.hydro_connect(valve.n, boil.p),
        # #         Steam.hydro_connect(boil.n, turbine.p),
        # #         Steam.hydro_connect(turbine.hp.n, openfw.p1),
        # #         Steam.hydro_connect(turbine.lp.n, condensor.p),
        # #         Steam.hydro_connect(condensor.n, iores.p),
        # #         Steam.hydro_connect(iores.n, pumpA.p),
        # #         Steam.hydro_connect(pumpA.n, openfw.p2),
        # #         Steam.hydro_connect(openfw.n, pumpB.p),
        # #         work_connect(WorkRes, turbine.lp.w, turbine.hp.w, pumpA.w, pumpB.w),
        # #         heat_connect(ColdUtil, condensor.q),
        # #         heat_connect(HotUtil, boil.q),
        # #     )

        # #     systems = [
        # #         valve,
        # #         boil,
        # #         turbine,
        # #         pumpB,
        # #         iores,
        # #         pumpA,
        # #         condensor,
        # #         openfw,
        # #         WorkRes,
        # #         ColdUtil,
        # #         HotUtil,
        # #     ]
        # #     ODESystem(connections, t; name = name, systems = systems)
        # # end

        # function ComplexBraytonRegen(; name)
        #     TminCycle = 300
        #     PminCycle = 15
        #     @named WorkRes = Steam.WorkPin()
        #     @named ColdUtil = Steam.HeatTransferPin()
        #     @named HotUtil = Steam.HeatTransferPin()
        #     @named res = Gas.TwoPortReservoir(P = PminCycle, T = TminCycle)
        #     @named valve = Gas.GasFlowValve()
        #     @named comp1 = Gas.ActiveThermoCompressor(rp = 1.7, η = 0.9)
        #     @named ic1 = Gas.Intercooler(Tout = TminCycle)
        #     @named comp2 = Gas.ActiveThermoCompressor(rp = 1.5, η = 0.95)
        #     @named ic2 = Gas.Intercooler(Tout = TminCycle)
        #     @named comp3 = Gas.ActiveThermoCompressor(rp = 1.5, η = 0.95)
        #     @named regen = Gas.Regenerator()
        #     @named HeatIn = Gas.ThermoHeatTransfer()
        #     @named turbine = Gas.PassiveThermoTurbine()
        #     @named cool = Gas.IdealCooler()

        #     connections = vcat(
        #         Gas.gas_connect(res.n, valve.p),
        #         Gas.gas_connect(valve.n, comp1.p),
        #         Gas.gas_connect(comp1.n, ic1.p),
        #         Gas.gas_connect(ic1.n, comp2.p),
        #         Gas.gas_connect(comp2.n, ic2.p),
        #         Gas.gas_connect(ic2.n, comp3.p),
        #         Gas.hx_connect(regen, comp3, HeatIn, turbine, cool),
        #         Gas.gas_connect(HeatIn.n, turbine.p),
        #         Gas.gas_connect(cool.n, res.p),
        #         work_connect(WorkRes, turbine.w, comp1.w, comp2.w, comp3.w),
        #         heat_connect(ColdUtil, cool.q, ic2, ic1),
        #         heat_connect(HotUtil, HeatIn.q),
        #     )


        #     systemNames = [
        #         res,
        #         valve,
        #         comp1,
        #         ic1,
        #         comp2,
        #         ic2,
        #         comp3,
        #         regen,
        #         HeatIn,
        #         turbine,
        #         cool,
        #         WorkRes,
        #         ColdUtil,
        #         HotUtil,
        #     ]
        #     ODESystem(connections, t; name = name, systems = systemNames)
        # end

        # function Water_loop(; name, Pmin = 32, Pmax = 40)
        #     @named WorkRes = Steam.WorkPin()
        #     @named ColdUtil = Steam.HeatTransferPin()

        #     @named res = Steam.ioReservoir(P = Pmin, fixboth = true)
        #     @named valve = Steam.SteamFlowValve()
        #     @named pump = Steam.AdiabaticPump(η = 0.9, setpressure = true, Pout = Pmax)
        #     @named pset = Steam.SetPressure(P = Pmin)
        #     @named HeatIn = Steam.SteamHeatTransfer()
        #     @named HeatTx = Steam.SteamHeatTransfer()
        #     @named HeatRej = Steam.IdealCondensor()
        #     @named throttle = Steam.throttle()

        #     connections = vcat(
        #         Steam.hydro_connect(res.n, valve.p),
        #         Steam.hydro_connect(valve.n, pump.p),
        #         Steam.hydro_connect(pump.n, HeatIn.p),
        #         Steam.hydro_connect(HeatIn.n, HeatTx.p),
        #         Steam.hydro_connect(HeatTx.n, throttle.p),
        #         Steam.hydro_connect(throttle.n, HeatRej.p),
        #         Steam.hydro_connect(HeatRej.n, res.p),
        #         work_connect(WorkRes, pump.w),
        #         heat_connect(ColdUtil, HeatRej.q),
        #     )

        #     systemNames = [res, valve, pump, HeatIn, HeatTx, HeatRej, throttle, WorkRes, ColdUtil]
        #     ODESystem(connections, t; name = name, systems = systemNames)
        # end

        # function He_loop(; name, Pmin = 80, Tmin = 300, Pmax = 85)
        #     @named WorkRes = Gas.WorkPin()
        #     @named ColdUtil = Gas.HeatTransferPin()

        #     @named res = Gas.TwoPortReservoir(P = Pmin, T = Tmin)
        #     @named valve = Gas.GasFlowValve()
        #     @named circulator = Gas.PassiveThermoCompressor(η = 0.9)
        #     @named pset = Gas.SetPressure(P = Pmax)
        #     @named HeatIn = Gas.ThermoHeatTransfer()
        #     @named HeatTx = Gas.ThermoHeatTransfer()
        #     @named HeatRej = Gas.ThermoHeatTransfer()
        #     @named throttle = Gas.throttle()

        #     connections = vcat(
        #         Gas.gas_connect(res.n, valve.p),
        #         Gas.gas_connect(valve.n, circulator.p),
        #         Gas.gas_connect(circulator.n, pset.p),
        #         Gas.gas_connect(pset.n, HeatIn.p),
        #         Gas.gas_connect(HeatIn.n, HeatTx.p),
        #         Gas.gas_connect(HeatTx.n, HeatRej.p),
        #         Gas.gas_connect(HeatRej.n, throttle.p),
        #         Gas.gas_connect(throttle.n, res.p),
        #         work_connect(WorkRes, circulator.w),
        #         heat_connect(ColdUtil, HeatRej.q),
        #     )

        #     systemNames =
        #         [res, valve, circulator, pset, HeatIn, HeatTx, HeatRej, throttle, WorkRes, ColdUtil]
        #     ODESystem(connections, t; name = name, systems = systemNames)
        # end

        # function He_inter_loop(; name, Pmin = 80, Tmin = 300, Pmax = 85)
        #     @named WorkRes = Gas.WorkPin()
        #     @named ColdUtil = Gas.HeatTransferPin()
        #     @named res = Gas.TwoPortReservoir(P = Pmin, T = Tmin)
        #     @named valve = Gas.GasFlowValve()
        #     @named circulator = Gas.PassiveThermoCompressor(η = 0.9)
        #     @named pset = Gas.SetPressure(P = Pmax)
        #     @named HeatInA = Gas.ThermoHeatTransfer(ΔP = 0.1)
        #     @named HeatInB = Gas.ThermoHeatTransfer(ΔP = 0.1)
        #     @named HeatInC = Gas.ThermoHeatTransfer(ΔP = 0.0)
        #     @named HeatTx = Gas.ThermoHeatTransfer(ΔP = 0.0)
        #     @named HeatRej = Gas.ThermoHeatTransfer()
        #     @named throttle = Gas.throttle()

        #     connections = vcat(
        #         Gas.gas_connect(res.n, valve.p),
        #         Gas.gas_connect(valve.n, circulator.p),
        #         Gas.gas_connect(circulator.n, pset.p),
        #         Gas.gas_connect(pset.n, HeatInA.p),
        #         Gas.gas_connect(HeatInA.n, HeatInB.p),
        #         Gas.gas_connect(HeatInB.n, HeatInC.p),
        #         Gas.gas_connect(HeatInC.n, HeatTx.p),
        #         Gas.gas_connect(HeatTx.n, HeatRej.p),
        #         Gas.gas_connect(HeatRej.n, throttle.p),
        #         Gas.gas_connect(throttle.n, res.p),
        #         heat_connect(ColdUtil, HeatRej.q),
        #         work_connect(WorkRes, circulator.w),
        #     )

        #     systemNames = [
        #         res,
        #         valve,
        #         circulator,
        #         pset,
        #         HeatInA,
        #         HeatInB,
        #         HeatInC,
        #         HeatTx,
        #         HeatRej,
        #         throttle,
        #         WorkRes,
        #         ColdUtil,
        #     ]
        #     ODESystem(connections, t; name = name, systems = systemNames)
        # end

        # function He_quad_inter_loop(; name, Pmin = 80, Tmin = 300, Pmax = 85)
        #     @named WorkRes = Gas.WorkPin()
        #     @named ColdUtil = Gas.HeatTransferPin()
        #     @named res = Gas.TwoPortReservoir(P = Pmin, T = Tmin)
        #     @named valve = Gas.GasFlowValve()
        #     @named circulator = Gas.PassiveThermoCompressor(η = 0.9)
        #     @named pset = Gas.SetPressure(P = Pmax)
        #     @named HeatInA = Gas.ThermoHeatTransfer(ΔP = 0.1)
        #     @named HeatInB = Gas.ThermoHeatTransfer(ΔP = 0.1)
        #     @named HeatInC = Gas.ThermoHeatTransfer(ΔP = 0.0)
        #     @named boilTx = Gas.ThermoHeatTransfer(ΔP = 0.0)
        #     @named reheatTx = Gas.ThermoHeatTransfer(ΔP = 0.0)
        #     @named HeatRej = Gas.ThermoHeatTransfer()
        #     @named throttle = Gas.throttle()

        #     connections = vcat(
        #         Gas.gas_connect(res.n, valve.p),
        #         Gas.gas_connect(valve.n, circulator.p),
        #         Gas.gas_connect(circulator.n, pset.p),
        #         Gas.gas_connect(pset.n, HeatInA.p),
        #         Gas.gas_connect(HeatInA.n, HeatInB.p),
        #         Gas.gas_connect(HeatInB.n, HeatInC.p),
        #         Gas.gas_connect(HeatInC.n, boilTx.p),
        #         Gas.gas_connect(boilTx.n, reheatTx.p),
        #         Gas.gas_connect(reheatTx.n, HeatRej.p),
        #         Gas.gas_connect(HeatRej.n, throttle.p),
        #         Gas.gas_connect(throttle.n, res.p),
        #         heat_connect(ColdUtil, HeatRej.q),
        #         work_connect(WorkRes, circulator.w),
        #     )

        #     systemNames = [
        #         res,
        #         valve,
        #         circulator,
        #         pset,
        #         HeatInA,
        #         HeatInB,
        #         HeatInC,
        #         boilTx,
        #         reheatTx,
        #         HeatRej,
        #         throttle,
        #         WorkRes,
        #         ColdUtil,
        #     ]
        #     ODESystem(connections, t; name = name, systems = systemNames)
        # end

        # function breeder_loop(; name, Pmin = 32, Pmax = 40, Tmin = 600)
        #     @named WorkRes = Gas.WorkPin()
        #     @named ColdUtil = Gas.HeatTransferPin()
        #     @named res = Liq.TwoPortReservoir(P = Pmin, T = Tmin)
        #     @named valve = Liq.IncompressibleFlowValve()
        #     @named pump = Liq.PassiveIncompressiblePump()
        #     @named pset = Liq.SetPressure(P = Pmax)
        #     @named HeatIn = Liq.IncompressibleHeatTransfer()
        #     @named HeatTx = Liq.IncompressibleHeatTransfer()
        #     @named HeatRej = Liq.IdealCooler()
        #     @named throttle = Liq.throttle()


        #     connections = vcat(
        #         Liq.incompressible_connect(res.n, valve.p),
        #         Liq.incompressible_connect(valve.n, pump.p),
        #         Liq.incompressible_connect(pump.n, pset.p),
        #         Liq.incompressible_connect(pset.n, HeatIn.p),
        #         Liq.incompressible_connect(HeatIn.n, HeatTx.p),
        #         Liq.incompressible_connect(HeatTx.n, HeatRej.p),
        #         Liq.incompressible_connect(HeatRej.n, throttle.p),
        #         Liq.incompressible_connect(throttle.n, res.p),
        #         heat_connect(ColdUtil, HeatRej.q),
        #         work_connect(WorkRes, pump.w),
        #     )

        #     systemNames =
        #         [res, valve, pump, pset, HeatIn, HeatTx, HeatRej, throttle, WorkRes, ColdUtil]
        #     ODESystem(connections, t; name = name, systems = systemNames)
        # end

        # function divertor_circuit(; max_pressure = 80, pressrue_drop = 5, Tmin = 450 +273.15, Tmax = 550+273.15, load = 100e6)
        #     params = @parameters Qdivertor = load
        #     pressure_max_divertor = max_pressure;  # bar
        #     pressure_drop_divertor = pressrue_drop;
        #     pressure_min_divertor = pressure_max_divertor-pressure_drop_divertor;
        #     Tmin_divertor = Tmin;
        #     Tmax_divertor = Tmax;

        #     @named divertor_supply          = Gas.SinglePortReservoir(P = pressure_min_divertor, T = Tmin_divertor)
        #     @named divertor_circulator      = Gas.PassiveThermoCompressor(η = 0.9)
        #     @named divertor_const_pressure  = Gas.SetPressure(P=pressure_max_divertor)
        #     @named divertor_heat            = Gas.FlowControlThermoHeatTransfer(ΔP = pressure_drop_divertor,Tout = Tmax_divertor)
        #     @named divertor_hx              = Gas.ThermoHeatTransfer()  
        #     @named divertor_relief          = Gas.ReliefElement()
        #     divertor_sys = [divertor_supply,divertor_circulator,divertor_const_pressure,divertor_heat,divertor_hx,divertor_relief]
        #     sysdict = sys2dict(divertor_sys)
        #     divertor_connections = vcat(Gas.gas_connect(divertor_supply.n,divertor_circulator.p,divertor_relief.n),
        #                             Gas.gas_connect(divertor_circulator.n,divertor_const_pressure.p),
        #                             Gas.gas_connect(divertor_const_pressure.n,divertor_heat.p),
        #                             Gas.gas_connect(divertor_heat.n,divertor_hx.p),
        #                             Gas.gas_connect(divertor_hx.n,divertor_relief.p),
        #                             divertor_heat.Q̇ ~ Qdivertor)

        #     return divertor_sys, divertor_connections, params, sysdict
        # end

        # function breeder_circuit(;
        #     max_pressure = 40,
        #     pressrue_drop = 8,
        #     Tmin = 750 + 273.15,
        #     Tmax = 900 + 273.15,
        # )
        #     params = @parameters Qbreeder = 100e6
        #     pressure_max_breeder = max_pressure  # bar
        #     pressure_drop_breeder = pressrue_drop
        #     pressure_min_breeder = pressure_max_breeder - pressure_drop_breeder
        #     Tmin_breeder = Tmin
        #     Tmax_breeder = Tmax

        #     @named breeder_supply =
        #         Liq.SinglePortReservoir(P = pressure_min_breeder, T = Tmin_breeder)
        #     @named breeder_circulator = Liq.PassiveIncompressiblePump2(η = 0.9)
        #     @named breeder_const_pressure = Liq.SetPressure2(P = pressure_max_breeder)
        #     @named breeder_heat = Liq.FlowControlIncompressibleHeatTransfer(
        #         ΔP = pressure_drop_breeder,
        #         Tout = Tmax_breeder,
        #     )
        #     @named breeder_hx = Liq.IncompressibleHeatTransfer()
        #     @named breeder_relief = Liq.ReliefElement()

        #     breeder_sys = [
        #         breeder_supply,
        #         breeder_circulator,
        #         breeder_const_pressure,
        #         breeder_heat,
        #         breeder_hx,
        #         breeder_relief,
        #     ]
        #     sysdict = sys2dict(breeder_sys)
        #     breeder_connections = vcat(
        #         Liq.incompressible_connect(
        #             breeder_supply.n,
        #             breeder_circulator.p,
        #             breeder_relief.n,
        #         ),
        #         Liq.incompressible_connect(breeder_circulator.n, breeder_const_pressure.p),
        #         Liq.incompressible_connect(breeder_const_pressure.n, breeder_heat.p),
        #         Liq.incompressible_connect(breeder_heat.n, breeder_hx.p),
        #         Liq.incompressible_connect(breeder_hx.n, breeder_relief.p),
        #         breeder_heat.Q̇ ~ Qbreeder,
        #     )


        #     return breeder_sys, breeder_connections, params, sysdict
        # end

        # function feedwater_rankine(;
        #     max_pressure = 150,
        #     mid_pressure = 25,
        #     min_pressure = 0.1,
        #     flowrate = 50,
        # )
        #     params = @parameters steam_ṁ = flowrate
        #     steam_pmid = mid_pressure
        #     steam_pmin = min_pressure
        #     steam_pmax = max_pressure
        #     @named steam_supply = Steam.ContinuityReservoir()
        #     @named steam_hp_pump =
        #         Steam.AdiabaticPump(Pout = steam_pmax, setpressure = true, η = 1.0)
        #     @named steam_boiler = Steam.SteamHeatTransfer()
        #     @named steam_turbine = Steam.SIMOAdiabaticTurbine(
        #         setpressure = true,
        #         Pyin = steam_pmid,
        #         Pzin = steam_pmin,
        #         ηin = 1.0,
        #     )
        #     @named steam_lp_pump =
        #         Steam.AdiabaticPump(Pout = steam_pmid, setpressure = false, η = 1.0)
        #     @named steam_condensor = Steam.IdealCondensor()
        #     @named steam_openfw = Steam.OpenFeedwaterHeater()

        #     steam_connections = vcat(
        #         Steam.hydro_connect(steam_supply.n, steam_boiler.p),
        #         Steam.hydro_connect(steam_boiler.n, steam_turbine.p),
        #         Steam.hydro_connect(steam_turbine.hp.n, steam_openfw.p1),
        #         Steam.hydro_connect(steam_turbine.lp.n, steam_condensor.p),
        #         Steam.hydro_connect(steam_condensor.n, steam_lp_pump.p),
        #         Steam.hydro_connect(steam_lp_pump.n, steam_openfw.p2),
        #         Steam.hydro_connect(steam_openfw.n, steam_hp_pump.p),
        #         Steam.hydro_connect(steam_hp_pump.n, steam_supply.p),
        #         steam_boiler.p.ṁ ~ steam_ṁ,
        #     )

        #     steam_systems = [
        #         steam_boiler,
        #         steam_turbine,
        #         steam_lp_pump,
        #         steam_condensor,
        #         steam_openfw,
        #         steam_hp_pump,
        #         steam_supply,
        #     ]
        #     sysdict = sys2dict(steam_systems)
        #     return steam_systems, steam_connections, params, sysdict
        # end

        # function feedwater_rankine2(;
        #     max_pressure = 150,
        #     mid_pressure = 10,
        #     min_pressure = 0.1,
        #     flowrate = 50,
        # )
        #     params = @parameters steam_ṁ = flowrate
        #     steam_pmid = mid_pressure
        #     steam_pmin = min_pressure
        #     steam_pmax = max_pressure
        #     # @named steam_supply         = Steam.Reservoir(P = steam_pmin) #Steam.ContinuityReservoir()
        #     # @named steam_valve          = Steam.SteamFlowSource(ṁ = flowrate)
        #     @named steam_hp_pump =
        #         Steam.AdiabaticPump(Pout = steam_pmax, setpressure = true, η = 1.0)
        #     @named steam_boiler = Steam.SteamHeatTransfer()
        #     @named steam_turbine = Steam.SIMOAdiabaticTurbine(
        #         setpressure = true,
        #         Pyin = steam_pmid,
        #         Pzin = steam_pmin,
        #         ηin = 1.0,
        #     )
        #     @named steam_lp_pump =
        #         Steam.AdiabaticPump(Pout = steam_pmid, setpressure = false, η = 1.0)
        #     @named steam_condensor = Steam.IdealCondensor()
        #     @named steam_openfw = Steam.OpenFeedwaterHeater()

        #     steam_connections = vcat(
        #         Steam.hydro_connect(steam_openfw.n, steam_hp_pump.p),
        #         Steam.hydro_connect(steam_hp_pump.n, steam_boiler.p),
        #         Steam.hydro_connect(steam_boiler.n, steam_turbine.p),
        #         Steam.hydro_connect(steam_turbine.hp.n, steam_openfw.p1),
        #         Steam.hydro_connect(steam_turbine.lp.n, steam_condensor.p),
        #         Steam.hydro_connect(steam_condensor.n, steam_lp_pump.p),
        #         Steam.hydro_connect(steam_lp_pump.n, steam_openfw.p2),
        #         steam_hp_pump.p.ṁ ~ steam_ṁ,
        #     )

        #     steam_systems = [
        #         steam_boiler,
        #         steam_turbine,
        #         steam_lp_pump,
        #         steam_condensor,
        #         steam_openfw,
        #         steam_hp_pump,
        #     ]
        #     sysdict = sys2dict(steam_systems)
        #     return steam_systems, steam_connections, params, sysdict
        # end


        # function intermediate_loop(;
        #     Pmax = 40,
        #     Pmin = 32,
        #     Nhx = 3,
        #     Tmin = 350 + 273.15,
        #     flowrate = 50,
        # )
        #     params = @parameters inter_loop_ṁ = flowrate
        #     pressure_max_loop = 40  # bar
        #     pressure_drop_loop = 8
        #     pressure_min_loop = pressure_max_loop - pressure_drop_loop
        #     Tmin_loop = 350 + 273.15

        #     pdrop_per = pressure_drop_loop / Nhx

        #     @named inter_loop_supply = Gas.SinglePortReservoir(P = pressure_min_loop, T = Tmin_loop)
        #     @named inter_loop_circulator = Gas.PassiveThermoCompressor(η = 0.9)
        #     @named inter_loop_const_pressure = Gas.SetPressure(P = pressure_max_loop)
        #     @named inter_loop_hx1 = Gas.ThermoHeatTransfer(ΔP = pdrop_per)
        #     @named inter_loop_relief = Gas.ReliefElement()

        #     inter_loop_sys = [
        #         inter_loop_supply,
        #         inter_loop_relief,
        #         inter_loop_circulator,
        #         inter_loop_const_pressure,
        #         inter_loop_hx1,
        #     ]

        #     inter_loop_connections = vcat(
        #         Gas.gas_connect(inter_loop_supply.n, inter_loop_circulator.p, inter_loop_relief.n),
        #         Gas.gas_connect(inter_loop_circulator.n, inter_loop_const_pressure.p),
        #         Gas.gas_connect(inter_loop_const_pressure.n, inter_loop_hx1.p),
        #         inter_loop_circulator.p.ṁ ~ inter_loop_ṁ,
        #     )
        #     for i = 2:Nhx
        #         push!(
        #             inter_loop_sys,
        #             Gas.ThermoHeatTransfer(name = Symbol("inter_loop_hx" * "$(i)"), ΔP = pdrop_per),
        #         )
        #         inter_loop_connections = vcat(
        #             inter_loop_connections,
        #             Gas.gas_connect(inter_loop_sys[end-1].n, inter_loop_sys[end].p),
        #         )
        #     end
        #     inter_loop_connections = vcat(
        #         inter_loop_connections,
        #         Gas.gas_connect(inter_loop_sys[end].n, inter_loop_relief.p),
        #     )

        #     sysdict = sys2dict(inter_loop_sys)
        #     return inter_loop_sys, inter_loop_connections, params, sysdict
        # end

function brayton_regenerator(; flowrate = 50, TminCycle = 300, PminCycle = 15)
    params = @parameters cycle_ṁ = flowrate
    @named cycle_supply        = Gas.TwoPortReservoir(P = PminCycle, T = TminCycle)
    @named cycle_compressor_lp = Gas.ActiveThermoCompressor(rp = 1.7, η = 0.9)
    @named cycle_intercooler_1 = Gas.Intercooler(Tout = TminCycle)
    @named cycle_compressor_mp = Gas.ActiveThermoCompressor(rp = 1.5, η = 0.95)
    @named cycle_intercooler_2 = Gas.Intercooler(Tout = TminCycle)
    @named cycle_compressor_hp = Gas.ActiveThermoCompressor(rp = 1.5, η = 0.95)
    @named cycle_heat          = Gas.ThermoHeatTransfer()
    @named cycle_turbine       = Gas.PassiveThermoTurbine()
    @named cycle_cooler        = Gas.IdealCooler()

    @named cycle_regenerator_A = Gas.ThermoHeatTransfer()
    @named cycle_regenerator_B = Gas.ThermoHeatTransfer()
    
    @named hx4 = Gen_HeatExchanger(
        A = cycle_regenerator_A,
        B = cycle_regenerator_B,
        returnmode = :eq,
    )

    connections = vcat(
        Gas.gas_connect(cycle_supply.n, cycle_compressor_lp.p),
        Gas.gas_connect(cycle_compressor_lp.n, cycle_intercooler_1.p),
        Gas.gas_connect(cycle_intercooler_1.n, cycle_compressor_mp.p),
        Gas.gas_connect(cycle_compressor_mp.n, cycle_intercooler_2.p),
        Gas.gas_connect(cycle_intercooler_2.n, cycle_compressor_hp.p),
        Gas.gas_connect(cycle_compressor_hp.n, cycle_regenerator_A.p),
        Gas.gas_connect(cycle_regenerator_A.n, cycle_heat.p),
        Gas.gas_connect(cycle_heat.n, cycle_turbine.p),
        Gas.gas_connect(cycle_turbine.n, cycle_regenerator_B.p),
        Gas.gas_connect(cycle_regenerator_B.n, cycle_cooler.p),
        Gas.gas_connect(cycle_cooler.n, cycle_supply.p),
        cycle_compressor_lp.p.ṁ ~ cycle_ṁ,
    )

    push!(connections, hx4...)

    cyclesys = [
        cycle_supply,
        cycle_compressor_lp,
        cycle_intercooler_1,
        cycle_compressor_mp,
        cycle_intercooler_2,
        cycle_compressor_hp,
        cycle_regenerator_A,
        cycle_regenerator_B,
        cycle_heat,
        cycle_turbine,
        cycle_cooler,
    ]
    # ODESystem(connections,t; name = name, systems = systemNames)
    sysdict = sys2dict(cyclesys)
    return cyclesys, connections, params, sysdict
end

"""
    hotswap!(G, n::Int64)

DOCSTRING

# Arguments:
- `G`: DESCRIPTION
- `n::Int64`: DESCRIPTION
"""
function hotswap!(G, n::Int64)
    in_v = inneighbors(G, n)
    ou_v = outneighbors(G, n)

    # @printf "Attempting to remove node %i , %s from G\n\t" n get_prop(G, n, :name)
    @assert length(in_v) == 1 && length(ou_v) == 1 "$in_v "
    in_v = in_v[1]
    ou_v = ou_v[1]

    # @printf "Incoming node %i , %s\n\t" in_v get_prop(G, in_v, :name)
    # @printf "Outgoing node %i , %s\n\t" ou_v get_prop(G, ou_v, :name)
    # in_props = props(G,Edge(in_v,n))
    ou_props = props(G, n, ou_v)

    add_edge!(G, in_v, ou_v)

    # @printf "Added Edge : (%i , %s) -> (%i , %s) \n\t" in_v get_prop(G, in_v, :name) ou_v get_prop(
    #     G,
    #     ou_v,
    #     :name,
    # )
    set_single_edge_props!(G, Edge(in_v, ou_v), ou_props)
    # @printf "with properties\n\t\t"
    # for p in collect(keys(ou_props))
    #     @printf "key: %s ,\t value: %s\n\t\t" p ou_props[p]
    # end
    rem_vertex!(G, n)
    # @printf("\n")
end

"""
    rem_all_edge!(G,nds; ignore = [],verbose = false)
    removes all edges connected to nds that connect nds to any element besides those in ignore
"""
function rem_all_edge!(G, nds; ignore = [], verbose = false)
    remd = 0
    for na in nds
        e_na = all_neighbors(G, na)
        for ena in e_na
            if ena ∈ ignore
                continue
            end
            verbose ? println("$(ena) is not in $(ignore)") : nothing

            if has_edge(G, na, ena)
                rem_edge!(G, na, ena)
                remd += 1
            end
            if has_edge(G, ena, na)
                rem_edge!(G, ena, na)
                remd += 1
            end
        end
    end
    verbose ? println("removed $(remd) edges") : nothing
    return G
end

"""
    reverse_edge!(G,src,dst)


"""
function reverse_edge!(G, src, dst)
    eprops = props(G, src, dst)
    add_edge!(G, dst, src)
    set_props!(G, dst, src, eprops)
    rem_edge!(G, src, dst)
end

"""
    directed2undir_adjacency!(adjmat)


    Inputs:
        adjmat : adjacency matrix of a directed graph
    Returns
        returns the adjacency matrix for an undirected graph
"""
function directed2undir_adjacency!(adjmat::Union{Matrix,Array})
    idx = findall(x -> x == 1, adjmat)
    for id in idx
        adjmat[id[2], id[1]] = 1
    end
end

function directed2undir_adjacency(G::AbstractGraph)
    adjmat = Array(adjacency_matrix(G))
    idx = findall(x -> x == 1, adjmat)
    for id in idx
        adjmat[id[2], id[1]] = 1
    end
    println("Graph")
    return adjmat
end
#
# Getters and setters

"""
    node_prop(G, prop)

DOCSTRING
    return vector of props, in order of nodes (first element corresponds to vertex/node 1)

# Arguments:
- `G`: DESCRIPTION
- `prop`: DESCRIPTION
"""
function node_prop(G, prop)
    return [get_prop(G, i, prop) for i = 1:nv(G)]
end

"""
    node_propdict(G, prop)

DOCSTRING
node_propdict(G,prop)
return prop dict
# Arguments:
- `G`: DESCRIPTION
- `prop`: DESCRIPTION
"""
function node_propdict(G, prop)
    return Dict([i => get_prop(G, i, prop) for i = 1:nv(G)])
end

"""
    edge_prop(G, prop)

DOCSTRING

# Arguments:
- `G`: DESCRIPTION
- `prop`: DESCRIPTION
"""
function edge_prop(G, prop)
    return [get_prop(G, e, prop) for e in collect(edges(G))]
end

"""
    edge_propdict(G, prop; edgekey = :tuple)

DOCSTRING

# Arguments:
- `G`: DESCRIPTION
- `prop`: DESCRIPTION
OPTIONAL INPUTS
- `edgekey = :tuple`: DESCRIPTION
"""
function edge_propdict(G, prop; edgekey = :tuple)
    if (edgekey) == :tuple
        return Dict([
            (e.src, e.dst) => get_prop(G, e.src, e.dst, prop) for e in collect(edges(G))
        ])
    else
        return Dict([e => get_prop(G, e, prop) for e in collect(edges(G))])
    end
end

"""
    set_single_edge_props!(G, e::Edge, pdict::Dict{Symbol, T})

DOCSTRING

# Arguments:
- `G`: DESCRIPTION
- `e::Edge`: DESCRIPTION
- `pdict::Dict{Symbol, T}`: DESCRIPTION
"""
function set_single_edge_props!(G, e::Edge, pdict::Dict{Symbol,T}) where {T<:Any}
    for k in collect(keys(pdict))
        set_prop!(G, e, k, pdict[k])
    end
end

"""
    set_default_edge_prop!(G, propname::Symbol, propval)

DOCSTRING

# Arguments:
- `G`: DESCRIPTION
- `propname::Symbol`: DESCRIPTION
- `propval`: DESCRIPTION
"""
function set_default_edge_prop!(G, propname::Symbol, propval)
    [set_prop!(G, e, propname, propval) for e in collect(edges(G))]
end

"""
    set_default_edge_prop!(G, pdict::Dict{Symbol, T})

DOCSTRING

# Arguments:
- `G`: DESCRIPTION
- `pdict::Dict{Symbol, T}`: DESCRIPTION
"""
function set_default_edge_prop!(G, pdict::Dict{Symbol,T}) where {T<:Any}
    for k in collect(keys(pdict))
        [set_prop!(G, e, pdict[k]) for e in collect(edges(G))]
    end
end

"""
    set_default_node_prop!(G,propname::Symbol,propval)
    DOCSTRING
    applys props to whole graph
"""
function set_default_node_prop!(G, propname::Symbol, propval)
    [set_prop!(G, n, propname, propval) for n in collect(vertices(G))]
end

"""
    set_default_node_prop!(G, pdict::Dict{Symbol, T})

DOCSTRING

# Arguments:
- `G`: DESCRIPTION
- `pdict::Dict{Symbol, T}`: DESCRIPTION
"""
function set_default_node_prop!(G, pdict::Dict{Symbol,T}) where {T<:Any}
    for k in collect(keys(pdict))
        [set_prop!(G, n, k, pdict[k]) for n in collect(vertices(G))]
    end
end

"""
    set_default_node_props!(G, propname::Vector{Symbol}, propval::Vector{Any})

DOCSTRING

# Arguments:
- `G`: DESCRIPTION
- `propname::Vector{Symbol}`: DESCRIPTION
- `propval::Vector{Any}`: DESCRIPTION
"""
function set_default_node_props!(G, propname::Vector{Symbol}, propval::Vector{Any})
    length(propname) == length(propval) ? nothing :
    throw(
        AssertionError(
            "length(propname) [$(length(propname))]  != length(propval) [$(length(val))]  in set_default_node_prop!s",
        ),
    )
    for i = 1:length(propname)
        [set_prop!(G, n, propname[i], propval[i]) for n in collect(vertices(G))]
    end
end

"""
    init_node_prop!(g, propname, deft)

DOCSTRING

# Arguments:
- `g`: DESCRIPTION
- `propname`: DESCRIPTION
- `deft`: DESCRIPTION
"""
function init_node_prop!(g, propname, deft)
    for i = 1:nv(g)
        set_prop!(g, i, propname, deft)
    end
end

"""
    edge_d_func(G::AbstractGraph, vj::Int64)

DOCSTRING

# Arguments:
- `G::AbstractGraph`: DESCRIPTION
- `vj::Int64`: DESCRIPTION
"""
function edge_d_func(G::AbstractGraph, vj::Int64)
    nbs = all_neighbors(G, vj)

    d = [
        edge_d_func(
            get_prop(G, vj, :pos)[1],
            get_prop(G, vj, :pos)[2],
            get_prop(G, vj, :width),
            get_prop(G, vj, :height),
            get_prop(G, vj, :normwidth),
            get_prop(G, vj, :normheight),
            get_prop(G, vi, :pos)[1],
            get_prop(G, vi, :pos)[2],
            get_prop(G, vi, :width),
            get_prop(G, vi, :height),
            get_prop(G, vi, :normwidth),
            get_prop(G, vi, :normheight),
        ) for vi in nbs
    ]

    return d
end

"""
    edge_d_func(G::AbstractGraph, vj::Int64, testpos::T)

DOCSTRING

# Arguments:
- `G::AbstractGraph`: DESCRIPTION
- `vj::Int64`: DESCRIPTION
- `testpos::T`: DESCRIPTION
"""
function edge_d_func(G::AbstractGraph, vj::Int64, testpos::T) where {T<:Point}
    nbs = all_neighbors(G, vj)

    d = [
        edge_d_func(
            testpos[1],
            testpos[2],
            get_prop(G, vj, :width),
            get_prop(G, vj, :height),
            get_prop(G, vj, :normwidth),
            get_prop(G, vj, :normheight),
            get_prop(G, vi, :pos)[1],
            get_prop(G, vi, :pos)[2],
            get_prop(G, vi, :width),
            get_prop(G, vi, :height),
            get_prop(G, vi, :normwidth),
            get_prop(G, vi, :normheight),
        ) for vi in nbs
    ]

    return d
end

"""
    setpos!(G, lay)

DOCSTRING

# Arguments:
- `G`: DESCRIPTION
- `lay`: DESCRIPTION
"""
function setpos!(G, lay)
    for i = 1:nv(G)
        set_prop!(G, i, :pos, lay[i])
    end
end
#
"""
    occupied_grid_index(G)

DOCSTRING

# Arguments:
- `G`: DESCRIPTION
"""
function occupied_grid_index(G)
    v_idx = []
    for i = 1:nv(G)
        push!(v_idx, get_prop(G, i, :grididx))
    end
    return v_idx
end

"""
    occupied_grid_positions(G)

DOCSTRING

# Arguments:
- `G`: DESCRIPTION
"""
function occupied_grid_positions(G)
    v_idx = Point2f[]
    for i = 1:nv(G)
        push!(v_idx, get_prop(G, i, :pos))
    end
    return v_idx
end

"""
    available_grid_positions(G)

DOCSTRING

# Arguments:
- `G`: DESCRIPTION
"""
function available_grid_positions(G)
    op = occupied_grid_positions(G)
    allpos = get_prop(G, :grid)
    return setdiff(allpos, op)
end

"""
    node2rect(pos, ht, wid)

DOCSTRING

# Arguments:
- `pos`: DESCRIPTION
- `ht`: DESCRIPTION
- `wid`: DESCRIPTION
"""
function node2rect(pos, ht, wid)
    x = pos[1]
    y = pos[2] - ht #shift to bottom left
    rt = GeometryBasics.Rect(Vec(x, y), Vec(wid, ht))
    return rt
end

"""
    graphBlocks(G)

DOCSTRING

# Arguments:
- `G`: DESCRIPTION
"""
function graphBlocks(G)
    lay = occupied_grid_positions(G) .* get_prop(G, :c)
    blocks = [
        node2rect(lay[i], get_prop(G, i, :height), get_prop(G, i, :width)) for
        i = 1:nv(G)
    ]
    return lay, blocks
end

"""
    edgeTuples(G)

DOCSTRING

# Arguments:
- `G`: DESCRIPTION
"""
function edgeTuples(G)
    return [(e.src, e.dst) for e in collect(edges(G))]
end

"""
    available_grid_position_with_size(G, vj)

DOCSTRING

# Arguments:
- `G`: DESCRIPTION
- `vj`: DESCRIPTION
"""
function available_grid_position_with_size(G, vj)

    c = get_prop(G, :c)   # normalized grid spacing
    available_pos = available_grid_positions(G)

    opos = node_prop(G, :pos)        #all pos
    owid = node_prop(G, :normwidth)
    ohgt = node_prop(G, :normheight)

    overflow_wid = findall(x -> x > c, owid)    # indices where the width > cell width
    overflow_ht = findall(x -> x > c, ohgt)    # indices where the height is greater than cell height

    all_indices = unique(union(overflow_ht, overflow_wid))
    # nodes which spill into the cell below
    for vremove in all_indices
        pos = opos[vremove]
        ht = ohgt[vremove]
        wd = owid[vremove]
        x = pos[1]
        y = pos[2] - ht #shift to bottom left
        rt = GeometryBasics.Rect(Vec(x, y), Vec(wd, ht))
        overflow_nodes = findall(x -> x ∈ rt, available_pos) # indexin within apos
        available_pos = setdiff(available_pos, available_pos[overflow_nodes])
    end

    hp = get_prop(G, vj, :normheight)
    wp = get_prop(G, vj, :normwidth)
    if hp <= 1 && wp <= 1
        # fits in every grid
        return available_pos
    else
        av_pos = deepcopy(available_pos)
        for ap in available_pos
            x = ap[1]
            y = ap[2] - hp #shift to bottom left
            rt = GeometryBasics.Rect(Vec(x, y), Vec(wp, hp))

            # checking if any adjacent cells are filled 
            overflow_nodes = findfirst(x -> x ∈ rt, opos) # indexin within apos
            if !isempty(overflow_nodes)
                av_pos = setdiff(av_pos, ap)
            end
        end
    end
    return av_pos
end

"""
    update_pos!(G::AbstractGraph)

DOCSTRING

# Arguments:
- `G::AbstractGraph`: DESCRIPTION
"""
function update_pos!(G::AbstractGraph)
    for vj = 1:nv(G)
        newpos = get_prop(G, vj, :pos)
        newidx = findfirst(x -> x == newpos, get_prop(G, :grid))
        set_prop!(G, vj, :pos, newpos)
        set_prop!(G, vj, :grididx, newidx)
    end
end

"""
    update_pos!(G::AbstractGraph, vj::Int64, newpos::T)

DOCSTRING

# Arguments:
- `G::AbstractGraph`: DESCRIPTION
- `vj::Int64`: DESCRIPTION
- `newpos::T`: DESCRIPTION
"""
function update_pos!(G::AbstractGraph, vj::Int64, newpos::T) where {T<:Point}
    newidx = findfirst(x -> x == newpos, get_prop(G, :grid))
    set_prop!(G, vj, :pos, newpos)
    set_prop!(G, vj, :grididx, newidx)
end

"""
    swap_pos!(G::AbstractGraph, vj::Int64, vi::Int64)

DOCSTRING

# Arguments:
- `G::AbstractGraph`: DESCRIPTION
- `vj::Int64`: DESCRIPTION
- `vi::Int64`: DESCRIPTION
"""
function swap_pos!(G::AbstractGraph, vj::Int64, vi::Int64)
    vj_pos = get_prop(G, vj, :pos)
    vj_idx = get_prop(G, vj, :grididx)
    set_prop!(G, vj, :pos, get_prop(G, vi, :pos))
    set_prop!(G, vj, :grididx, get_prop(G, vi, :grididx))
    set_prop!(G, vi, :pos, vj_pos)
    set_prop!(G, vi, :grididx, vj_idx)
end

"""
    swap_adjacents(G::AbstractGraph, vj::Int64; verbose = false, maxiter = 10)

DOCSTRING

# Arguments:
- `G::AbstractGraph`: DESCRIPTION
- `vj::Int64`: DESCRIPTION
OPTIONAL INPUTS
- `verbose = false`: DESCRIPTION
- `maxiter = 10`: DESCRIPTION
"""
function swap_adjacents(G::AbstractGraph, vj::Int64; verbose = false, maxiter = 10)
    vj_pos = get_prop(G, vj, :pos)
    vj_edge_len = sum(edge_d_func(G, vj))
    ocpos = occupied_grid_positions(G)

    mht_dst = [manhattan_distance(vj_pos, op) for op in ocpos]
    # node indexes for the closest nodes
    check_idx = sortperm(mht_dst)

    count = 1.0

    for vi in check_idx
        if vi == vj
            continue
        end

        vi_edge_length = sum(edge_d_func(G, vi))
        og_tot_length = vj_edge_len + vi_edge_length
        new_tot_length =
            sum(edge_d_func(G, vi, get_prop(G, vj, :pos))) +
            sum(edge_d_func(G, vj, get_prop(G, vi, :pos)))

        if new_tot_length < og_tot_length
            swap_pos!(G, vi, vj)
            verbose ? println("Swapped nodes $(vi) and $(vj) at iteration $(count)") :
            nothing
            return true
        end

        if count > maxiter
            verbose ? println("no matches found") : nothing
            return false
        end
        count = count + 1
    end
end

"""
    swap_adjacents_check_size(G::AbstractGraph, vj::Int64; verbose = false, maxiter = 10)

DOCSTRING

# Arguments:
- `G::AbstractGraph`: DESCRIPTION
- `vj::Int64`: DESCRIPTION
OPTIONAL INPUTS
- `verbose = false`: DESCRIPTION
- `maxiter = 10`: DESCRIPTION
"""
function swap_adjacents_check_size(
    G::AbstractGraph,
    vj::Int64;
    verbose = false,
    maxiter = 10,
)


    vj_pos = get_prop(G, vj, :pos)
    vj_edge_len = sum(edge_d_func(G, vj))
    ocpos = occupied_grid_positions(G)
    vj_wid = get_prop(G, vj, :normwidth)
    vj_ht = get_prop(G, vj, :normheight)


    mht_dst = [manhattan_distance(vj_pos, op) for op in ocpos]
    # node indexes for the closest nodes
    check_idx = sortperm(mht_dst)
    count = 1.0

    for vi in check_idx
        if vi == vj
            continue
        end

        allowable = true

        vi_wid = get_prop(G, vi, :normwidth)
        vi_ht = get_prop(G, vi, :normheight)
        vi_pos = get_prop(G, vi, :pos)

        if vi_wid != vj_wid || vi_ht != vj_ht
            # ronfirming swap is okay
            # postitions are swapped on purpose
            rect_i = node2rect(vj_pos, vi_ht, vi_wd)
            rect_j = node2rect(vi_pos, vj_ht, vj_wid)

            tochech = setdiff(ocpos, [vj_pos, vi_pos])

            ck = findall(x -> (x ∈ rect_i) || (x ∈ rect_j), ocp)

            if !isempty(ck)
                continue
            end
        end

        vi_edge_length = sum(edge_d_func(G, vi))
        og_tot_length = vj_edge_len + vi_edge_length
        new_tot_length =
            sum(edge_d_func(G, vi, get_prop(G, vj, :pos))) +
            sum(edge_d_func(G, vj, get_prop(G, vi, :pos)))

        if new_tot_length < og_tot_length
            swap_pos!(G, vi, vj)
            verbose ? println("Swapped nodes $(vi) and $(vj) at iteration $(count)") :
            nothing
            return true
        end

        if count > maxiter
            verbose ? println("no matches found") : nothing
            return false
        end
        count = count + 1
    end
end

"""
    vlinecheck(y1::Real, h1::Real, y2::Real, h2::Real)

DOCSTRING

# Arguments:
- `y1::Real`: DESCRIPTION
- `h1::Real`: DESCRIPTION
- `y2::Real`: DESCRIPTION
- `h2::Real`: DESCRIPTION
"""
function vlinecheck(y1::Real, h1::Real, y2::Real, h2::Real)
    ylow = y1 - h1
    if y2 == y1
        return true
    elseif y2 < y1 && y2 > ylow
        #case |
        #     | |
        #       |
        return true
    elseif (y2 - h2) > ylow && (y2 - h2) < y1
        return true
    elseif y2 > y1 && (y2 - h2) < ylow
        return true
    end
    return false
end

"""
    hlinecheck(x1::Real, w1::Real, x2::Real, w2::Real)

DOCSTRING

# Arguments:
- `x1::Real`: DESCRIPTION
- `w1::Real`: DESCRIPTION
- `x2::Real`: DESCRIPTION
- `w2::Real`: DESCRIPTION
"""
function hlinecheck(x1::Real, w1::Real, x2::Real, w2::Real)
    xhi = x1 + w1
    if x1 == x2
        return true
    elseif x2 > x1 && x2 < xhi
        #case x1-----xhi
        #       x2---
        return true
    elseif (x2 + w2) > x1 && x2 < x1
        #case x1-----xhi
        #   x2---  
        return true
    elseif x2 < x1 && (x2 + w2) > xhi
        return true
    end
    return false
end

"""
    vlinecheck(srcl::Vector{Line}, dst::Line)

DOCSTRING

# Arguments:
- `srcl::Vector{Line}`: DESCRIPTION
- `dst::Line`: DESCRIPTION
"""
function vlinecheck(srcl::Vector{Line}, dst::Line)
    for sl in srcl
        sldc = decompose(Point2f, sl)
        bool = vlinecheck(
            sl[1][2],
            (sl[2][2] - sl[1][2]),
            dst[1][2],
            (dst[2][2] - dst[1][2]),
        )
        if bool == true
            return true
        end
    end
    return false
end

"""
    vLineOverlap(yvec::Vector{<:Real}, y2::Real, h2::Real)

DOCSTRING

# Arguments:
- `yvec::Vector{<:Real}`: DESCRIPTION
- `y2::Real`: DESCRIPTION
- `h2::Real`: DESCRIPTION
"""
function vLineOverlap(yvec::Vector{<:Real}, y2::Real, h2::Real)
    yinter = findall(y -> (y <= y2) && (y >= y2 - h2), yvec)
    if y2 > yvec[1] && (y2 + h2) <= yvec[end]
        yinter = [1:length(xvec)...]
    end
    return yinter
end

"""
    hLineOverlap(xvec::Vector{<:Real}, x2::Real, w2::Real)

DOCSTRING

# Arguments:
- `xvec::Vector{<:Real}`: DESCRIPTION
- `x2::Real`: DESCRIPTION
- `w2::Real`: DESCRIPTION
"""
function hLineOverlap(xvec::Vector{<:Real}, x2::Real, w2::Real)
    xinter = findall(x -> (x >= x2) && (x <= x2 + w2), xvec)
    if x2 <= xvec[1] && xvec[end] <= (x2 + w2)
        xinter = [1:length(xvec)...]
    end
    return xinter
end

"""
    horizontal_visibility_graph(G::AbstractGraph)

DOCSTRING

# Arguments:
- `G::AbstractGraph`: DESCRIPTION
"""
function horizontal_visibility_graph(G::AbstractGraph)
    gg = DiGraph()
    add_vertices!(gg, nv(G))
    npos = occupied_grid_positions(G)

    xpos = [n[1] for n in npos]
    ypos = [n[2] for n in npos]
    hts = [get_prop(G, i, :height) for i = 1:nv(G)]

    for i = 1:nv(gg)
        node_x = xpos[i]
        node_y = ypos[i]
        node_h = hts[i]

        # println("Node $(i) from $(node_y), $(node_h)")

        src_pts = collect(LinRange(node_y - node_h, node_y, 100))
        # xcriteria = findall(x -> (xpos[x] > node_x) , [1:nv(G)...])
        allowable_pts = findall(
            x ->
                (xpos[x] > node_x) &&
                    vlinecheck(node_y, node_h, ypos[x], hts[x]) == true,
            [1:nv(G)...],
        )
        # @show allowable_pts
        xpts = xpos[allowable_pts]
        idx = allowable_pts[sortperm(xpts)]

        for ix in idx
            # destination line
            # println("   node $(i) connections to $ix to $(ypos[ix]),$(hts[ix])")
            yinter = vLineOverlap(src_pts, ypos[ix], hts[ix])
            if !isempty(yinter)
                # println("   Edge added $(i) => $(ix)")
                add_edge!(gg, i, ix)
                src_pts = setdiff(src_pts, src_pts[yinter])
            end
            if isempty(src_pts)
                break
            end
        end
    end
    return gg
end

"""
    vertical_visibility_graph(G::AbstractGraph; accountForWidth = true)

DOCSTRING

# Arguments:
- `G::AbstractGraph`: DESCRIPTION
OPTIONAL INPUTS
- `accountForWidth = true`: DESCRIPTION
"""
function vertical_visibility_graph(G::AbstractGraph; accountForWidth = true)
    gg = DiGraph()
    add_vertices!(gg, nv(G))
    npos = occupied_grid_positions(G)

    xpos = [n[1] for n in npos]
    ypos = [n[2] for n in npos]
    wds = 0.1 .* ones(nv(G))
    if accountForWidth == true
        wds = [get_prop(G, i, :width) for i = 1:nv(G)]
    end

    for i = 1:nv(gg)
        node_x = xpos[i]
        node_y = ypos[i]
        node_w = wds[i]

        # println("Node $(i) from $(node_x), $(node_w)")
        src_pts = collect(LinRange(node_x, node_x + node_w, 100))
        allowable_pts = findall(
            x ->
                (ypos[x] < node_y) && (hlinecheck(node_x, node_w, xpos[x], wds[x])),
            1:nv(G),
        )

        ypts = ypos[allowable_pts]
        idx = allowable_pts[reverse(sortperm(ypts))]

        for ix in idx
            # destination line
            # println("   node $(i) connections to $ix to $(xpos[ix]),$(wds[ix])")
            xinter = hLineOverlap(src_pts, xpos[ix], wds[ix])
            if !isempty(xinter)
                add_edge!(gg, i, ix)
                src_pts = setdiff(src_pts, src_pts[xinter])
            end
            if isempty(src_pts)
                break
            end
        end
    end
    return gg
end

"""
    visibility_graph(G::AbstractGraph, dir = 1)

DOCSTRING

# Arguments:
- `G::AbstractGraph`: DESCRIPTION
- `dir = 1`: DESCRIPTION
"""
function visibility_graph(G::AbstractGraph, dir = 1)
    if dir == 1
        return horizontal_visibility_graph(G)
    end
    return vertical_visibility_graph(G)
end

"""
    shiftlayout!(lay::Vector{T})

DOCSTRING

# Arguments:
- `lay::Vector{T}`: DESCRIPTION
"""
function shiftlayout!(lay::Vector{T};) where {T<:Point}
    # # full sankey
    xlay = [xy[1] for xy in lay]
    ylay = [xy[2] for xy in lay]
    xoff = -minimum(xlay)
    yoff = -minimum(ylay)
    shft = Point2f(xoff, yoff)
    lay = lay .+ shft
    return lay
end

"""
    find_common_name(stringvec::Vector{String}; verbose = false)
    returns the string shared by all members of a vector
    for example: find_common_name(["steam_boiler","steam_turbine","steam_pump"...]) returns "steam_"
        returns "" if not found
"""
function find_common_name(stringvec::Vector{String}; verbose = false)
    length(stringvec) == 1 && verbose == true ? (println("1 el"); return stringvec[1]) : nothing

    check_str = stringvec[1]
    ref_str = stringvec[2]

    common_string = "-1"
    # looping through characters of an aribitrary string element
    # finding point at which ref_str no longer contains substring of check_str
    for i = 1:length(check_str)
        ss = SubString(check_str, 1, i)
        verbose ? println(ss) : nothing

        if contains(ref_str, ss)
            common_string = ss
            continue
        else
            verbose ? println("last match = $(common_string)") : nothing
            break
        end
    end

    if length(stringvec) > 2 && common_string != "-1"
        ref_str = stringvec[3]
        while !contains(ref_str, common_string) && length(common_string) > 1
            common_string = SubString(common_string, 1, length(common_string) - 1)
            verbose ? println("revised common string = $(common_string)") : nothing
            if length(common_string) == 1
                common_string = "-1"
                break
            end
        end
    end
    verbose ? println("Returned final string: $(common_string)") : nothing
    return common_string

end

printnames(G) = [println("$i  $(get_prop(G,i,:name))") for i = 1:nv(G)]
xpnt(x) = [xinner[1] for xinner in x]
ypnt(x) = [xinner[2] for xinner in x]
ck(x) = collect(keys(x))
ce(x) = collect(edges(x))
F32r(x, n) = Float32.(round.(x, digits = n))

"""
    create_plot_graph(GG; verbose = false, removeleafs = true, remove_from_plot, toignore, words_to_flag)

DOCSTRING

# Arguments:
- `GG`: DESCRIPTION
OPTIONAL INPUTS
- `verbose = false`: DESCRIPTION
- `removeleafs = true`: DESCRIPTION
- `remove_from_plot`: DESCRIPTION
- `toignore`: DESCRIPTION
- `words_to_flag`: DESCRIPTION
"""
function create_plot_graph(
    GG;
    verbose = false,
    removeleafs = true,
    remove_from_plot = [:Electric],
    toignore = [:steam_condensor],
    words_to_flag = ["supply" "const"])


    #   removeleafs : remove all nodes with only 1 edge connected ?
    # copy graph for safe keeping
    G = deepcopy(GG)
    ## orignal graph = GG, G = copy of GG

    # Removing unnecessary comoponents FROM "G" (that don't need to be plotted)
    #   This includes all components with 1 edge only
    idx = Int[]
    
    removeleafs && (idx = reverse(findall(x -> length(all_neighbors(G, x)) == 1, [1:nv(G)...])))                # components with 1 edge only
    
    idx = reverse(sort(vcat(idx, [G[flag, :name] for flag in remove_from_plot]...)))                           # reverse sort so you can remove them in a loop
    
    for id in idx
        verbose ? println("Removing node: $(get_prop(GG,id,:name))") : nothing
        rem_vertex!(G, id)
    end

    #hot swapping
    name_search_string = string.(node_prop(G, :name))
    removedIdx = []
    for i in eachindex(words_to_flag)
        flagword = words_to_flag[i]
        push!(removedIdx, findall(x -> occursin(flagword, x), name_search_string)...)
    end
    removedNames = name_search_string[removedIdx]
    for trm in reverse(sort(removedIdx))
        hotswap!(G, trm)
    end

    verbose ? println("Components removed: $(removedNames)") : nothing

    # removing all edges (FROM "G") from the cold utility, ignoring steam condensor
    # i.e. making it so only the steam condensor connects to the cold utility (cycle heat rejection element) 
    # toig = G[toignore, :name]    #to_ignore
    # toig = [G[:cycle_cooler,:name]] #,G[:cycle_intercooler_1,:name],G[:cycle_intercooler_2,:name]]
    torm = G[:ColdUtility, :name]        #to_remove
    rem_all_edge!(G, [torm]; ignore = [G[toignore[toig], :name] for toig = eachindex(toignore)], verbose = false)

    # (FROM "G") not needed
    # all_neighbors(G,torm)
    # verbose = false
    # torem = get_prop(gg,:utility_vector)

    # finding all cyclic instances within system graph
    cycles = simplecycles_hawick_james(G)
    if verbose
        println("EXISTING CYCLES")
        for c in cycles
            println("CYCLE $c")
            for ci in c
                println(get_prop(G, ci, :name))
            end
        end
    end

    # Creating "gcopy" 
    # STEP 1) of plotting process finding independent inter_loop_supply
    gcopy = deepcopy(G)
    erem = []   # Logging array for storing the edges from gcopy removed
    epropsave = Dict()
    # STEP 2) removing all edges that are not associated with a mass flow rate variable
    # meaning flow cyles are no longer connected since the exchenge is through a heat flow edge
    # all edges that are removed are pushed to "erem"
    begin
        for e in edges(gcopy)
            if get_prop(gcopy, e, :etype) != :flow
                verbose ?
                println(
                    "Removing edge: ($(e.src),$(get_prop(gcopy,e.src,:name))) -> ($(e.dst),$(get_prop(gcopy,e.dst,:name))) [$(get_prop(gcopy,e,:etype))]",
                ) : nothing
                eprops = props(gcopy, e)
                epropsave[e] = eprops
                rem_edge!(gcopy, e)
                push!(erem, e)
            end
        end
        for e in edges(gcopy)
            if get_prop(gcopy, e, :etype) != :flow
                verbose ?
                println(
                    "Removing edge: ($(e.src),$(get_prop(gcopy,e.src,:name))) --[$(get_prop(gcopy,e,:etype))]-->($(e.dst),$(get_prop(gcopy,e.dst,:name)))",
                ) : nothing
                eprops = props(gcopy, e)
                epropsave[e] = eprops
                rem_edge!(gcopy, e)
                push!(erem, e)
            end
        end
    end

    # Printing edge types to confirm only flow edges are left
    # if verbose 
    #     for e = edges(gcopy)
    #         println(get_prop(gcopy,e,:etype))
    #     end
    # end

    # STEP 3) Creating "G2", which is the condensed graph where each "cycle" becomes a node
    # so the graph would look like : Heat source -> cycle 1 -> cycle 2 -> heat sink
    # G2 vertices have property :node_group which contains the inner nodes of each cycle

    # creating auxilary graph with only flow elements
    flow_groups = connected_components(gcopy)           # vector of connected components, 
    G2 = MetaDiGraph()                                 # graph object
    add_vertices!(G2, length(flow_groups))               # adding a vertex for each unique cycle

    for v in vertices(G2)

        set_prop!(G2, v, :node_group, flow_groups[v])

        if length(flow_groups[v]) == 1
            set_prop!(G2, v, :name, get_prop(gcopy, flow_groups[v][1], :name))
            set_prop!(
                gcopy,
                flow_groups[v][1],
                :parent,
                get_prop(gcopy, flow_groups[v][1], :name),
            )
        else
            # also attempting assign names based on shared component substrings
            flow_group_names = [
                string(get_prop(gcopy, flow_groups[v][i], :name)) for
                i = 1:length(flow_groups[v])
            ]
            verbose ? println("Flow Group $(v):", flow_group_names) : nothing
            sysname = find_common_name(flow_group_names; verbose = false)
            sysname == "-1" ? sysname = Symbol("System" * string(v)) :
            sysname = Symbol(sysname)
            set_prop!(G2, v, :name, sysname)
        end

        for childV in flow_groups[v]
            set_prop!(gcopy, childV, :parent, get_prop(G2, v, :name))
        end
    end
    verbose = false
    # STEP 4) identifing connection scheme, Note G vs gcopy, have same nodes, not edges
    embedded_dict = Dict()      # mapping G -> G2, by means of ∀v1 ∈ G, f(v1)=>(v2 ∈ G2) if v1 ∈ v2.flow_groups ) 
    # in other words embedded_dict[steam_boiler] => steam cycle in G2
    default_node_props = Dict([:input => -1, :output => -1])

    set_default_node_prop!(gcopy, Dict([:input => -1, :output => -1]))
    set_default_node_prop!(G2, default_node_props)

    # STEP 5) identifying source and sink nodes within each cycle
    for v in vertices(gcopy)                                    # for vertex "v" in gcopy (main graph)
        name = get_prop(gcopy, v, :name)                          # name = component name (symbol)
        fgv = findfirst(x -> v ∈ x, flow_groups)                # fgv = flow group that contains v
        v_flow_connections = flow_groups[fgv]                  # other nodes in the same cycle as v
        g_in_connections = inneighbors(G, v)                  # incoming adjacent nodes to v within orignal graph
        g_out_connections = outneighbors(G, v)                 # outgoing adjecent nodes to v within orignal graph

        # difference in [incoming neigbors] and [cycle components]
        # any difference in these values means that the component -> v is from another cycle 
        in_diff = setdiff(g_in_connections, v_flow_connections)    # difference in incom
        # out_diff = setdiff(g_out_connections,v_flow_connections)  

        #if length(in_diff) != 0                 # if there are components connected to vertex v, that are not within the flow group
        for id in in_diff
            # heat input
            verbose ? println("--") : nothing
            verbose && println("Cycle  $fgv , Recieving Node $(v) , $(get_prop(gcopy,v,:name))")
            # src_v = in_diff[1]
            src_v = id
            src_flow_group_idx = findfirst(x -> src_v ∈ x, flow_groups) # cycle 2
            # println(src_v)
            add_edge!(G2, src_flow_group_idx, fgv)
            set_prop!(G2, Edge(src_flow_group_idx, fgv), :referenced_edge, Edge(src_v, v))
            # println("Cycle  $fgv connected to cycle $(src_flow_group_idx),
            if verbose
            @printf "Cycle %i has input from Cycle %i \n" fgv src_flow_group_idx

            @printf "Connection details:\n   Source:\t(Cycle: %i, Node: %i, component: %s) \n   Reciever:\t(Cycle: %i, Node: %i, component: %s)\n" src_flow_group_idx src_v get_prop(
                gcopy,
                src_v,
                :name,
            ) fgv v get_prop(gcopy, v, :name)
            end
            # (Node $(src_v), $(get_prop(gcopy,src_v,:name))) -> (Node $(v) , $(get_prop(gcopy,v,:name)))")
            # println("$src_v => $src_flow_group_idx")
            set_prop!(G2, src_flow_group_idx, :output, src_v)
            set_prop!(G2, fgv, :input, v)
            verbose ?
            println("   Node group $(fgv) Heat reciever = $(get_prop(gcopy,v,:name))") :
            nothing
            verbose ?
            println(
                "   Node group $(src_flow_group_idx) Heat rejector = $(get_prop(gcopy,src_v,:name))",
            ) : nothing
            verbose ?
            println(
                "   G has edge  $(get_prop(gcopy,src_v,:name)) => $(get_prop(gcopy,v,:name)) ? $(has_edge(gcopy,src_v,v))",
            ) : nothing
        end
        embedded_dict[v] = fgv
        # println("Node $(v) in group $(fgv)")
    end
    # end

    cycles = simplecycles_hawick_james(gcopy)
    xs, ys, paths = solve_positions(Zarate(), G2)
    set_default_edge_prop!(G2, :weight, 1)
    # quick_plot(G2,xs,ys,paths)
    # quick_plot(gcopy,xs,ys,paths)
    # begin
    # handeling cycles with more than 1 input or output by removing interconnecting edges (so they can be in the same layer)
    #     input_cycles = [(v == -1 ? -1 : embedded_dict[v]) for v in node_prop(G2,:input)]
    #     output_cycles =  [(v == -1 ? -1 : embedded_dict[v]) for v in node_prop(G2,:output)]
    #     ##
    #     G2_edges = collect(edges(G2))
    #     G2_src = [e.src for e in G2_edges]
    #     G2_dst = [e.dst for e in G2_edges]

    #     # input = more than one oncoming
    #     # output = more than one outgoin
    #     multi_output_cycles  = [G2_src[z] for z in findall(x -> length(findall(y -> y == x, G2_src)) > 1, unique(G2_src))]
    #     multi_input_cycles   = [G2_dst[z] for z in findall(x -> length(findall(y -> y == x, G2_dst)) > 1, unique(G2_dst))]

    #     get_prop(G2, multi_input_cycles[1], :name)
    #     get_prop(G2, multi_output_cycles[1], :name)
    # #

    # this is temp solution, it would be better to remove the edges and add later so they can be in same layer
    for n = 1:nv(G2)
        inn = inneighbors(G2, n)
        if length(inn) > 2

            incoming_edges = [Edge(n_in, n) for n_in in inn]
            recieving_nodes =
                [get_prop(G2, e, :referenced_edge).dst for e in incoming_edges] # edges in gcopy
            # println(recieving_nodes)

            # finding, for each recieving node, the number of neighbors which are also recieving nodes
            f = [
                length(intersect(all_neighbors(gcopy, rn), recieving_nodes)) for
                rn in recieving_nodes
            ]
            set_prop!(G2, n, :input, recieving_nodes[findfirst(x -> x == maximum(f), f)])
        end
        otn = outneighbors(G2, n)
        if length(otn) > 2
            outgoing_edges = [Edge(n, n_out) for n_out in otn]
            outgoing_nodes = [get_prop(G2, e, :referenced_edge).src for e in outgoing_edges] # edges in gcopy
            f = [
                length(intersect(all_neighbors(gcopy, ogn), outgoing_nodes)) for
                ogn in outgoing_nodes
            ]
            set_prop!(G2, n, :input, outgoing_nodes[findfirst(x -> x == maximum(f), f)])
        end
    end

    ## 
    for v in vertices(G2)
        # fgv = findfirst(x -> v ∈ x, flow_groups)
        inn = get_prop(G2, v, :input)
        verbose && println(inn)
        outt = get_prop(G2, v, :output)
        if inn != -1
            verbose && println("   Node group $(v) Heat reciever = $(get_prop(G,inn,:name))")
        end
        if outt != -1
            verbose && println("                Heat Rejector = $(get_prop(G,outt,:name))")
        end
        # println("   Node group $(src_flow_group_idx) Heat rejector = $(get_prop(G2,src_v,:name))")
    end

    # changing edge directions # IMPORTANT
    set_default_edge_prop!(gcopy, :isReversed, false)
    for v in vertices(G2)
        verbose ?
        println(
            "NODE $v Input = $(get_prop(G2,v,:input)), Output = $(get_prop(G2,v,:output))",
        ) : nothing
        inn = get_prop(G2, v, :input)
        outt = get_prop(G2, v, :output)
        count = 0
        while inn != -1 && outt != -1 && has_path(gcopy, outt, inn) && (inn != outt)
            # while there exists a path from outlet to inlet
            intr_nodes = enumerate_paths(dijkstra_shortest_paths(gcopy, outt), inn)
            for i = 2:length(intr_nodes)
                ndst = intr_nodes[i]
                nsrc = intr_nodes[i-1]

                eprops = props(gcopy, nsrc, ndst)
                eprops[:isReversed] = !eprops[:isReversed]

                rem_edge!(gcopy, nsrc, ndst)
                add_edge!(gcopy, ndst, nsrc)

                set_single_edge_props!(gcopy, Edge(ndst, nsrc), eprops)
            end
            count = count + 1
            if count > 5000
                break
            end
        end
        verbose && println("Took $count itterations to reverse route path $(outt) => $(inn)")
    end
    ##
    non_flow_edges = Pair{Int64,Int64}[]
    for e in erem
        add_edge!(gcopy, e.src, e.dst)
        set_props!(gcopy, e.src, e.dst, epropsave[e])
        push!(non_flow_edges, e.src => e.dst)
        set_prop!(gcopy, e, :isReversed, false)
        # println("($(e.src),$(get_prop(gcopy,e.src,:name))) -> ($(e.dst),$(get_prop(gcopy,e.dst,:name))) [$(get_prop(gcopy,e,:etype))]")
    end

    cycles = simplecycles_hawick_james(gcopy)
    # println(cycles)

    if verbose
        for e in edges(gcopy)
            println("$(get_prop(G,e.src,:name)) => $(get_prop(G,e.dst,:name))")
        end

        for n = 1:nv(gcopy)
            println("$n = $(get_prop(gcopy,n,:name))")
        end
    end
    defaultweight!(gcopy, 1)
    weightfield(gcopy)
    set_default_edge_prop!(gcopy, :weight, 1)
    set_prop!(gcopy, :reduced_graph, G2)
    util_labels = string.(get_prop(gcopy, :utility_vector))
    system_labels = string.(node_prop(G2, :name))
    set_prop!(gcopy, :system_labels, system_labels)
    set_prop!(gcopy, :util_labels, util_labels)
    return gcopy
    # edge_prop(gcopy,weightfield(gcopy))
    # weights_mat = Dict([e.src,e.dst] => get_prop(gcopy,e,:weight) for e in collect(edges(gcopy)))
end


"""
    layers_to_force!(gcopy; maxAlignCnt = 3, doplot = false, plotattr)

DOCSTRING
    layers_to_force!(
        gcopy;
        maxAlignCnt = 3,
        doplot = false,
        plotattr = (
            linewidth = 0.5,
            marker = :circle,
            markersize = 1.5,
            markercolor = :red,
            color = :black,
            alpha = 0.7,
            legend = false,
        ),
    )
    takes in graph, returns 
    maxAlignCnt specifies the largest groupings to requires equal layering, sometimes zarate struggles with more than 3
    requirementsf ror force_equal_layers = xLayReqs,
                force_order = vSortReqs,
    from recent output xs,ys,paths
# Arguments:
- `gcopy: graph
OPTIONAL INPUTS
- `maxAlignCnt = 3`: DESCRIPTION
- `doplot = false`: DESCRIPTION
- `plotattr`: DESCRIPTION
"""
function layers_to_force!(
    gcopy;
    maxAlignCnt = 3,
    doplot = false,
    plotattr = (
        linewidth = 0.5,
        marker = :circle,
        markersize = 1.5,
        markercolor = :red,
        color = :black,
        alpha = 0.7,
        legend = false,
    ),
)
    xs, ys, paths = solve_positions(Zarate(), gcopy)#force_layer=[3=>14], force_order =[10=>6, 7=>3, 3=>2],
    # doplot ? quick_plot(gcopy, xs, ys, paths) : nothing
    path_pts = collect(values(paths))

    # node positions
    xs = F32r(xs, 3)
    ys = F32r(ys, 3)

    # path positions
    xpath = F32r.(xpnt(path_pts), 3)
    ypath = F32r.(ypnt(path_pts), 3)

    # Indexes with more than 2 pts (center point), applys to both x and p paths
    inter_pt_inx = findall(x -> length(x) > 2, xpath)

    # Vector of unique points (counted correctly)
    xLayActuals = xs
    xLayInters = vcat([x[2:end-1] for x in xpath[inter_pt_inx]]...)
    # yLayInters  = vcat([y[2:end-1] for y in ypath[inter_pt_inx]]...)
    xLayPts = vcat(xLayActuals, xLayInters)
    #y LayPts = vcat(ys,vcat(ypath[inter_pt_inx]...))
    # unique x points
    xUnqLay = sort(unique(xLayPts))

    # count of how many points are positioned at the different unique x values
    xLayRealCnt = [count(==(xp), xLayActuals, dims = 1) for xp in xUnqLay]
    xLayInterCnt = [count(==(xp), xLayInters, dims = 1) for xp in xUnqLay]

    xLayCnt = [count(==(xp), xLayPts, dims = 1) for xp in xUnqLay]
    # xlaydict    = [i => xLayCnt[i] for i in eachindex(xLayCnt)]

    xLayReqs = Pair{Int64,Int64}[]
    vSortReqs = Pair{Int64,Int64}[]    # alphabetical
    # xLayCnt = [count(==(xp),xpath,dims=1)-count(==(xp),xs,dims=1) for xp in xUnqLay]    
    nprop_parent_lookup = node_propdict(gcopy, :parent)
    npl = nprop_parent_lookup
    lay2node = Dict()
    for xul in eachindex(xUnqLay)                       # xul = index from layer 1 to end
        laynodes = findall(x -> x == xul, xs)            #nodes in that layer
        lay2node[xul] = laynodes
        if length(laynodes) > 1
            laynodeparents = [npl[ln] for ln in laynodes]
            laynodes1 = laynodes[sortperm(string.(laynodeparents))]
            laynodeparents1 = sort(laynodeparents)
            for i = 1:length(laynodes1)-1
                if laynodeparents[i] != laynodeparents[i+1]
                    push!(vSortReqs, (laynodes1[i] => laynodes1[i+1]))
                end
            end
        end
        for ln in laynodes
            set_prop!(gcopy, ln, :layer, xul)
        end
    end

    # unique counts 
    xUnqCnt = reverse(sort(unique(xLayCnt)))
    rule2cnt = Dict()

    for xuc in xUnqCnt
        if xuc[1] != 1
            # all layers with the same amount of points

            shared_layers = findall(x -> x == xuc, xLayCnt)
            real_ct = xLayRealCnt[shared_layers]
            fake_ct = xLayInterCnt[shared_layers]
            rating = fake_ct - real_ct    # shared_layers = shared_layers[findall(x -> x[1] < 2 ,rating)]
            acceptable_idx = findall(x -> x[1] < 2, fake_ct)
            applyto = shared_layers[acceptable_idx]
            if length(applyto) > 1 && xuc[1] < maxAlignCnt
                ref_lay = applyto[1]
                for forced_lay in applyto[2:end]
                    # @show typeof(ref_lay => forced_lay)
                    push!(xLayReqs, (ref_lay => forced_lay))
                    rule2cnt[(ref_lay=>forced_lay)] = xuc[1]
                    ref_lay = forced_lay
                end
            end
        end
    end
    xLayReqs
    # shuffle!(xLayReqs)# didnt do anyting
    xLayApplied = Pair{Int64,Int64}[]
    xs, ys, paths = solve_positions(
        Zarate(),
        gcopy;
        force_equal_layers = xLayReqs,
        force_order = vSortReqs,
    )
    # xs, ys, paths = solve_positions(Zarate(), gcopy;  force_order = vSortReqs)
    # if doplot
    #     quick_plot(gcopy, xs, ys, paths; nametext = true, plotattr = plotattr)
    # end
    return xLayReqs, vSortReqs, xs, ys, paths, lay2node
end

"""
    initialize_plot_props!(gcopy, lay2node, xs, ys, paths)

DOCSTRING

# Arguments:
- `gcopy`: DESCRIPTION
- `lay2node`: DESCRIPTION
- `xs`: DESCRIPTION
- `ys`: DESCRIPTION
- `paths`: DESCRIPTION
"""
function initialize_plot_props!(gcopy, lay2node,xs,ys,paths)
    default_plot_properties =
        Dict([:normheight => 1, :normwidth => 1, :height => 1, :width => 1])
    G2 = get_prop(gcopy, :reduced_graph)
    set_default_node_prop!(gcopy, default_plot_properties)
    lay = [Point2f(xs[i], ys[i]) for i = 1:nv(gcopy)]
    setpos!(gcopy, lay)
    maxh = 5
    for k in collect(keys(lay2node))
        nodes = lay2node[k]
        numnodes = length(nodes)
        hts = maxh / numnodes
        for n in nodes
            set_prop!(gcopy, n, :height, hts)
            set_prop!(gcopy, n, :normheight, hts / maxh)
        end
    end

    set_prop!(gcopy, :paths, paths)
    [set_prop!(gcopy, e, :path, paths[e]) for e in collect(edges(gcopy))]


    util_labels = string.(get_prop(gcopy, :utility_vector))
    system_labels = string.(node_prop(G2, :name))
    component_full_name = string.(node_prop(gcopy, :name))
    component_display_name = deepcopy(component_full_name)

    set_prop!(gcopy, :system_labels, system_labels)
    set_prop!(gcopy, :util_labels, util_labels)

    for i in eachindex(component_display_name)
        cfn = component_display_name[i]
        groupidx = findfirst(x -> occursin(x, cfn), system_labels)
        if !isnothing(groupidx) && !(cfn ∈ util_labels)
            component_display_name[i] = replace(cfn, (system_labels[groupidx] .=> ""))
        end
        set_prop!(gcopy, i, :displayName, component_display_name[i])
    end

    rev_status = edge_prop(gcopy, :isReversed)
    flipped_edges = findall(x -> x == false, rev_status)
    edge_ids = [[c.src, c.dst] for c in ce(gcopy)]
    # edge_ids[flipped_edges] .= reverse!.(edge_ids[flipped_edges])

    # this is dope, I should do this all the time, create a checker function 
    checker(x) = x == 1 ? 1 : -1
    directional_value = checker.(rev_status)
    # quick_plotG(gcopy)
end


"""
    add_plot_elments(gcopy; verbose = false, default_plot_properties = Dict([:displayName => "", :normheight => 1, :normwidth => 1, :height => 1, :width => 1, :nodeType => :fake]))

DOCSTRING

# Arguments:
- `gcopy`: DESCRIPTION
OPTIONAL INPUTS
- `verbose = false`: DESCRIPTION
- `default_plot_properties = Dict([:displayName => "", :normheight => 1, :normwidth => 1, :height => 1, :width => 1, :nodeType => :fake])`: DESCRIPTION
"""
function add_plot_elments(
    gcopy;
    verbose = false,
    default_plot_properties = Dict([
        :displayName => "",
        :normheight => 1,
        :normwidth => 1,
        :height => 1,
        :width => 1,
        :nodeType => :fake,
    ]),
)
    gc = deepcopy(gcopy)
    xs, ys, paths = solve_positions(Zarate(), gc)#force_layer=[3=>14], force_order =[10=>6, 7=>3, 3=>2],
    # quick_plot(gc,xs,ys,paths)
    path_pts = collect(values(paths))
    # node positions
    xs = F32r(xs, 3)
    ys = F32r(ys, 3)
    lay = [Point2f(xs[i], ys[i]) for i = 1:nv(gc)]
    setpos!(gc, lay)
    set_default_edge_prop!(gc, :sectionId, 1)
    set_default_edge_prop!(gc, :nSections, 1)
    set_default_node_prop!(gc, :nodeType, :real)
    default_plot_properties = Dict()
    default_plot_properties = Dict([
        :displayName => "",
        :normheight => 1,
        :normwidth => 1,
        :height => 1,
        :width => 1,
        :nodeType => :fake,
    ])
    countt = 1
    for e in collect(edges(gc))
        eprops = props(gc, e)
        epathx, epathy = paths[e]
        refnode_props = props(gc, e.src)
        refnode_keys = collect(keys(refnode_props))
        if length(epathx) > 2 #removing edge and replacing with fake nodes
            nsections = length(epathx) - 1
            nnodes = nsections - 1
            startnode = e.src
            for i = 1:nnodes
                add_vertex!(gc)
                set_props!(
                    gc,
                    nv(gc),
                    Dict([
                        :normheight => 1,
                        :normwidth => 1,
                        :height => 1,
                        :width => 1,
                        :nodeType => :fake,
                        :name => Symbol("plotnode" * "$countt"),
                        :displayName => "",
                    ]),
                )
                set_prop!(gc, nv(gc), :pos, Point2f(epathx[i+1], epathy[i+1]))
                for rpk in refnode_keys
                    if has_prop(gc, nv(gc), rpk) == false
                        verbose ? println("$rpk") : nothing
                        set_prop!(gc, nv(gc), rpk, refnode_props[rpk])
                    end
                end
                add_edge!(gc, startnode, nv(gc))
                set_single_edge_props!(gc, Edge(startnode, nv(gc)), eprops)
                set_single_edge_props!(
                    gc,
                    Edge(startnode, nv(gc)),
                    Dict([:sectionId => i, :nSections => nsections]),
                )
                startnode = nv(gc)
                countt = countt + 1
            end
            add_edge!(gc, nv(gc), e.dst)
            set_single_edge_props!(gc, Edge(nv(gc), e.dst), eprops)
            set_single_edge_props!(
                gc,
                Edge(nv(gc), e.dst),
                Dict([:sectionId => nsections, :nSections => nsections]),
            )
            rem_edge!(gc, e.src, e.dst)
        end
    end
    return gc
    #
        # # path positions
        # xpath = F32r.(xpnt(path_pts),3)
        # ypath = F32r.(ypnt(path_pts),3)

        # # Indexes with more than 2 pts (center point), applys to both x and p paths
        # inter_pt_inx = findall(x -> length(x) > 2, xpath)

        # # Vector of unique points (counted correctly)
        # xLayActuals = xs
        # xLayInters  = vcat([x[2:end-1] for x in xpath[inter_pt_inx]]...)
        # yLayInters  = vcat([y[2:end-1] for y in ypath[inter_pt_inx]]...)
        # xLayPts     = vcat(xLayActuals,xLayInters)
        # #y LayPts = vcat(ys,vcat(ypath[inter_pt_inx]...))
        # # unique x points
        # xUnqLay = sort(unique(xLayPts)) 

        # # count of how many points are positioned at the different unique x values
        # xLayRealCnt     = [count(==(xp),xLayActuals,dims=1) for xp in xUnqLay]     
        # xLayInterCnt    = [count(==(xp),xLayInters,dims=1) for xp in xUnqLay]   

        # xLayCnt     = [count(==(xp),xLayPts,dims=1) for xp in xUnqLay]      
        # xlaydict    = [i => xLayCnt[i] for i in eachindex(xLayCnt)]
end

"""
    add_plot_elments!(gcopy; verbose = false, default_plot_properties = Dict([:displayName => "", :normheight => 1, :normwidth => 1, :height => 1, :width => 1, :nodeType => :fake]))

DOCSTRING

# Arguments:
- `gcopy`: DESCRIPTION
OPTIONAL INPUTS
- `verbose = false`: DESCRIPTION
- `default_plot_properties = Dict([:displayName => "", :normheight => 1, :normwidth => 1, :height => 1, :width => 1, :nodeType => :fake])`: DESCRIPTION
"""
function add_plot_elments!(
    gcopy;
    verbose = false,
    default_plot_properties = Dict([
        :displayName => "",
        :normheight => 1,
        :normwidth => 1,
        :height => 1,
        :width => 1,
        :nodeType => :fake,
    ]),
)
    xs, ys, paths = solve_positions(Zarate(), gcopy)#force_layer=[3=>14], force_order =[10=>6, 7=>3, 3=>2],
    # quick_plot(gcopy,xs,ys,paths)
    path_pts = collect(values(paths))
    # node positions
    xs = F32r(xs, 3)
    ys = F32r(ys, 3)
    lay = [Point2f(xs[i], ys[i]) for i = 1:nv(gcopy)]
    setpos!(gcopy, lay)
    set_default_edge_prop!(gcopy, :sectionId, 1)
    set_default_edge_prop!(gcopy, :nSections, 1)
    set_default_node_prop!(gcopy, :nodeType, :real)
    default_plot_properties = Dict()
    default_plot_properties = Dict([
        :displayName => "",
        :normheight => 1,
        :normwidth => 1,
        :height => 1,
        :width => 1,
        :nodeType => :fake,
    ])
    countt = 1
    for e in collect(edges(gcopy))
        eprops = props(gcopy, e)
        epathx, epathy = paths[e]
        refnode_props = props(gcopy, e.src)
        refnode_keys = collect(keys(refnode_props))
        if length(epathx) > 2 #removing edge and replacing with fake nodes
            nsections = length(epathx) - 1
            nnodes = nsections - 1
            startnode = e.src
            for i = 1:nnodes
                add_vertex!(gcopy)
                set_props!(
                    gcopy,
                    nv(gcopy),
                    Dict([
                        :normheight => 1,
                        :normwidth => 1,
                        :height => 1,
                        :width => 1,
                        :nodeType => :fake,
                        :name => Symbol("plotnode" * "$countt"),
                        :displayName => "",
                    ]),
                )
                set_prop!(gcopy, nv(gcopy), :pos, Point2f(epathx[i+1], epathy[i+1]))
                for rpk in refnode_keys
                    if has_prop(gcopy, nv(gcopy), rpk) == false
                        verbose ? println("$rpk") : nothing
                        set_prop!(gcopy, nv(gcopy), rpk, refnode_props[rpk])
                    end
                end
                add_edge!(gcopy, startnode, nv(gcopy))
                set_single_edge_props!(gcopy, Edge(startnode, nv(gcopy)), eprops)
                set_single_edge_props!(
                    gcopy,
                    Edge(startnode, nv(gcopy)),
                    Dict([:sectionId => i, :nSections => nsections]),
                )
                startnode = nv(gcopy)
                countt = countt + 1
            end
            add_edge!(gcopy, nv(gcopy), e.dst)
            set_single_edge_props!(gcopy, Edge(nv(gcopy), e.dst), eprops)
            set_single_edge_props!(
                gcopy,
                Edge(nv(gcopy), e.dst),
                Dict([:sectionId => nsections, :nSections => nsections]),
            )
            rem_edge!(gcopy, e.src, e.dst)
        end
    end
    #
        # # path positions
        # xpath = F32r.(xpnt(path_pts),3)
        # ypath = F32r.(ypnt(path_pts),3)

        # # Indexes with more than 2 pts (center point), applys to both x and p paths
        # inter_pt_inx = findall(x -> length(x) > 2, xpath)

        # # Vector of unique points (counted correctly)
        # xLayActuals = xs
        # xLayInters  = vcat([x[2:end-1] for x in xpath[inter_pt_inx]]...)
        # yLayInters  = vcat([y[2:end-1] for y in ypath[inter_pt_inx]]...)
        # xLayPts     = vcat(xLayActuals,xLayInters)
        # #y LayPts = vcat(ys,vcat(ypath[inter_pt_inx]...))
        # # unique x points
        # xUnqLay = sort(unique(xLayPts)) 

        # # count of how many points are positioned at the different unique x values
        # xLayRealCnt     = [count(==(xp),xLayActuals,dims=1) for xp in xUnqLay]     
        # xLayInterCnt    = [count(==(xp),xLayInters,dims=1) for xp in xUnqLay]   

        # xLayCnt     = [count(==(xp),xLayPts,dims=1) for xp in xUnqLay]      
        # xlaydict    = [i => xLayCnt[i] for i in eachindex(xLayCnt)]
end


"""
    setVerticalSpacing!(gc; vspan = 5.0, pad = 2.0, doplot = true, plotattr)

DOCSTRING

# Arguments:
- `gc`: DESCRIPTION
OPTIONAL INPUTS
- `vspan = 5.0`: DESCRIPTION
- `pad = 2.0`: DESCRIPTION
- `doplot = true`: DESCRIPTION
- `plotattr`: DESCRIPTION
"""
function setVerticalSpacing!(
    gc;
    vspan = 5.0,
    pad = 2.0,
    doplot = true,
    plotattr = (
        linewidth = 1,
        marker = :circle,
        markersize = 1.5,
        markercolor = :red,
        color = :black,
        alpha = 0.7,
        legend = false,
    ),
)

    gvert = vertical_visibility_graph(gc; accountForWidth = false)
    vstacks = connected_components(gvert)
    npos = node_prop(gc, :pos)
    nposy = ypnt(npos)
    for i in eachindex(vstacks)
        compgroup = vstacks[i]
        nInLayer = length(compgroup)
        layerorder = sortperm(nposy[compgroup])
        yspacing = vspan / (nInLayer)
        new_ypts = -(vspan / 2.0) + yspacing / 2.0
        for idx = 1:nInLayer
            layerorder[idx]                     # lowest node first
            comp = compgroup[layerorder[idx]]
            oldpos = npos[comp]
            set_prop!(gc, comp, :pos, Point2f([oldpos[1], new_ypts]))
            new_ypts += yspacing
        end
    end
    npos = node_prop(gc, :pos)
    xs = [nn[1] for nn in npos]
    ys = [nn[2] for nn in npos]
    # set_prop!(gc,:vertical_visibility_graph,gvert)
    # quick_plot(gc, xs, ys, paths; nametext = true, mode = :patha, plotattr = plotattr)
    return xs, ys
end

"""
    setLayerWidth!(gc; pad = 0.5, verbose = false)

DOCSTRING

# Arguments:
- `gc`: DESCRIPTION
OPTIONAL INPUTS
- `pad = 0.5`: DESCRIPTION
- `verbose = false`: DESCRIPTION
"""
function setLayerWidth!(gc; pad = 0.5, verbose = false)
    gvert = vertical_visibility_graph(gc)
    vstacks = connected_components(gvert)
    nwids = node_prop(gc, :width)
    maxwids = [maximum(nwids[compgroup]) for compgroup in vstacks]


    npos = node_prop(gc, :pos)
    y = [n[2] for n in npos]
    x = [n[1] for n in npos]
    xcurrent = [x[vx[1]] for vx in vstacks]
    layorder = sortperm(xcurrent)

    ordered_xpos = xcurrent[layorder]
    ordered_maxwids = maxwids[layorder]
    ordered_maxwids_with_pad = ordered_maxwids .+ vcat(0, pad .* ones(length(maxwids) - 1))
    xshifted = cumsum(ordered_maxwids_with_pad)
    verbose ? println("Position  : ", F32r(ordered_xpos, 2)) : nothing
    verbose ? println("Widths    : ", F32r(ordered_maxwids, 2)) : nothing
    verbose ? println("Wdth+Pad  : ", F32r(ordered_maxwids_with_pad, 2)) : nothing
    verbose ? println("New x pos : ", F32r(xshifted, 2)) : nothing
    node_propdict(gc, :layer)

    for i in eachindex(layorder)
        lay_i = layorder[i]         # index of i'th closest layer
        verbose ?
        println(
            "Layer i: ",
            i,
            ", index",
            lay_i,
            ",  x = $(xcurrent[lay_i]), x new = $(xshifted[i])",
        ) : nothing
        if i == 1
            continue
        end
        cgroup = vstacks[lay_i]
        # println("LAYER ",lay_i,"  compontnts = $(cgroup)")
        for cg in cgroup
            set_prop!(gc, cg, :pos, Point2f([xshifted[i], y[cg]]))
        end
    end
    # plotplant(gc)
end

"""
    edgeroute(p1::T, p2::T)

DOCSTRING

# Arguments:
- `p1::T`: DESCRIPTION
- `p2::T`: DESCRIPTION
"""
function edgeroute(p1::T, p2::T) where {T<:Point}
    x1 = p1[1]
    x2 = p2[1]
    y1 = p1[2]
    y2 = p2[2]

    dx = abs(x2 - x1)
    dy = abs(y2 - y1)
    if dx == 0 || dy == 0
        return [x1, x2], [y1, y2]
    end
    if dx > dy
        return [x1, x2, x2], [y1, y1, y2]
    else #dy > dx
        return [x1, x1, x2], [y1, y2, y2]
    end
end

"""
    edgeroute_nodes(gc; voff = 0.025)

DOCSTRING

# Arguments:
- `gc`: DESCRIPTION
OPTIONAL INPUTS
- `voff = 0.025`: DESCRIPTION
"""
function edgeroute_nodes(gc; voff = 0.025)
    # first setting the outbound params
    npos = node_prop(gc, :pos)
    xpts = [n[1] for n in npos]
    ypts = [n[2] for n in npos]
    for i = 1:nv(gc)
        outn = outneighbors(gc, i)
        inpos = npos[i]
        xi = inpos[1]
        yi = inpos[2]

        yext = yi
        if length(outn) > 1
            yext = yi .+ LinRange(-voff, voff, length(outn))
        end

        yout = ypts[outn]
        xout = xpts[outn]
        sortidx = sortperm(yout)      # indexes 1 - nout ordered from bottom to top

        for k = 1:length(outn)
            idx = sortidx[k]
            n2 = outn[idx]       # idx so lowest first

            y1 = yext[k]

            y = yout[idx]
            x = xout[idx]

            xmid = (xi + x) / 2
            if y == yi
                xpath = [xi, x]
                ypath = [yi, yi]
                set_prop!(gc, i, n2, :path, (xpath, ypath))
                # println("edge $(i) to $(n2) direct")
            else
                xpath = [xi, xmid, xmid, x]
                ypath = [y1, y1, y, y]
                set_prop!(gc, i, n2, :path, (xpath, ypath))
                # println("edge $(i) to $(n2) $xpath")
            end
        end
    end


    # now doing recieving nodes
    for i = 1:nv(gc)
        inn = inneighbors(gc, i)
        inpos = npos[i]

        xi = inpos[1]
        yi = inpos[2]

        if length(inn) == 1
            continue
        end

        yext = yi .+ LinRange(-voff, voff, length(inn))

        yin = ypts[inn]
        xin = xpts[inn]
        sortidx = sortperm(yin)      # indexes 1 - nout ordered from bottom to top

        for k = 1:length(inn)
            idx = sortidx[k]

            nk = inn[idx]            # idx so lowest first

            y1 = yext[k]             # lowest recieving point

            pthx, pthy = get_prop(gc, nk, i, :path)
            # println("Original edge $(nk) to $(i) $pthx $pthy")
            idx2fix = findall(x -> x == yi, pthy)
            pthy[idx2fix] .= y1
            set_prop!(gc, nk, i, :path, (pthx, pthy))
            # println("edge $(nk) to $(i) $pthx $pthy")
        end
    end
end

"""
    set_plot_props!(gc)

DOCSTRING

# Arguments:
- `gc`: DESCRIPTION
"""
function set_plot_props!(gc)
    etd = edge_propdict(gc, :etype; edgekey = :edge)
    nodeTypes = node_propdict(gc, :nodeType)
    set_default_node_prop!(gc, :nshape, :rect)
    set_default_node_prop!(gc, :markeroutline, :black)
    set_default_edge_prop!(gc, :color, :black)
    for e in collect(edges(gc))
        if etd[e] == :flow
            set_prop!(gc, e.src, e.dst, :color, :black)
        elseif etd[e] == :externalheat
            set_prop!(gc, e.src, e.dst, :color, :red)
            nodeTypes[e.dst] == :real ?
            set_props!(gc, e.dst, Dict(:nshape => :hexagon, :markeroutline => :red)) :
            nothing
            nodeTypes[e.src] == :real ?
            set_props!(gc, e.src, Dict(:nshape => :hexagon, :markeroutline => :red)) :
            nothing
        elseif etd[e] == :externalcool
            set_prop!(gc, e.src, e.dst, :color, :blue)
            nodeTypes[e.dst] == :real ? set_props!(gc, e.dst, Dict(:nshape => :octagon, :markeroutline => :blue)) : nothing
            nodeTypes[e.src] == :real ? set_props!(gc, e.src, Dict(:nshape => :octagon, :markeroutline => :blue)) :
            nothing
        elseif etd[e] == :transferheat
            set_prop!(gc, e.src, e.dst, :color, :magenta2)
            nodeTypes[e.dst] == :real ? set_props!(gc, e.dst, Dict(:nshape => :circle, :markeroutline => :magenta2)) : nothing
            nodeTypes[e.src] == :real ? set_props!(gc, e.src, Dict(:nshape => :circle, :markeroutline => :magenta2)) : nothing
        end
    end
end

"""
    plotnode(pos::Union{Point2f, Vector}, plotshape::Symbol; r = 0.5, w = 1.0, h = 1.0, ngon = 5)

DOCSTRING

# Arguments:
- `pos::Union{Point2f, Vector}`: DESCRIPTION
- `plotshape::Symbol`: DESCRIPTION
OPTIONAL INPUTS
- `r = 0.5`: DESCRIPTION
- `w = 1.0`: DESCRIPTION
- `h = 1.0`: DESCRIPTION
- `ngon = 5`: DESCRIPTION
"""
function plotnode(
    pos::Union{Point2f,Vector},
    plotshape::Symbol;
    r = 0.5,
    w = 1.0,
    h = 1.0,
    ngon = 5,
)

    xcenter = pos[1]
    ycenter = pos[2]


    ngondict = Dict([:pentagon => 5, :hexagon => 6, :octagon => 8])


    if plotshape ∈ collect(keys(ngondict))
        ngon = ngondict[plotshape]
        angles = LinRange(0, 2 * pi, ngon + 1)
        xcoords = r .* cos.(angles)
        ycoords = r .* sin.(angles)
        return Plots.Shape(xcoords .+ xcenter, ycoords .+ ycenter)
    end

    if plotshape == :circle
        pts = Plots.partialcircle(0, 2π, 100, r)
        x, y = Plots.unzip(pts)
        x .+= xcenter
        y .+= ycenter
        return Plots.Shape(x, y)
    end


    y1 = ycenter + h / 2
    y2 = ycenter - h / 2
    x1 = xcenter + w / 2
    x2 = xcenter - w / 2
    return Plots.Shape([x1, x1, x2, x2], [y1, y2, y2, y1])
end

"""
    plotplant(g; mode = :patho, nsize = 1.0, compnamesubs, compnameattr, compnamerot = 45, numbering = false, sysnamedict = Dict(), legpad = 0.2, legheight = 0.5, legoffset = 0.5, nodesize = 6, legwid = 2, pathattr, figattr)

DOCSTRING

# Arguments:
- `g`: DESCRIPTION
OPTIONAL INPUTS
- `mode = :patho`: DESCRIPTION
- `nsize = 1.0`: DESCRIPTION
- `compnamesubs`: DESCRIPTION
- `compnameattr`: DESCRIPTION
- `compnamerot = 45`: DESCRIPTION
- `numbering = false`: DESCRIPTION
- `sysnamedict = Dict()`: DESCRIPTION
- `legpad = 0.2`: DESCRIPTION
- `legheight = 0.5`: DESCRIPTION
- `legoffset = 0.5`: DESCRIPTION
- `nodesize = 6`: DESCRIPTION
- `legwid = 2`: DESCRIPTION
- `pathattr`: DESCRIPTION
- `figattr`: DESCRIPTION
"""
function plotplant(
    g;
    mode = :patho,
    nsize = 1.0,
    compnamesubs = ("enfw" => "en\nfw", "_" => "\n", "circulator" => "pump"),
    compnameattr = (:right, :top, 5),
    compnamerot = 45,
    numbering = false,
    sysnamedict = Dict(),
    legpad = 0.2,
    legheight = 0.5,
    legoffset = 0.5,
    nodesize = 6,
    legwid = 2,
    pathattr = (
        linewidth = 1,
        marker = :circle,
        markersize = 0.0,
        markercolor = :red,
        color = :black,
        alpha = 0.7,
        legend = false,
    ),
    figattr = (grid = true, aspect_ratio = :equal, showaxis = false),
)

    GR.setarrowsize(1)

    gc = g

    npos = node_prop(gc, :pos)
    xs = [nn[1] for nn in npos]
    ys = [nn[2] for nn in npos]
    node_prop(gc, :displayName)

    nv(gc) == length(xs) == length(ys) || error("need 1 position per vertex")

    # reference vectors to prevent functino call backs
    nodenames = node_propdict(gc, :name)
    nodeDisplayNames = node_propdict(gc, :displayName)
    nodeTypes = node_propdict(gc, :nodeType)
    nodePos = node_propdict(gc, :pos)
    nodeHts = node_propdict(gc, :height)
    nodeWds = node_propdict(gc, :width)
    edgeRev = edge_propdict(gc, :isReversed; edgekey = :edge)
    edgetype = edge_propdict(gc, :etype; edgekey = :edge)
    nodeShps = node_propdict(gc, :nshape)

    soln = get_prop(gc, :soln)

    ecolor      = edge_propdict(gc, :color; edgekey = :edge)
    fn(gc, x)   = get_prop(gc, x, :nodeType) == :real
    real_nodes  = filter_vertices(gc, fn)

    # dots
    p = Plots.scatter(xs, ys; markeralpha = 0)
    xreal = [xs[rn] for rn in real_nodes]
    yreal = [ys[rn] for rn in real_nodes]
    mflowvals = String[]
    ymax = maximum(ys)
    ymin = minimum(ys)
    colopts = [
        "mintcream",
        "lightseagreen",
        "beige",
        "darkseagreen1",
        "ivory3",
        "cadetblue1",
        "tomato1",
        "darkorange1",
        "lightyellow",
        "darkslategray2",
        "gray14",
        "lemonchiffon4",
        "mediumblue",
        "mintcream",
        "purple3",
    ]
    syslabs = get_prop(gc, :system_labels)
    numsys = length(syslabs)
    legheight = minimum([(ymax - ymin) / numsys, legheight])
    legspacing = legheight + legpad
    totheight = legspacing * numsys
    ybot_leg = -totheight / 2
    cix = 1

    # EDGE PLOTS
    for edge in edges(gc)
        if mode == :path
            lxs, lys = get_prop(gc, src(edge), dst(edge), :path)
            if edgeRev[edge] == true
                lxs = reverse(lxs)
                lys = reverse(lys)
            end
            defcolor = ecolor[edge]
            if length(lxs) == 4
                lxarrow = lxs[2:3]
                lyarrow = lys[2:3]
                lxn = [lxarrow[1], (lxarrow[1] + lxarrow[2]) / 2]
                lyn = [lyarrow[1], (lyarrow[1] + lyarrow[2]) / 2]

                # Plots.plot!(lxn, lyn; color = defcolor, arrow = (:closed), linewidth=0.25)

                ul = lxn[2] - lxn[1]
                vl = lyn[2] - lyn[1]

                Plots.quiver!(
                    [lxn[1]],
                    [lyn[1]],
                    quiver = ([ul], [vl]);
                    color = defcolor,
                    arrow = Plots.arrow(:closed, :head, 0, 0),
                    linewidth = 0.25,
                )

                par = get_prop(gc, edge.dst, :parent)
                edn = get_prop(gc, edge, :sectionId)
                ned = get_prop(gc, edge, :nSections)

                var = get_prop(gc, edge, :evar)
                unstr = "\nMW"
                txt = string(par) * string(round((soln(var) / 10^6); digits = 0)) * unstr
                showtxt = string(Int(round((soln(var) / 10^6); digits = 0))) * unstr

                if edgetype[edge] == :flow
                    var = get_prop(gc, edge, :mflow)
                    unstr = "\nkg/s"
                    txt =
                        string(par) *
                        string(Int(round(abs(soln(var)); digits = 0))) *
                        "\nkg/s"
                    showtxt = string(round(abs(soln(var)); digits = 0)) * "\nkg/s"
                end

                # if (txt ∈ mflowvals) == false
                if edn < ned || abs(lxn[2] - lxn[1]) > nsize || abs(lyn[2] - lyn[1]) > nsize
                    # if (txt ∈ mflowvals) == false
                        Plots.plot!(;
                            annotations = (
                                [lxn[2] + 0.8],
                                [lyn[2] - nsize/1.5],
                                Plots.text(showtxt, :left, 4, rotation = 0),
                            ),
                        )
                        push!(mflowvals, txt)
                    # end
                end
            end
            Plots.plot!(lxs, lys; color = defcolor, pathattr...)
        else
            if edgeRev[edge] == true
                lxs = reverse(lxs)
                lys = reverse(lys)
            end
            lxs = [xs[src(edge)], xs[dst(edge)]]
            lys = [ys[src(edge)], ys[dst(edge)]]
            Plots.plot!(lxs, lys; pathattr...)
        end
    end

    # fned(g,e) = (get_prop(g,e,:etype) == :transferheat || get_prop(g,e,:etype) == :externalcool)
    # for edge in filter_edges(gc,fned)
    #     edn = get_prop(gc,edge,:sectionId)
    #     if edn == 1
    #         lxs,lys = get_prop(gc,src(edge),dst(edge),:path)
    #         var = get_prop(gc,edge,:evar)
    #         unstr = "\nMW"
    #         showtxt = string(round((soln(var)/10^6); digits = 0)) * unstr
    #         lxM = mean(lxs)
    #         lyM = mean(lys)
    #         Plots.plot!(; annotations = ([lxM], [lyM],Plots.text(showtxt, :center, :center, 4,rotation =0)))
    #     end
    # end


    ## plotting components 
    for sl in syslabs
        fnn(g, n) = string(get_prop(g, n, :parent)) == sl
        nvert = filter_vertices(gc, fnn)
        color = Symbol(colopts[cix])
        for n in nvert
            if nodeTypes[n] == :real
                outlinecolor = get_prop(gc, n, :markeroutline)
                outlinestroke = outlinecolor == :black ? 1.0 : 3.0
                # println(nodeShps[n],get_prop(gc, n, :markeroutline))
                shp = plotnode(nodePos[n], nodeShps[n]; w = nsize, r = nsize/2, h = nsize)
                Plots.plot!(shp;
                    color = color,
                    linecolor = get_prop(gc, n, :markeroutline),
                    stroke = outlinestroke,
                )
                #    0 Plots.scatter!(nodePos[n]; markeralpha=1.0, shape=nodeShps[n],markersize = nodesize, markercolor = color, markerstrokecolor = get_prop(gc,n,:markeroutline), markerstroke = outlinestroke)
            end
        end
        ylow_box = ybot_leg + (cix - 1) * legspacing
        yhi_box = ybot_leg + (cix - 1) * legspacing + legheight
        legbox = Plots.Shape([
            (-legwid - legoffset, ylow_box),
            (-legoffset, ylow_box),
            (-legoffset, yhi_box),
            (-legwid - legoffset, yhi_box),
        ])



        Plots.plot!(legbox; markeralpha = 1.0, fill = color)
        ytext = (ylow_box + yhi_box) / 2

        tstr =
            titlecase(replace(replace(lowercase(string(sl)), compnamesubs...), "\n" => " "))
        if !isempty(sysnamedict)
            tstr = titlecase(replace(sysnamedict[sl], "\n" => " "))
        end
        Plots.scatter!(
            [-legwid + legwid / 2 - legoffset],
            [ytext];
            markersize = markeralpha = 0.0,
            text = [(tstr, :center, :center, 5)],
        )
        cix = cix + 1
    end
    ncols = cix - 1
    if numbering
        Plots.scatter!(
            xs,
            ys;
            markersize = markeralpha = 0.0,
            text = [(i, :center, :center, :red, 5) for i = 1:nv(gc)],
        )
    end
    # labs = Shape([(0.3, 0.0), pts[95], pts[50], (0.3, 0.0)])
    # Plots.scatter!(xreal, yreal; markeralpha=1.0, shape=:rect,markersize = 7.5, markercolor = :white)
    # Plots.scatter!(xreal, yreal.-0.1; text=[(titlecase(replace(get_prop(gc,i,:displayName) , "enfw" => "en\nfw", "_" => "\n", "circulator" => "pump")), :center,:top,  5) for i in real_nodes], markeralpha = 0.0)
    # annotations
    Plots.plot!(;
        annotations = ([
            (
                [xs[i]],
                ys[i] - nsize/1.5,
                Plots.text(
                    titlecase(
                        replace(lowercase(get_prop(gc, i, :displayName)), compnamesubs...),
                    ),
                    rotation = compnamerot,
                    compnameattr...,
                ),
            ) for i in real_nodes
        ]),
    )

    tempLabs = String[]
    sys_pd = node_propdict(gc, :sys)

    for n = 1:nv(gc)
        sys = sys_pd[n]
        if hasproperty(sys, :q) && nodeTypes[n] == :real
            sysq = getproperty(sys, :q)
            if hasproperty(sysq, :Q̇)
                strn = soln(getproperty(sysq, :Q̇))
                unstr = "\nMW"
                showtxt = string(Int(round((strn / 10^6); digits = 0))) * unstr
                Plots.plot!(;
                    annotations = (
                        [xs[n]],
                        [ys[n]],
                        Plots.text(
                            showtxt,
                            :center,
                            :center,
                            Int(maximum([1, nodesize / 2])),
                            rotation = 0,
                        ),
                    ),
                )
            end
        end
        if hasproperty(sys, :w) && nodeTypes[n] == :real
            sysq = getproperty(sys, :w)
            if hasproperty(sysq, :Ẇ)
                strn = soln(getproperty(sysq, :Ẇ))
                unstr = "\nMW"
                showtxt = string(Int(round((strn / 10^6); digits = 0))) * unstr
                Plots.plot!(;
                    annotations = (
                        [xs[n]],
                        [ys[n]],
                        Plots.text(
                            showtxt,
                            :center,
                            :center,
                            Int(maximum([1, nodesize / 2])),
                            rotation = 0,
                        ),
                    ),
                )
            end
        end
        if hasproperty(sys, :hp) && hasproperty(sys, :lp) && nodeTypes[n] == :real
            syshp = getproperty(sys, :hp)
            syslp = getproperty(sys, :lp)

            hpw = getproperty(syshp, :w)
            lpw = getproperty(syslp, :w)

            strn = soln(getproperty(hpw, :Ẇ)) + soln(getproperty(lpw, :Ẇ))
            showtxt = string(Int(round((strn / 10^6); digits = 0))) * "\nMW"
            Plots.plot!(;
                annotations = (
                    [xs[n]],
                    [ys[n]],
                    Plots.text(
                        showtxt,
                        :center,
                        :center,
                        Int(maximum([1, nodesize / 2])),
                        rotation = 0,
                    ),
                ),
            )
        end
        if hasproperty(sys, :ic) && nodeTypes[n] == :real
            sysic = getproperty(sys, :ic)
            sysq = getproperty(sysic, :q)
            strn = soln(getproperty(sysq, :Q̇))
            unstr = "\nMW"
            showtxt = string(Int(round((strn / 10^6); digits = 0))) * unstr
            Plots.plot!(;
                annotations = (
                    [xs[n]],
                    [ys[n]],
                    Plots.text(
                        showtxt,
                        :center,
                        :center,
                        Int(maximum([1, nodesize / 2])),
                        rotation = 0,
                    ),
                ),
            )
        end

        if nodeTypes[n] == :fake

            inn = outneighbors(gc, n)[1]
            if hasproperty(sys, :n) && get_prop(gc, n, inn, :etype) == :flow
                par = get_prop(gc, n, :parent)
                sysp = getproperty(sys, :n)
                sysT = soln(getproperty(sysp, :T))
                sysP = soln(getproperty(sysp, :P))
                showtxt =
                    "T = $(Int(round(sysT-273.18;digits =0 ))) C" *
                    "\n" *
                    "P = $(Int(round(sysP; digits = 0))) Bar"
                txt = string(par) * string(sysT) * string(sysP)

                if (txt ∈ tempLabs) == false
                    Plots.plot!(;
                        annotations = (
                            [xs[n]],
                            [ys[n]],
                            Plots.text(
                                showtxt,
                                :right,
                                :left,
                                Int(maximum([1, nodesize / 2])),
                                rotation = 0,
                            ),
                        ),
                    )
                    push!(tempLabs, txt)
                end
            end
        end
    end

    Plots.plot!(; figattr...)
    display(p)
    return p
end

end
