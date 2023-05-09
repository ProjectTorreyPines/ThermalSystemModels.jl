using ModelingToolkit, Plots, DifferentialEquations, Revise, Unitful, CoolProp,  Logging,NonlinearSolve, Printf
using Revise
using Symbolics

Revise.retry()
@variables t
module Gas
    include("01-ThermoGas.jl")
end

module Steam
    include("01-MultiPhase.jl")
end

module Liq
    include("01-Incompressible.jl")
end

begin
    # Heat Exchangers
    # Liquid to gas
    @component function L2G_HeatExchanger(;name, ϵ = 0.95, A = Liq.IncompressibleHeaatTransfer(), B = Gas.ThermoHeatTransfer())

        ps = @parameters ϵ = ϵ

        @variables Q̇(t)=0.0 Cmin(t)=0.0
        
        eqs = [
            Q̇ ~ ϵ * (A.p.T - B.p.T) * ((A.C < B.C) *  A.C + (A.C >= B.C) * B.C)      #   heat transfer out of A -> B , if A/T > B/T 
            0 ~ A.Q̇ + B.Q̇
            A.Q̇ ~ -Q̇
        ]
        ODESystem(eqs, t, [Q̇], ps; name = name, systems = [A,B], defaults = [ϵ => 0.95] )
    end

    #   Steam to Gas
    @component function S2G_HeatExchanger(;name, ϵ = 0.95, A = Steam.SteamHeatTransfer(), B = Gas.ThermoHeatTransfer() )

        ps = @parameters ϵ = ϵ

        @variables Q̇(t)=0.0 Cmin(t)=0.0
        
        eqs = [
            Q̇ ~ ϵ * (A.p.T - B.p.T) * ((A.C < B.C) *  A.C + (A.C >= B.C) * B.C)      #   heat transfer out of A -> B , if A/T > B/T 
            0 ~ A.Q̇ + B.Q̇
            A.Q̇ ~ -Q̇
        ]
        ODESystem(eqs, t, [Q̇], ps; name = name, systems = [A,B], defaults = [ϵ => 0.95] )
    end

    #   Liquid To Steam
    @component function L2S_HeatExchanger(;name, ϵ = 0.95, A = Liq.IncompressibleHeaatTransfer(), B = Steam.SteamHeatTransfer())

        ps = @parameters ϵ = ϵ

        @variables Q̇(t)=0.0 Cmin(t)=0.0
        
        eqs = [
            Q̇ ~ ϵ * (A.p.T - B.p.T) * ((A.C < B.C) *  A.C + (A.C >= B.C) * B.C)      #   heat transfer out of A -> B , if A/T > B/T 
            0 ~ A.Q̇ + B.Q̇
            A.Q̇ ~ -Q̇
        ]
        ODESystem(eqs, t, [Q̇], ps; name = name, systems = [A,B], defaults = [ϵ => 0.95] )
    end

    @component function Gen_HeatExchanger(;name, ϵ = 0.95, A, B, returnmode = :ode)

        ps = @parameters ϵ = ϵ
        
        eqs = [
            -A.Q̇ ~ ϵ * (A.p.T - B.p.T) * ((A.C < B.C) *  A.C + (A.C >= B.C) * B.C)      #   heat transfer out of A -> B , if A/T > B/T 
            0 ~ A.Q̇ + B.Q̇
        ]
        if returnmode == :ode
            ODESystem(eqs, t, [Q̇], ps; name = name, systems = [A,B], defaults = [ϵ => 0.95] )
        else
            return eqs
        end
    end
end
begin
    @component function FeedwaterRankine(; name, Pmin = 0.1, Pmid = 10, Pmax = 150)
            @named iores        = Steam.ioReservoir(P=Pmin, fixboth = false)
            @named valve        = Steam.SteamFlowValve()
            @named pumpB        = Steam.AdiabaticPump(Pout = Pmax,setpressure = true)
            @named boil         = Steam.SteamHeatTransfer()
            @named turbine      = Steam.SIMOAdiabaticTurbine(setpressure = true, Pyin = Pmid, Pzin = Pmin,ηin = 1.0) 
            @named pumpA        = Steam.AdiabaticPump(Pout = 10, setpressure = true)
            @named condensor         = Steam.IdealCondensor()
            @named openfw       = Steam.OpenFeedwaterHeater()

            connections = vcat(Steam.hydro_connect(pumpB.n,valve.p),
                                Steam.hydro_connect(valve.n, boil.p),
                                Steam.hydro_connect(boil.n, turbine.p),
                                Steam.hydro_connect(turbine.hp.n, openfw.y),
                                Steam.hydro_connect(turbine.lp.n, condensor.p),
                                Steam.hydro_connect(condensor.n,iores.p),
                                Steam.hydro_connect(iores.n,pumpA.p),
                                Steam.hydro_connect(pumpA.n,openfw.z),
                                Steam.hydro_connect(openfw.n,pumpB.p))

            systems =[valve,boil,turbine,pumpB,iores,pumpA,condensor,openfw]
            ODESystem(connections,t;name = name, systems = systems)
            # smodel = substitute(odemodel,hydro_prop_dict)
            # simple_mod = structural_simplify(smodel)
            # prob = ODAEProblem(simple_mod,Pair[],(0.0,1.0))

            # plotmod(smodel)
            # prob.f.observed(openfw.n.P,prob.u0,prob.p,1.0)

            # pf(xin) = prob.f.observed(xin,prob.u0,prob.p,1.0)
            # nodes = [turbine.p, turbine.yn, turbine.zn, turbine.hp.p, turbine.lp.p,turbine.hp.n, turbine.lp.n,  cond.p, cond.n, pumpA.n,pumpB.p, openfw.z, openfw.y]
            # for n in nodes
            #     @printf "\t %s ṁ = %.2f \t P = % .2f \t h = %.2f\n" n.name pf(n.ṁ) pf.(n.P) pf(n.h)
            # end
    end

    @component function ComplexBraytonRegen(;name)
            TminCycle = 300
            PminCycle = 15
            @named res     = Gas.TwoPortReservoir(P = PminCycle, T = TminCycle)
            @named valve   = Gas.GasFlowValve()
            @named comp1   = Gas.ActiveThermoCompressor(rp = 1.7, η = 0.9)
            @named ic1     = Gas.Intercooler(Tout = TminCycle)
            @named comp2   = Gas.ActiveThermoCompressor(rp = 1.5, η = 0.95)
            @named ic2     = Gas.Intercooler(Tout = TminCycle)
            @named comp3   = Gas.ActiveThermoCompressor(rp = 1.5, η = 0.95)
            @named regen   = Gas.Regenerator()
            @named HeatIn  = Gas.ThermoHeatTransfer()
            @named turbine = Gas.PassiveThermoTurbine()
            @named cool    = Gas.IdealCooler()
            
            connections = vcat(Gas.gas_connect(res.n,valve.p),
                Gas.gas_connect(valve.n,comp1.p),
                Gas.gas_connect(comp1.n,ic1.p),
                Gas.gas_connect(ic1.n,comp2.p),
                Gas.gas_connect(comp2.n,ic2.p),
                Gas.gas_connect(ic2.n,comp3.p),
                Gas.hx_connect(regen,comp3,HeatIn,turbine,cool),
                Gas.gas_connect(HeatIn.n,turbine.p),
                Gas.gas_connect(cool.n,res.p));

        
        systemNames = [res,valve,comp1,ic1,comp2,ic2,comp3,regen,HeatIn,turbine,cool];
        ODESystem(connections,t; name = name, systems = systemNames)
    end
    @component function Water_loop(; name, Pmin=32, Pmax=40)
        @named res          = Steam.ioReservoir(P = Pmin, fixboth = true)
        @named valve        = Steam.SteamFlowValve()
        @named pump         = Steam.AdiabaticPump(η = 0.9, setpressure = true, Pout = Pmax)
        @named pset         = Steam.SetPressure(P=Pmin)
        @named HeatIn       = Steam.SteamHeatTransfer()
        @named HeatTx       = Steam.SteamHeatTransfer()
        @named HeatRej      = Steam.IdealCondensor()
        @named throttle     = Steam.throttle()

        connections = vcat(Steam.hydro_connect(res.n,valve.p),
                        Steam.hydro_connect(valve.n, pump.p),
                        Steam.hydro_connect(pump.n,HeatIn.p),
                        Steam.hydro_connect(HeatIn.n,HeatTx.p),
                        Steam.hydro_connect(HeatTx.n,throttle.p),
                        Steam.hydro_connect(throttle.n,HeatRej.p),
                        Steam.hydro_connect(HeatRej.n,res.p))

        systemNames = [res,valve,pump,HeatIn,HeatTx,HeatRej,throttle]
        ODESystem(connections,t; name = name, systems =  systemNames)
    end

    @component function He_loop(; name, Pmin=80,Tmin = 300, Pmax=85)
        @named res          = Gas.TwoPortReservoir(P = Pmin, T = Tmin)
        @named valve        = Gas.GasFlowValve()
        @named circulator   = Gas.PassiveThermoCompressor(η = 0.9)
        @named pset         = Gas.SetPressure(P=Pmax)
        @named HeatIn       = Gas.ThermoHeatTransfer()
        @named HeatTx       = Gas.ThermoHeatTransfer()
        @named HeatRej      = Gas.ThermoHeatTransfer()
        @named throttle     = Gas.throttle()

        connections = vcat(Gas.gas_connect(res.n,valve.p),
                        Gas.gas_connect(valve.n, circulator.p),
                        Gas.gas_connect(circulator.n,pset.p),
                        Gas.gas_connect(pset.n,HeatIn.p),
                        Gas.gas_connect(HeatIn.n,HeatTx.p),
                        Gas.gas_connect(HeatTx.n,HeatRej.p),
                        Gas.gas_connect(HeatRej.n,throttle.p),
                        Gas.gas_connect(throttle.n,res.p))

        systemNames = [res,valve,circulator,pset,HeatIn,HeatTx,HeatRej,throttle]
        ODESystem(connections,t;name =name, systems =  systemNames)
    end

    @component function He_inter_loop(; name, Pmin=80,Tmin = 300, Pmax=85)
        @named res          = Gas.TwoPortReservoir(P = Pmin, T = Tmin)
        @named valve        = Gas.GasFlowValve()
        @named circulator   = Gas.PassiveThermoCompressor(η = 0.9)
        @named pset         = Gas.SetPressure(P=Pmax)
        @named HeatInA      = Gas.ThermoHeatTransfer(ΔP = 0.1)
        @named HeatInB      = Gas.ThermoHeatTransfer(ΔP = 0.1)
        @named HeatInC      = Gas.ThermoHeatTransfer(ΔP = 0.0)
        @named HeatTx       = Gas.ThermoHeatTransfer(ΔP = 0.0)
        @named throttle     = Gas.throttle()
        @named HeatRej      = Gas.ThermoHeatTransfer()

        connections = vcat(Gas.gas_connect(res.n,valve.p),
                        Gas.gas_connect(valve.n, circulator.p),
                        Gas.gas_connect(circulator.n,pset.p),
                        Gas.gas_connect(pset.n,HeatInA.p),
                        Gas.gas_connect(HeatInA.n, HeatInB.p),
                        Gas.gas_connect(HeatInB.n, HeatInC.p),
                        Gas.gas_connect(HeatInC.n, HeatTx.p),
                        Gas.gas_connect(HeatTx.n,   throttle.p),
                        Gas.gas_connect(throttle.n, HeatRej.p),
                        Gas.gas_connect(HeatRej.n, res.p))

        systemNames = [res,valve,circulator,pset,HeatInA,HeatInB,HeatInC,HeatTx,HeatRej,throttle]
        ODESystem(connections,t;name =name, systems =  systemNames)
    end

    @component function breeder_loop(; name, Pmin=32,Pmax=40, Tmin=600)
            @named res      = Liq.TwoPortReservoir(P=Pmin,T=Tmin)
            @named valve    = Liq.IncompressibleFlowValve()
            @named pump     = Liq.PassiveIncompressiblePump()
            @named pset     = Liq.SetPressure(P=Pmax)
            @named HeatIn   = Liq.IncompressibleHeaatTransfer()
            @named HeatTx   = Liq.IncompressibleHeaatTransfer()
            @named HeatRej  = Liq.IdealCooler()
            @named throttle = Liq.throttle()

        
            connections = vcat(Liq.incompressible_connect(res.n,valve.p),
                            Liq.incompressible_connect(valve.n, pump.p),
                            Liq.incompressible_connect(pump.n,pset.p),
                            Liq.incompressible_connect(pset.n,HeatIn.p),
                            Liq.incompressible_connect(HeatIn.n,HeatTx.p),
                            Liq.incompressible_connect(HeatTx.n,HeatRej.p),
                            Liq.incompressible_connect(HeatRej.n,throttle.p),
                            Liq.incompressible_connect(throttle.n,res.p))
            

            systemNames = [res,valve,pump,pset,HeatIn,HeatTx,HeatRej,throttle]
            ODESystem(connections,t;name = name, systems =  systemNames)
    end
end


# HELIUM INTER LOOP
begin
    @named odesys = He_inter_loop(Pmin = 32);

    masfcn(t) = 100 + 15 * sin(t)
    Qf(t) = 100e6 + 35e6 * cos(t)  

    # @unpack valve, HeatIn, HeatTx, HeatRej,circulator = odesys
    @unpack valve, HeatInA, HeatInC,HeatInB, HeatTx, HeatRej, circulator= odesys
    controlled_eqs = [valve.p.ṁ ~ masfcn(t),
                        HeatInA.Q̇ ~ Qf(t),
                        HeatInB.Q̇ ~ .5*Qf(t),
                        HeatInC.Q̇ ~ 0.25*Qf(t),
                        HeatTx.Q̇ ~ 0]

    @named aux_sys = ODESystem(controlled_eqs,t)
    @named odesys = extend(odesys, aux_sys)

    odesys = substitute(odesys,Gas.propDict)

    simple_sys = structural_simplify(odesys)

    saveat = LinRange(1.0,9.0,25)

    kwargs = (abstol=1e-10, reltol=1e-2, saveat = saveat)
    tspan = (0.0,10.0)
    states(simple_sys)


    prob    = ODEProblem(simple_sys,[],tspan,)
    sol   = solve(prob, Rodas4(); kwargs...);

    # c   = [res,valve,comp1,reg.A,heat,turbine,reg.B,cool];
    p1 = plot(sol, vars = [valve.p.ṁ], title = "Cycle Flow Rate")
    p2 = plot(sol, vars = [HeatInA.Q̇, HeatInB.Q̇, HeatInC.Q̇], title = "Heat Loading")
    p3 = plot(sol,vars = [circulator.w.Ẇ])

    #     # c   = [res,valve,comp1,reg.A,heat,turbine,reg.B,cool];

    #     # p3 = plot(sol,vars = [pump.w.Ẇ])
    p=plot(p1,p2,p3,layout = (3,1))
end

# HELIUM LOOP
begin
    @named odesys = He_loop(Pmin = 32);

    masfcn(t) = 100 + 15 * sin(t)
    Qf(t) = 100e6 + 35e6 * cos(t)  

    # @unpack valve, HeatIn, HeatTx, HeatRej,circulator = odesys
    @unpack valve, HeatIn, HeatTx, HeatRej, circulator= odesys
    controlled_eqs = [valve.p.ṁ ~ masfcn(t),
                        HeatIn.Q̇ ~ Qf(t),
                        HeatTx.Q̇ ~ 0]

    @named aux_sys = ODESystem(controlled_eqs,t)
    @named odesys = extend(odesys, aux_sys)

    odesys = substitute(odesys,Gas.propDict)

    simple_sys = structural_simplify(odesys)

    saveat = LinRange(1.0,9.0,25)

    kwargs = (abstol=1e-10, reltol=1e-2, saveat = saveat)
    tspan = (0.0,10.0)
    states(simple_sys)


    prob    = ODEProblem(simple_sys,[],tspan,)
    sol   = solve(prob, Rodas4(); kwargs...);

    # c   = [res,valve,comp1,reg.A,heat,turbine,reg.B,cool];
    p1 = plot(sol, vars = [valve.p.ṁ], title = "Cycle Flow Rate")
    p2 = plot(sol, vars = [HeatIn.Q̇, HeatRej.Q̇], title = "Heat Loading")
    p3 = plot(sol,vars = [circulator.w.Ẇ])

    #     # c   = [res,valve,comp1,reg.A,heat,turbine,reg.B,cool];

    #     # p3 = plot(sol,vars = [pump.w.Ẇ])
    p=plot(p1,p2,p3,layout = (3,1))
end

#       WATER LOOP TEST
begin
    @named odesys = Water_loop(Pmin = 32);

    masfcn(t) = 100 + 15 * sin(t)
    Qf(t) = 100e6 + 35e6 * cos(t)  

    # @unpack valve, HeatIn, HeatTx, HeatRej,circulator = odesys
    @unpack valve, HeatIn, HeatTx, HeatRej, pump = odesys
    controlled_eqs = [valve.p.ṁ ~ masfcn(t),
                        HeatIn.Q̇ ~ Qf(t),
                        HeatTx.Q̇ ~ 0]

    @named aux_sys = ODESystem(controlled_eqs,t)
    @named odesys = extend(odesys, aux_sys)

    odesys = substitute(odesys,Steam.hydro_prop_dict)

    simple_sys = structural_simplify(odesys)

    saveat = LinRange(1.0,9.0,25)

    kwargs = (abstol=1e-10, reltol=1e-2, saveat = saveat)
    tspan = (0.0,10.0)
    states(simple_sys)


    prob    = ODEProblem(simple_sys,[],tspan,)
    sol   = solve(prob, Rodas4(); kwargs...);

    # c   = [res,valve,comp1,reg.A,heat,turbine,reg.B,cool];
    # p1 = plot(sol, vars = [valve.p.ṁ], title = "Cycle Flow Rate")
    # p2 = plot(sol, vars = [HeatIn.Q̇, HeatRej.Q̇], title = "Heat Loading")
    # p3 = plot(sol,vars = [pump.w.Ẇ])

    #     # c   = [res,valve,comp1,reg.A,heat,turbine,reg.B,cool];
    p1 = plot(sol.t, sol[valve.p.ṁ], title = "Cycle Flow Rate")
    p2 = plot(sol.t, sol[HeatIn.Q̇], title = "Heat Loading")
    p3 = plot(sol.t, sol[HeatRej.q.Q̇], title = "Heat Waste")

    #     # p3 = plot(sol,vars = [pump.w.Ẇ])
    p=plot(p1,p2,p3,layout = (3,1))
end

#      Breeder LOOP
begin
    @named odesys = breeder_loop(Pmin = 32);

    masfcn(t) = 100 + 15 * sin(t)
    Qf(t) = 100e6 + 35e6 * cos(t)  

    # @unpack valve, HeatIn, HeatTx, HeatRej,circulator = odesys
    @unpack valve, HeatIn, HeatTx, HeatRej, pump = odesys
    controlled_eqs = [valve.p.ṁ ~ masfcn(t),
                        HeatIn.Q̇ ~ Qf(t),
                        HeatTx.Q̇ ~ 0]

    @named aux_sys = ODESystem(controlled_eqs,t)
    @named odesys = extend(odesys, aux_sys)

    odesys = substitute(odesys,Liq.propDict)

    simple_sys = structural_simplify(odesys)

    saveat = LinRange(1.0,9.0,25)

    kwargs = (abstol=1e-10, reltol=1e-2, saveat = saveat)
    tspan = (0.0,10.0)
    states(simple_sys)


    prob    = ODEProblem(simple_sys,[],tspan,)
    sol   = solve(prob, Rodas4(); kwargs...);

    # c   = [res,valve,comp1,reg.A,heat,turbine,reg.B,cool];
    # p1 = plot(sol, vars = [valve.p.ṁ], title = "Cycle Flow Rate")
    # p2 = plot(sol, vars = [HeatIn.Q̇, HeatRej.Q̇], title = "Heat Loading")
    # p3 = plot(sol,vars = [pump.w.Ẇ])

    #     # c   = [res,valve,comp1,reg.A,heat,turbine,reg.B,cool];
    p1 = plot(sol.t, sol[valve.p.ṁ], title = "Cycle Flow Rate")
    p2 = plot(sol.t, sol[HeatIn.Q̇], title = "Heat Loading")
    p3 = plot(sol.t, sol[HeatRej.q.Q̇], title = "Heat Waste")

    #     # p3 = plot(sol,vars = [pump.w.Ẇ])
    p=plot(p1,p2,p3,layout = (3,1))
end

#     Rankine
begin
    @named odesys = FeedwaterRankine(Pmin = 32);

    masfcn(t) = 100 + 15 * sin(t)
    Qf(t) = 100e6 + 35e6 * cos(t)  

    # @unpack valve, HeatIn, HeatTx, HeatRej,circulator = odesys
    @unpack valve, boil, condensor, pumpA = odesys
    controlled_eqs = [valve.p.ṁ ~ masfcn(t),
                        boil.Q̇ ~ Qf(t)]

    @named aux_sys = ODESystem(controlled_eqs,t)
    @named odesys = extend(odesys, aux_sys)

    odesys = substitute(odesys,Steam.hydro_prop_dict)

    simple_sys = structural_simplify(odesys)

    saveat = LinRange(1.0,9.0,25)

    kwargs = (abstol=1e-10, reltol=1e-2, saveat = saveat)
    tspan = (0.0,10.0)
    states(simple_sys)


    prob    = ODEProblem(simple_sys,[],tspan,)
    sol   = solve(prob, Rodas4(); kwargs...);

    # c   = [res,valve,comp1,reg.A,heat,turbine,reg.B,cool];
    p1 = plot(sol, vars = [valve.p.ṁ], title = "Cycle Flow Rate")
    p2 = plot(sol, vars = [boil.Q̇, condensor.Q̇], title = "Heat Loading")
    p3 = plot(sol,vars = [pumpA.w.Ẇ])

    #     # p3 = plot(sol,vars = [pump.w.Ẇ])
    p=plot(p1,p2,p3,layout = (3,1))
end

#   Brayton
begin
    @named odesys = ComplexBraytonRegen()

    masfcn(t) = 100 + 15 * sin(t)
    Qf(t) = 100e6 + 35e6 * cos(t)  

    @unpack valve,turbine,regen, HeatIn = odesys
    controlled_eqs = [valve.p.ṁ ~ masfcn(t),
                    HeatIn.Q̇ ~ Qf(t)]

    @named aux_sys = ODESystem(controlled_eqs,t)
    @named odesys = extend(odesys, aux_sys)

    odesys = substitute(odesys, Gas.propDict)
    simple_sys = structural_simplify(odesys)


    kwargs = (abstol=1e-10, reltol=1e-2)
    tspan = (0.0,10.0)

    u0 = [valve.p.ṁ => 100
        turbine.p.T=> 1000
        HeatIn.p.cp => 5192
        turbine.rp=>3.5]

    prob = ODEProblem(simple_sys,u0,tspan,)
    sol   = solve(prob, Rodas4(); kwargs...);

    # c   = [res,valve,comp1,reg.A,heat,turbine,reg.B,cool];
    p1 = plot(sol,vars = [valve.p.ṁ], title = "Cycle Flow Rate")
    p2 = plot(sol, vars = [ HeatIn.Q̇], title = "Heat Loading")
    p3 = plot(sol,vars = [HeatIn.n.T, regen.B.n.T, regen.B.p.T, HeatIn.p.T,HeatIn.n.T])
    p=plot(p1,p2,p3,layout = (3,1))
end


##  MONEY -------------BRAYTON IS WORKING
begin
    Revise.retry()
    @named breeder_circuit  = breeder_loop(Tmin = 700);
    @named divertor_circuit = He_loop(Tmin = 300);
    @named wall_circuit     = He_loop(Tmin = 500);
    @named intloop          = He_inter_loop(Tmin = 250);
    @named powercycle       = ComplexBraytonRegen()

    massfcn(t) = 140 + 15 * sin(t)
    Qloads = [700e6 + 100e6*cos(2t), 150e6 + 10e6*cos(2t), 140e6 + 30e6 * sin(.5*t)]
    auxeq = Equation[]

    for (i,sys) in enumerate([divertor_circuit, wall_circuit, breeder_circuit])
        push!(auxeq,sys.valve.p.ṁ ~ 2*i*massfcn(t))
        push!(auxeq,sys.HeatIn.Q̇ ~ Qloads[i])
    end

    for (i,sys) in enumerate([powercycle.valve, intloop.valve])
        push!(auxeq,sys.p.ṁ ~ massfcn(t))
    end

    @named breeder_hx   = Gen_HeatExchanger(A = breeder_circuit.HeatTx, B = intloop.HeatInC, returnmode = :eq)
    @named wall_hx      = Gen_HeatExchanger(A = wall_circuit.HeatTx, B = intloop.HeatInB, returnmode =:eq)
    @named divertor_hx  = Gen_HeatExchanger(A = divertor_circuit.HeatTx, B = intloop.HeatInA,returnmode =:eq)
    @named primary_hx   = Gen_HeatExchanger(A = powercycle.HeatIn, B = intloop.HeatTx, returnmode =:eq)

    for eq in [breeder_hx, wall_hx, divertor_hx, primary_hx]
        push!(auxeq,eq...)
    end

    @named Plant = ODESystem(auxeq,t; name = :Plant, systems = [powercycle,wall_circuit,breeder_circuit,divertor_circuit,intloop])

    for pd in [Gas.propDict,Steam.hydro_prop_dict,Liq.propDict]
        Plant = substitute(Plant, pd)
    end
    Plant = alias_elimination(Plant)
    SimplePlant = structural_simplify(Plant)
    # kwargs = (abstol=1e-10, reltol=1e-2)
    # tspan = (0.0,10.0)
    prob = ODEProblem(SimplePlant,[],tspan,)
    sol   = solve(prob, Rodas4(); kwargs...);
    refarr = [powercycle,wall_circuit,breeder_circuit,divertor_circuit,intloop];
    fusarr = [breeder_circuit,divertor_circuit,wall_circuit]
    hxarr = [breeder_circuit,divertor_circuit,wall_circuit, intloop]
    # c   = [res,valve,comp1,reg.A,heat,turbine,reg.B,cool];
    # p1 = plot(sol,vars = [powercycle.valve.p.ṁ, breeder_circuit.valve.p.ṁ, wall_circuit.valve.p.ṁ, divertor_circuit.valve.p.ṁ], title = "Cycle Flow Rate")
    
    
    # [powercycle,wall_circuit,breeder_circuit,divertor_circuit,intloop]
    p1 = plot(sol, vars = [r.valve.p.ṁ for r in refarr], title = "Mass Flow Rate")
    p2 = plot(sol, vars = [r.HeatIn.Q̇ for r in fusarr], title = "Heat Load (MW)")
    p3 = plot(sol, vars = [r.HeatTx.Q̇ for r in hxarr], title = "Heat Coupling Rate")
    p=plot(p1,p2,p3,layout = (3,1))
end


##  MONEY -------------RANKINE IS WORKING
begin
    Revise.retry()
    @named breeder_circuit  = breeder_loop(Tmin = 700);
    @named divertor_circuit = He_loop(Tmin = 300);
    @named wall_circuit     = He_loop(Tmin = 500);
    @named intloop          = He_inter_loop(Tmin = 250);
    @named powercycle       = FeedwaterRankine()

    massfcn(t) = 140 + 15 * sin(t)
    Qloads = [700e6 + 100e6*cos(2t), 150e6 + 10e6*cos(2t), 140e6 + 30e6 * sin(.5*t)]
    auxeq = Equation[]

    for (i,sys) in enumerate([divertor_circuit, wall_circuit, breeder_circuit])
        push!(auxeq,sys.valve.p.ṁ ~ 2*i*massfcn(t))
        push!(auxeq,sys.HeatIn.Q̇ ~ Qloads[i])
    end

    for (i,sys) in enumerate([powercycle.valve, intloop.valve])
        push!(auxeq,sys.p.ṁ ~ massfcn(t))
    end

    @named breeder_hx   = Gen_HeatExchanger(A = breeder_circuit.HeatTx, B = intloop.HeatInC, returnmode = :eq)
    @named wall_hx      = Gen_HeatExchanger(A = wall_circuit.HeatTx, B = intloop.HeatInB, returnmode =:eq)
    @named divertor_hx  = Gen_HeatExchanger(A = divertor_circuit.HeatTx, B = intloop.HeatInA,returnmode =:eq)
    @named primary_hx   = Gen_HeatExchanger(A = powercycle.boil, B = intloop.HeatTx, returnmode =:eq)

    for eq in [breeder_hx, wall_hx, divertor_hx, primary_hx]
        push!(auxeq,eq...)
    end

    @named Plant = ODESystem(auxeq,t; name = :Plant, systems = [powercycle,wall_circuit,breeder_circuit,divertor_circuit,intloop])

    for pd in [Gas.propDict,Steam.hydro_prop_dict,Liq.propDict]
        Plant = substitute(Plant, pd)
    end

    Plant = alias_elimination(Plant)
    SimplePlant = structural_simplify(Plant)

    tspan = (0.0,10.0)

    saveat = LinRange(1.0,9.0,25)

    kwargs = (abstol=1e-10, reltol=1e-2, saveat =saveat)
    prob = ODAEProblem(SimplePlant,[],tspan,)
    sol   = solve(prob, Rodas4(); kwargs...);

    refarr = [powercycle,wall_circuit,breeder_circuit,divertor_circuit,intloop];
    fusarr = [breeder_circuit,divertor_circuit,wall_circuit]
    hxarr   = [breeder_circuit,divertor_circuit,wall_circuit, intloop]

    
    # [powercycle,wall_circuit,breeder_circuit,divertor_circuit,intloop]
    p1 = plot(sol, vars = [r.valve.p.ṁ for r in refarr], title = "Mass Flow Rate")
    p2 = plot(sol, vars = [r.HeatIn.Q̇ for r in fusarr], title = "Heat Load (MW)")
    p3 = plot(sol, vars = [powercycle.boil.Q̇], title = "Boiler Qin")

    p4 = plot(sol, vars = [powercycle.turbine.w.Ẇ], title = "Turbine ")
    p=plot!(p1,p2,p3,layout = (3,1))
end

