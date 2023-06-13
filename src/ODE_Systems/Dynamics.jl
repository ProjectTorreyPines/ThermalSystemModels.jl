module Dynamics
    using ModelingToolkit, Plots, DifferentialEquations, Revise, Unitful, CoolProp,  Logging,NonlinearSolve, Printf
    using Revise, OrdinaryDiffEq
    using Symbolics
    @variables t
    Logging.disable_logging(Logging.Warn)
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
    
    @component function L2G_HeatExchanger(;name, ϵ = 0.95, A = Liq.IncompressibleHeaatTransfer(), B = Gas.ThermoHeatTransfer())

        ps = @parameters ϵ = ϵ
    
        @variables Q̇(t)=0.0 Cmin(t)=0.0
        
        eqs = [
            Q̇ ~ ϵ * (A.p.T - B.p.T) * ((A.C < B.C) *  A.C + (A.C >= B.C) * B.C)      #heat transfer out of A -> B , if A/T > B/T 
            0 ~ A.Q̇ + B.Q̇
            A.Q̇ ~ -Q̇
        ]
        ODESystem(eqs, t, [Q̇], ps; name = name, systems = [A,B], defaults = [ϵ => 0.95])
    end
    
    #   Steam to Gas
    @component function S2G_HeatExchanger(;name, ϵ = 0.95, A = Steam.SteamHeatTransfer(), B = Gas.ThermoHeatTransfer(),returnmode = :eqs)
    
        ps = @parameters ϵ = ϵ
        
        eqs = [
            ΔTdes = ((B.p.T - A.p.T) * ϵ)                                                   # max temperature change by gas, + if heat from gas to steam
            A.Q̇ ~ ((A.p.ṁ * (Steam.stm_hptfunc(A.p.P, A.p.T + ((B.p.T - A.p.T) * ϵ))-A.p.h)) <= (B.C * ((B.p.T - A.p.T) * ϵ))) *(A.p.ṁ * (Steam.stm_hptfunc(A.p.P, A.p.T + ((B.p.T - A.p.T) * ϵ)  )-A.p.h)) + ((B.C * ((B.p.T - A.p.T) * ϵ)  ) <= (A.p.ṁ * (Steam.stm_hptfunc(A.p.P, A.p.T + ((B.p.T - A.p.T) * ϵ)  )-A.p.h))) * (B.C * ((B.p.T - A.p.T) * ϵ)  )     # heat transfer out of A -> B , if A/T > B/T 
            0 ~ A.Q̇ + B.Q̇
        ]

        if returnmode == :ode
            ODESystem(eqs, t, [], ps; name = name, systems = [A,B], defaults = [ϵ => 0.95] )
        else
            return eqs
        end
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
            ODESystem(eqs, t, [], ps; name = name, systems = [A,B], defaults = [ϵ => 0.95] )
        else
            return eqs
        end
    end
    
    @component function FeedwaterRankine2(; name, Pmin = 0.1, Pmid = 10, Pmax = 150)
        # Control elements
        @named gnd          = Steam.ContinuityReservoir()
        @named valve        = Steam.SteamFlowValve()
        @named pumpA        = Steam.AdiabaticPump(Pout = Pmid, setpressure = false)
        @named pumpB        = Steam.AdiabaticPump(Pout = Pmax, setpressure = true)
    
        #Boiler
        # @named boil         = IdealBoiler(Tout = 600+273)
        @named boil         = Steam.SteamHeatTransfer()
        @named turbn        = Steam.SIMOAdiabaticTurbine(setpressure = true,Pyin = Pmid, Pzin = Pmin, ηin = 1.0) 
    
        @named cndnsr       = Steam.IdealCondensor()
        @named openfw       = Steam.OpenFeedwaterHeater()
    
        @named WorkRes     = Steam.WorkPin()
        @named ColdUtil    = Steam.HeatTransferPin()
        @named HotUtil     = Steam.HeatTransferPin()
    
        connections = vcat(Steam.hydro_connect(openfw.n,valve.p),
                        Steam.hydro_connect(valve.n, pumpB.p),          # pump -> boilder
                        Steam.hydro_connect(pumpB.n,gnd.p),
                        Steam.hydro_connect(gnd.n,boil.p),
                        Steam.hydro_connect(boil.n, turbn.p),             # boiler -> turbine
                        Steam.hydro_connect(turbn.hp.n, openfw.p1),        # turbine -> openfw y
                        Steam.hydro_connect(turbn.lp.n, cndnsr.p),        #   '------> condensor
                        Steam.hydro_connect(cndnsr.n, pumpA.p),
                        Steam.hydro_connect(pumpA.n,openfw.p2),                                
                        work_connect(WorkRes, turbn.lp.w,turbn.hp.w, pumpA.w, pumpB.w),
                        heat_connect(ColdUtil, cndnsr.q),
                        heat_connect(HotUtil, boil.q));
        
        systems = [valve, pumpA, pumpB, boil, turbn, cndnsr, openfw, gnd, WorkRes, ColdUtil, HotUtil]
    
        ODESystem(connections,t; name = name, systems = systems)
    end
    
    @component function FwRankine3(;name)
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
        @named valve        = Steam.SteamFlowValve()
        @named openfw       = Steam.OpenFeedwaterHeater()
        @named closedfw     = Steam.ClosedFeedwaterHeater()
        @named mixer        = Steam.MixingChamber()

        @named pump1        = Steam.AdiabaticPump(Pout = Pmid1, setpressure = false, η=1.0)
        @named pump2        = Steam.AdiabaticPump(Pout = Pmax, setpressure = false, η=1.0)
        @named pump3        = Steam.AdiabaticPump(Pout = Pmax, setpressure = true, η=1.0, controlinlet = true)

        @named turbine1    = Steam.AdiabaticTurbine(setpressure = true, Pout = Pmid2)
        @named turbine2    = Steam.SIMOAdiabaticTurbine(setpressure = true, Pyin = Pmid1, Pzin = Pmin, ηin = 1.0) 

        @named boiler       = Steam.SteamHeatTransfer()
        @named reheat       = Steam.SteamHeatTransfer()
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
                            work_connect(ElectricGen, turbine1.w, turbine2.hp.w, turbine2.lp.w),
                            work_connect(ElectricUtil,pump1.w, pump2.w, pump3.w),
                            heat_connect(ColdUtil, condensor.q))
                

        systems = [reservoir, valve, openfw, closedfw, mixer, pump1, pump2,pump3,turbine1,turbine2, boiler,reheat,condensor, ElectricGen, ElectricUtil,ColdUtil]
        ODESystem(connections, t;name = name, systems = systems)
    end
    @component function FeedwaterRankine(; name, Pmin = 0.1, Pmid = 10, Pmax = 150)
            @named iores        = Steam.ioReservoir(P=Pmin, fixboth = false)
            @named valve        = Steam.SteamFlowValve()
            @named pumpB        = Steam.AdiabaticPump(Pout = Pmax,setpressure = true)
            @named boil         = Steam.SteamHeatTransfer()
            @named turbine      = Steam.SIMOAdiabaticTurbine(setpressure = true, Pyin = Pmid, Pzin = Pmin,ηin = 1.0) 
            @named pumpA        = Steam.AdiabaticPump(Pout = 10, setpressure = true)
            @named condensor    = Steam.IdealCondensor()
            @named openfw       = Steam.OpenFeedwaterHeater()
    
            @named WorkRes     = Steam.WorkPin()
            @named ColdUtil    = Steam.HeatTransferPin()
            @named HotUtil     = Steam.HeatTransferPin()
    
            connections = vcat(Steam.hydro_connect(pumpB.n,valve.p),
                                Steam.hydro_connect(valve.n, boil.p),
                                Steam.hydro_connect(boil.n, turbine.p),
                                Steam.hydro_connect(turbine.hp.n, openfw.p1),
                                Steam.hydro_connect(turbine.lp.n, condensor.p),
                                Steam.hydro_connect(condensor.n,iores.p),
                                Steam.hydro_connect(iores.n,pumpA.p),
                                Steam.hydro_connect(pumpA.n,openfw.p2),
                                Steam.hydro_connect(openfw.n,pumpB.p),
                                work_connect(WorkRes, turbine.lp.w,turbine.hp.w, pumpA.w, pumpB.w),
                                heat_connect(ColdUtil, condensor.q),
                                heat_connect(HotUtil, boil.q));
    
            systems =[valve,boil,turbine,pumpB,iores,pumpA,condensor,openfw, WorkRes, ColdUtil, HotUtil]
            ODESystem(connections,t;name = name, systems = systems)
    end
    
    @component function ComplexBraytonRegen(;name)
            TminCycle = 300
            PminCycle = 15
            @named WorkRes     = Steam.WorkPin()
            @named ColdUtil    = Steam.HeatTransferPin()
            @named HotUtil     = Steam.HeatTransferPin()
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
                Gas.gas_connect(cool.n,res.p),
                work_connect(WorkRes, turbine.w,comp1.w, comp2.w, comp3.w),
                heat_connect(ColdUtil, cool.q, ic2, ic1),
                heat_connect(HotUtil, HeatIn.q));
    
        
        systemNames = [res,valve,comp1,ic1,comp2,ic2,comp3,regen,HeatIn,turbine,cool, WorkRes, ColdUtil, HotUtil];
        ODESystem(connections,t; name = name, systems = systemNames)
    end
    
    @component function Water_loop(; name, Pmin=32, Pmax=40)
        @named WorkRes     = Steam.WorkPin()
        @named ColdUtil    = Steam.HeatTransferPin()
    
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
                        Steam.hydro_connect(HeatRej.n,res.p),
                        work_connect(WorkRes, pump.w),
                        heat_connect(ColdUtil, HeatRej.q))
    
        systemNames = [res,valve,pump,HeatIn,HeatTx,HeatRej,throttle, WorkRes, ColdUtil]
        ODESystem(connections,t; name = name, systems =  systemNames)
    end
    
    @component function He_loop(; name, Pmin=80,Tmin = 300, Pmax=85)
        @named WorkRes     = Gas.WorkPin()
        @named ColdUtil    = Gas.HeatTransferPin()
    
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
                        Gas.gas_connect(throttle.n,res.p),
                        work_connect(WorkRes, circulator.w),
                        heat_connect(ColdUtil, HeatRej.q))
    
        systemNames = [res,valve,circulator,pset,HeatIn,HeatTx,HeatRej,throttle, WorkRes, ColdUtil]
        ODESystem(connections,t;name =name, systems =  systemNames)
    end
    
    @component function He_inter_loop(; name, Pmin=80,Tmin = 300, Pmax=85)
        @named WorkRes     = Gas.WorkPin()
        @named ColdUtil    = Gas.HeatTransferPin()
        @named res          = Gas.TwoPortReservoir(P = Pmin, T = Tmin)
        @named valve        = Gas.GasFlowValve()
        @named circulator   = Gas.PassiveThermoCompressor(η = 0.9)
        @named pset         = Gas.SetPressure(P=Pmax)
        @named HeatInA      = Gas.ThermoHeatTransfer(ΔP = 0.1)
        @named HeatInB      = Gas.ThermoHeatTransfer(ΔP = 0.1)
        @named HeatInC      = Gas.ThermoHeatTransfer(ΔP = 0.0)
        @named HeatTx       = Gas.ThermoHeatTransfer(ΔP = 0.0)
        @named HeatRej      = Gas.ThermoHeatTransfer()
        @named throttle     = Gas.throttle()
    
        connections = vcat(Gas.gas_connect(res.n,valve.p),
                        Gas.gas_connect(valve.n, circulator.p),
                        Gas.gas_connect(circulator.n,pset.p),
                        Gas.gas_connect(pset.n,HeatInA.p),
                        Gas.gas_connect(HeatInA.n, HeatInB.p),
                        Gas.gas_connect(HeatInB.n, HeatInC.p),
                        Gas.gas_connect(HeatInC.n, HeatTx.p),
                        Gas.gas_connect(HeatTx.n,   HeatRej.p),
                        Gas.gas_connect(HeatRej.n, throttle.p),
                        Gas.gas_connect(throttle.n, res.p),
                        heat_connect(ColdUtil, HeatRej.q),
                        work_connect(WorkRes, circulator.w))
    
        systemNames = [res,valve,circulator,pset,HeatInA,HeatInB,HeatInC,HeatTx,HeatRej,throttle,WorkRes, ColdUtil]
        ODESystem(connections,t;name =name, systems =  systemNames)
    end
    
    @component function He_quad_inter_loop(; name, Pmin=80,Tmin = 300, Pmax=85)
        @named WorkRes     = Gas.WorkPin()
        @named ColdUtil    = Gas.HeatTransferPin()
        @named res          = Gas.TwoPortReservoir(P = Pmin, T = Tmin)
        @named valve        = Gas.GasFlowValve()
        @named circulator   = Gas.PassiveThermoCompressor(η = 0.9)
        @named pset         = Gas.SetPressure(P=Pmax)
        @named HeatInA      = Gas.ThermoHeatTransfer(ΔP = 0.1)
        @named HeatInB      = Gas.ThermoHeatTransfer(ΔP = 0.1)
        @named HeatInC      = Gas.ThermoHeatTransfer(ΔP = 0.0)
        @named boilTx       = Gas.ThermoHeatTransfer(ΔP = 0.0)
        @named reheatTx     = Gas.ThermoHeatTransfer(ΔP = 0.0)
        @named HeatRej      = Gas.ThermoHeatTransfer()
        @named throttle     = Gas.throttle()
    
        connections = vcat(Gas.gas_connect(res.n,valve.p),
                        Gas.gas_connect(valve.n, circulator.p),
                        Gas.gas_connect(circulator.n,pset.p),
                        Gas.gas_connect(pset.n,HeatInA.p),
                        Gas.gas_connect(HeatInA.n, HeatInB.p),
                        Gas.gas_connect(HeatInB.n, HeatInC.p),
                        Gas.gas_connect(HeatInC.n, boilTx.p),
                        Gas.gas_connect(boilTx.n,   reheatTx.p),
                        Gas.gas_connect(reheatTx.n, HeatRej.p),
                        Gas.gas_connect(HeatRej.n, throttle.p),
                        Gas.gas_connect(throttle.n, res.p),
                        heat_connect(ColdUtil, HeatRej.q),
                        work_connect(WorkRes, circulator.w))
    
        systemNames = [res,valve,circulator,pset,HeatInA,HeatInB,HeatInC,boilTx,reheatTx,HeatRej,throttle,WorkRes, ColdUtil]
        ODESystem(connections,t;name =name, systems =  systemNames)
    end

    @component function breeder_loop(; name, Pmin=32,Pmax=40, Tmin=600)
            @named WorkRes     = Gas.WorkPin()
            @named ColdUtil    = Gas.HeatTransferPin()
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
                            Liq.incompressible_connect(throttle.n,res.p),
                            heat_connect(ColdUtil, HeatRej.q),
                            work_connect(WorkRes, pump.w))
    
            
    
            systemNames = [res,valve,pump,pset,HeatIn,HeatTx,HeatRej,throttle, WorkRes, ColdUtil]
            ODESystem(connections,t;name = name, systems =  systemNames)
    end
    
    # HELIUM INTER LOOP
    function test_inter_loop()
        @named odesys = He_inter_loop(Pmin = 32);
    
        masfcn(t) = 100 + 15 * sin(t)
        Qf(t) = 100e6 + 35e6 * cos(t)  
    
        # @unpack valve, HeatIn, HeatTx, HeatRej,circulator = odesys
        @unpack valve, HeatInA, HeatInC,HeatInB, HeatTx, HeatRej, circulator, WorkRes= odesys
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
        p3 = plot(sol,vars = [WorkRes.Ẇ])
    
        #     # c   = [res,valve,comp1,reg.A,heat,turbine,reg.B,cool];
        #     # p3 = plot(sol,vars = [pump.w.Ẇ])
        p=plot(p1,p2,p3,layout = (3,1))
    end
    
    # HELIUM LOOP
    function test_he_loop()
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
        p3 = plot(sol,vars = [WorkRes.Ẇ])
    
        #     # c   = [res,valve,comp1,reg.A,heat,turbine,reg.B,cool];
        #     # p3 = plot(sol,vars = [pump.w.Ẇ])
        p=plot(p1,p2,p3,layout = (3,1))
    end
    
    #       WATER LOOP TEST
    function test_water_loop()
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
    function test_breeder_loop()
        @named odesys = breeder_loop(Pmin = 32);
    
        masfcn(t) = 100 + 15 * sin(t)
        Qf(t) = 100e6 + 35e6 * cos(t)  
    
        # @unpack valve, HeatIn, HeatTx, HeatRej,circulator = odesys
        @unpack valve, HeatIn, HeatTx, HeatRej, pumpA = odesys
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
    function test_feedwater_ranking()
        @named odesys = FeedwaterRankine2(Pmin = 0.1);
        refsys =odesys
    
        masfcn(t) = 100 + 15 * sin(t)
        Qf(t) = 100e6 + 35e6 * cos(t)  
    
        # @unpack valve, HeatIn, HeatTx, HeatRej,circulator = odesys
        @unpack valve, boil,gnd, cndnsr,openfw, pumpA, WorkRes, turbn = odesys
        controlled_eqs = [valve.p.ṁ ~ masfcn(t),
                            boil.Q̇ ~ Qf(t)]
    
        @named aux_sys = ODESystem(controlled_eqs,t)
        @named odesys = extend(odesys, aux_sys)
    
        # odesys = substitute(odesys,Steam.hydro_prop_dict)
        simple_sys = structural_simplify(odesys)
        saveat = LinRange(1.0,9.0,25)
        kwargs = (abstol=1e-10, reltol=1e-2, saveat = saveat)
        tspan = (0.0,10.0)
        states(simple_sys)
        prob    = ODEProblem(simple_sys,[],tspan,)
        sol   = solve(prob, Rodas4(); kwargs...);
        # c   = [res,valve,comp1,reg.A,heat,turbine,reg.B,cool];
        p1 = plot(sol, vars = [valve.p.ṁ], title = "Cycle Flow Rate")
        p2 = plot(sol, vars = [boil.Q̇, cndnsr.Q̇], title = "Heat Loading")
        p3 = plot(sol,vars = [WorkRes.Ẇ], title = "Net Work")
        p4 = plot(sol.t, sol[openfw.n.Φ], title = "PumpA work")
        #     # p3 = plot(sol,vars = [pump.w.Ẇ])
        p=plot(p1,p2,p3,p4,layout = (4,1))
        showsol([valve, boil,gnd, cndnsr,openfw, pumpA, WorkRes, turbn],sol)
    end
    
    function test_brayton()
        @named odesys = ComplexBraytonRegen()
    
        masfcn(t) = 100 + 15 * sin(t)
        Qf(t) = 100e6 + 35e6 * cos(t)  
    
        @unpack res,valve,comp1,ic1,comp2,ic2,comp3,regen,HeatIn,turbine,cool, WorkRes, ColdUtil, HotUtil = odesys
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
    
        # c = [res,valve,comp1,reg.A,heat,turbine,reg.B,cool];
        p1 = plot(sol,vars = [valve.p.ṁ], title = "Cycle Flow Rate")
        p2 = plot(sol, vars = [ HeatIn.Q̇], title = "Heat Loading")
        p3 = plot(sol,vars = [HeatIn.n.T, regen.B.n.T, regen.B.p.T, HeatIn.p.T,HeatIn.n.T])
        p4 = plot(sol, vars = [WorkRes.Ẇ, ColdUtil.Q̇], title = "Net Work")
        p=plot(p1,p2,p3,p4,layout = (4,1))

        p2 = plot(sol, vars = [turbine.p.P, turbine.n.P, turbine.rp])
        display(p2)

        showsol([res,valve,comp1,ic1,comp2,ic2,comp3,regen,HeatIn,turbine,cool, WorkRes, ColdUtil, HotUtil], sol)

    end

    function para_full_system_brayton()
        Revise.retry()
        @named breeder_circuit  = breeder_loop(Tmin = 700);
        @named divertor_circuit = He_loop(Tmin = 300);
        @named wall_circuit     = He_loop(Tmin = 500);
        @named intloop          = He_inter_loop(Tmin = 250);
        @named powercycle       = ComplexBraytonRegen()

        sts = @parameters Qdiv(t) Qwall(t) Qbrd(t) cycle_mass_fcn(t) div_mass_fcn(t) breeder_mass_fcn(t) wall_mass_fcn(t) loop_mass_fcn(t) Qdiv(t) Qwall(t) Qbrd(t)  loop_mass_fcn(t)
        auxeq = [divertor_circuit.valve.p.ṁ ~ div_mass_fcn
                breeder_circuit.valve.p.ṁ ~ breeder_mass_fcn
                powercycle.valve.p.ṁ ~ cycle_mass_fcn
                intloop.valve.p.ṁ ~ loop_mass_fcn
                wall_circuit.valve.p.ṁ ~ wall_mass_fcn
                divertor_circuit.HeatIn.Q̇ ~ Qdiv
                breeder_circuit.HeatIn.Q̇ ~ Qbrd
                wall_circuit.HeatIn.Q̇ ~ Qwall]

        @named breeder_hx   = Gen_HeatExchanger(A = breeder_circuit.HeatTx, B = intloop.HeatInC, returnmode = :eq)
        @named wall_hx      = Gen_HeatExchanger(A = wall_circuit.HeatTx, B = intloop.HeatInB, returnmode =:eq)
        @named divertor_hx  = Gen_HeatExchanger(A = divertor_circuit.HeatTx, B = intloop.HeatInA,returnmode =:eq)
        @named primary_hx   = Gen_HeatExchanger(A = powercycle.HeatIn, B = intloop.HeatTx, returnmode =:eq)
    
        for eq in [breeder_hx, wall_hx, divertor_hx, primary_hx]
            push!(auxeq,eq...)
        end
        # @named cnrl_sys = ODESystem(auxeq,t)
        plant_sys = compose(ODESystem(auxeq,t,[],sts; name = :plant_sys),powercycle,wall_circuit,breeder_circuit,divertor_circuit,intloop)
        return plant_sys,sts
    end

    function para_full_system_rankine()
        Revise.retry()
        @named breeder_circuit  = breeder_loop(Tmin = 700);
        @named divertor_circuit = He_loop(Tmin = 300);
        @named wall_circuit     = He_loop(Tmin = 500);
        @named intloop          = He_quad_inter_loop(Tmin = 250);
        @named powercycle       = FwRankine3()

        sts = @parameters Qdiv=150e6 Qwall=100e6 Qbrd=700e6 cycle_mass_fcn=20 div_mass_fcn=50 breeder_mass_fcn=100 wall_mass_fcn=50 loop_mass_fcn=40
        auxeq = [divertor_circuit.valve.p.ṁ ~ div_mass_fcn
                breeder_circuit.valve.p.ṁ ~ breeder_mass_fcn
                powercycle.valve.p.ṁ ~ cycle_mass_fcn
                intloop.valve.p.ṁ ~ loop_mass_fcn
                wall_circuit.valve.p.ṁ ~ wall_mass_fcn
                divertor_circuit.HeatIn.Q̇ ~ Qdiv
                breeder_circuit.HeatIn.Q̇ ~ Qbrd
                wall_circuit.HeatIn.Q̇ ~ Qwall]

        @named breeder_hx   = Gen_HeatExchanger(A = breeder_circuit.HeatTx, B = intloop.HeatInC, returnmode = :eq)
        @named wall_hx      = Gen_HeatExchanger(A = wall_circuit.HeatTx, B = intloop.HeatInB, returnmode =:eq)
        @named divertor_hx  = Gen_HeatExchanger(A = divertor_circuit.HeatTx, B = intloop.HeatInA,returnmode =:eq)
        @named primary_hx   = Gen_HeatExchanger(A = powercycle.boiler, B = intloop.boilTx, returnmode =:eq)
        @named reheat_hx   = Gen_HeatExchanger(A = powercycle.reheat, B = intloop.reheatTx, returnmode =:eq)
        for eq in [breeder_hx, wall_hx, divertor_hx, primary_hx, reheat_hx]
            push!(auxeq,eq...)
        end
        # @named cnrl_sys = ODESystem(auxeq,t)
        Plant = compose(ODESystem(auxeq,t,[],sts; name = :Plant),powercycle,wall_circuit,breeder_circuit,divertor_circuit,intloop)

        for pd in [Gas.propDict,Liq.propDict]
            Plant = substitute(Plant, pd)
        end
        return Plant,sts
    end

    @component function full_system_brayton(; name)
        Revise.retry()
        @named breeder_circuit  = breeder_loop(Tmin = 700);
        @named divertor_circuit = He_loop(Tmin = 300);
        @named wall_circuit     = He_loop(Tmin = 500);
        @named intloop          = He_inter_loop(Tmin = 250);
        @named powercycle       = ComplexBraytonRegen()

        @variables cycle_mass_fcn(t) div_mass_fcn(t) breeder_mass_fcn(t) wall_mass_fcn(t)
        @variables loop_mass_fcn(t)
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
        ODESystem(auxeq,t; name = name, systems = [powercycle,wall_circuit,breeder_circuit,divertor_circuit,intloop])
    end

    ##  -------------BRAYTON IS WORKING
    function test_full_system_brayton()
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
        
        copyplant = deepcopy(Plant)
        Psys = plantsys(Plant);

        for pd in [Gas.propDict,Steam.hydro_prop_dict,Liq.propDict]
            Plant = substitute(Plant, pd)
        end

        SimplePlant = structural_simplify(Plant)

        kwargs = (abstol=1e-10, reltol=1e-2)
        tspan = (0.0,10.0)
        u0 = [powercycle.turbine.p.T=> 1000
            powercycle.turbine.rp => 3.5]

        prob    = ODEProblem(SimplePlant,u0,tspan,)
        sol     = solve(prob, Rodas4(); kwargs...);

        refarr  = [powercycle,wall_circuit,breeder_circuit,divertor_circuit,intloop];
        fusarr  = [breeder_circuit,divertor_circuit,wall_circuit]
        hxarr   = [breeder_circuit,divertor_circuit,wall_circuit, intloop]
        # c   = [res,valve,comp1,reg.A,heat,turbine,reg.B,cool];
        # p1 = plot(sol,vars = [powercycle.valve.p.ṁ, breeder_circuit.valve.p.ṁ, wall_circuit.valve.p.ṁ, divertor_circuit.valve.p.ṁ], title = "Cycle Flow Rate")
        Wnet = sum(s -> sol[s], [powercycle.turbine.w.Ẇ + powercycle.comp1.w.Ẇ, powercycle.comp2.w.Ẇ, powercycle.comp3.w.Ẇ])
        Wnet = sum(s -> sol[s], [powercycle.turbine.w.Ẇ + powercycle.comp1.w.Ẇ, powercycle.comp2.w.Ẇ, powercycle.comp3.w.Ẇ])
        # # [powercycle,wall_circuit,breeder_circuit,divertor_circuit,intloop]
        p1 = plot(sol, vars = [r.valve.p.ṁ for r in refarr], title = "Mass Flow Rate")
        p2 = plot(sol, vars = [r.HeatIn.Q̇ for r in fusarr], title = "Heat Load (MW)")
        p3 = plot(sol, vars = [r.HeatTx.Q̇ for r in hxarr], title = "Heat Exchanger Rate")
        p4 = plot(sol, vars = [-powercycle.turbine.w.Ẇ, powercycle.comp1.w.Ẇ, powercycle.comp2.w.Ẇ, powercycle.comp3.w.Ẇ], title = "Cycle work")
        p=plot(p1,p2,p3,p4,layout = (2,2))
        display(p)
        return SimplePlant,Plant,copyplant,sol
        # Gas.showsol([(r.s.name for s in r.systems) for r in fusarr],sol)
    end
    
    function test_full_system_rankine()
        @named breeder_circuit  = breeder_loop(Tmin = 700);
        @named divertor_circuit = He_loop(Tmin = 300);
        @named wall_circuit     = He_loop(Tmin = 500);
        @named intloop          = He_inter_loop(Tmin = 250);
        # @named powercycle       = FeedwaterRankine2()
        @named powercycle = FwRankine3()
    
        # Setting random f(t) to test dynamics
        massfcn(t) = 140 + 15 * sin(t)
        Qloads = [700e6 + 100e6*cos(2t), 
                    150e6 + 10e6*cos(2t), 
                    140e6 + 30e6 * sin(.5*t)]
    
        auxeq = Equation[]
    
        for (i,sys) in enumerate([divertor_circuit, wall_circuit, breeder_circuit])
            push!(auxeq,sys.valve.p.ṁ ~ 2*i*massfcn(t))
            push!(auxeq,sys.HeatIn.Q̇ ~ Qloads[i])
        end
    
        for (i,sys) in enumerate([powercycle.valve, intloop.valve])
            push!(auxeq,sys.p.ṁ ~ massfcn(t))
        end
    
        # Adding heat exchanger relations
        @named breeder_hx   = Gen_HeatExchanger(A = breeder_circuit.HeatTx, B = intloop.HeatInC, returnmode = :eq)
        @named wall_hx      = Gen_HeatExchanger(A = wall_circuit.HeatTx, B = intloop.HeatInB, returnmode =:eq)
        @named divertor_hx  = Gen_HeatExchanger(A = divertor_circuit.HeatTx, B = intloop.HeatInA,returnmode =:eq)
        @named primary_hx   = Gen_HeatExchanger(A = powercycle.boiler, B = intloop.HeatTx, returnmode =:eq)

        for eq in [breeder_hx, wall_hx, divertor_hx, primary_hx]
            push!(auxeq,eq...)
        end


        @named Plant = ODESystem(auxeq,t; name = :Plant, systems = [powercycle,wall_circuit,breeder_circuit,divertor_circuit,intloop])
    
        copyplant = deepcopy(Plant)
        

    
        for pd in [Gas.propDict,Steam.hydro_prop_dict,Liq.propDict]
            Plant = substitute(Plant, pd)
        end
    
        # Plant = alias_elimination(Plant)
        SimplePlant = structural_simplify(Plant)
    
        tspan = (0.0,10.0)
    
        saveat = LinRange(1.0,9.0,25)
    
        kwargs = (abstol=1e-10, reltol=1e-2, saveat = saveat)
    
        prob = ODAEProblem(SimplePlant,[],tspan,)
    
        sol   = solve(prob, Rodas4(); kwargs...);

        # refarr  = [powercycle,wall_circuit,breeder_circuit,divertor_circuit,intloop];
        # fusarr  = [breeder_circuit,divertor_circuit,wall_circuit]
        # hxarr   = [breeder_circuit,divertor_circuit,wall_circuit, intloop]
    
        # # [powercycle,wall_circuit,breeder_circuit,divertor_circuit,intloop]
        # p1 = plot(sol, vars = [r.valve.p.ṁ for r in refarr], title = "Mass Flow Rate")
        # p2 = plot(sol, vars = [r.HeatIn.Q̇ for r in fusarr], title = "Heat Load (MW)")
        # p3 = plot(sol, vars = [powercycle.boil.Q̇], title = "Boiler Qin")
        # p4 = plot(sol, vars = [intloop.WorkRes.Ẇ, breeder_circuit.WorkRes.Ẇ, divertor_circuit.WorkRes.Ẇ, powercycle.WorkRes.Ẇ], title = "Circulator work ")
        # p=plot!(p1,p2,p3,p4,layout = (4,1), size = (800,800))
        return copyplant,Plant, sol
    end
end