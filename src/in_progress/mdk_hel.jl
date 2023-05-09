
using ModelingToolkit, Plots, DifferentialEquations, Revise, Unitful, CoolProp, Symbolics, Logging
using NonlinearSolve
Logging.disable_logging(Logging.Warn)
@variables t

begin
    hpt_he(x,y) = PropsSI("HMASS","P",x*u"bar","T",y*u"K","He").val;
    tph_he(x,y) = PropsSI("T","P",x*u"bar","h",y*u"J/kg","He").val;
end

#   PROPERTIES DICT
begin
    cphe(x) = 5192.0    #J/kg/K
    cvhe(x) = 3115.8    #J/kg/K
    khe(x) = 1.667
    @variables x
    @register_symbolic cpfunc(x)
    @register_symbolic cvfunc(x)
    @register_symbolic kfunc(x)
    propDict = Dict(cpfunc => cphe, kfunc => khe)
end

# PINS AND CONNECTORS
begin
    # display(methods(cpfunc))
    @connector function ThermoPin(; name, Pdef = 10.0, Tdef = 300, ṁdef = 0.0)
        @variables P(t) = Pdef              #[unit = u"bar"] 
        @variables ṁ(t) = ṁdef              #[unit = u"kg/s"] 
        @variables T(t)  = Tdef             #[unit = u"K"]
        @variables Φ(t) = 1.0
        @variables cp(t) = 5192
        @variables k(t)  = 1.667
        eq = [
            cp ~ cpfunc(T)
            k ~ kfunc(T)
        ]
        sts = [T, P, ṁ, cp, k, Φ]
        ODESystem(eq, t, sts, []; name = name, defaults = [P => Pdef, T => Tdef, ṁ => ṁdef, Φ => 0.0])
    end

    @connector function WorkPin(; name)
        sts = @variables Ẇ(t) = 0.0
        ODESystem(Equation[], t, sts, []; name = name)
    end

    # convention = +Q = Qin
    @connector function FixedHeatFlowPin(; name, Qin = 1e6)  #input in W
        sts = @variables Q̇(t)=Qin 
        ps = @parameters Qin = Qin
        eqs=[Q̇ ~ Qin]
        ODESystem(eqs,t, sts, ps; name = name, defaults = [Q̇ => Qin])
    end

    @connector function HeatTransferPin(; name)  #input in W
        @variables Q̇(t) =0.0
        ODESystem(Equation[],t, [Q̇], []; name = name, defaults = [Q̇ => 0.0])
    end
end

## Connector functions
begin
    """
        hydro_connect(pins...)
        Adds mass and energy balance eqs to all nodes

    """
    function gas_connect(pins...)

        eqs = [
            sum(pin -> pin.ṁ, pins) ~ 0.0, # Mass
            sum(pin -> pin.Φ, pins) ~ 0.0, # Energy
        ]

        for pin in pins[2:end]
            push!(eqs, pins[1].P ~ pin.P)   # p_1 = p_2, p_1 = p_3,...
            push!(eqs, pins[1].T ~ pin.T)   # p_1 = p_2, p_1 = p_3,...
        end
        return eqs
    end
    """
        gas_mirror_pin(pinA,pinB) 
        Used for internal pins to switch mass flow rate and energy flow rate
    """
    function gas_CopyConnect(pins...)   # creates thermoPin connections  P = P, T = T, mflow = mflow
        eqs = Equation[]

        for pin in pins[2:end]
            push!(eqs, pins[1].P ~ pin.P)   # p_1 = p_2, p_1 = p_3,...
            push!(eqs, pins[1].T ~ pin.T)   # p_1 = p_2, p_1 = p_3,...
            push!(eqs, pins[1].ṁ ~ pin.ṁ)   # p_1 = p_2, p_1 = p_3,...
            push!(eqs, pins[1].Φ ~ pin.Φ)   # p_1 = p_2, p_1 = p_3,...
        end

        return eqs
    end


    """
        Dont use this for flowsource
    """
    function series_connect(pins)
        # @assert length(pins)>2
        eqs =  basic_connect(pins[1].n,pins[2].p)   # n -> p

        if length(pins)>2
            for i = 2:length(pins)-1
                eqs = vcat(eqs,  basic_connect(pins[i].n,pins[i+1].p))
            end
        end
        return eqs
    end

    function hx_connect(hx, compAin, compAout, compBin, compBout)

        eqs = vcat(gas_connect(compAin.n,hx.A.p),
                    gas_connect(hx.A.n,compAout.p),
                    gas_connect(compBin.n,hx.B.p),
                    gas_connect(hx.B.n,compBout.p))

        return eqs
    end
end

## RESERVOIRs - TEMPERATURE SETTERS
begin
    @component function SetPressure(;name, P = 0.1)
        @named p = ThermoPin(Pdef = P)
        @named n = ThermoPin(Pdef = P)
        ps = @parameters P = P 
        eqs = [
            p.P ~ P
            n.P ~ p.P
            n.T ~ p.T
            0 ~ p.Φ + n.Φ             # conservation of energy
            0 ~ p.ṁ + n.ṁ
        ]
        ODESystem(eqs, t,[], ps; name = name, systems = [p,n], defaults = [P => 0.1])
    end

    @component function SetTemperature(;name, T = 300)
        @named p = ThermoPin(Tdef = T)
        @named n = ThermoPin(Tdef = T)
        ps = @parameters T = T 
        eqs = [
            p.T ~ T
            n.T ~ T
            p.P ~ n.P
            0 ~ p.Φ + n.Φ             # conservation of energy
            0 ~ p.ṁ + n.ṁ
        ]
        ODESystem(eqs, t,[], ps; name = name, systems = [p,n], defaults = [T => 300])
    end
end

#   HELIUM COMPONENTS
begin
    @component function SinglePortReservoir(;name, P = 0.1, T = 300)
        @named n = ThermoPin(Pdef = P, Tdef = T)
        ps = @parameters P=P 
        eqs = [
            n.P ~ P
            n.T ~ T
            n.Φ ~ 0
        ]
        ODESystem(eqs, t,[], ps; name = name, systems = [n], defaults = [P => 0.1, T =>300])
    end

    @component function TwoPortReservoir(;name, P = 0.1, T = 300)
        @named p = ThermoPin(Pdef = P, Tdef = T)
        @named n = ThermoPin(Pdef = P, Tdef = T)
        ps = @parameters P=P T=T
        sts = @variables ΔΦ(t) = 0.0
        eqs = [
            0 ~ p.ṁ + n.ṁ 
            n.P ~ P
            n.T ~ T
            p.T ~ n.T
            p.P ~ n.P
            n.Φ ~ 0
            ΔΦ ~ p.Φ - n.Φ
        ]
        ODESystem(eqs, t,sts, ps; name = name, systems = [p,n], defaults = [P => 0.1, T =>300, ΔΦ => 0.0])
    end

    @component function GasFlowSource(;name, Ṁ = 1.0)
        @named n = ThermoPin(ṁdef = -Ṁ)
        @named p = ThermoPin(ṁdef = Ṁ)

        ps = @parameters Ṁ = Ṁ
        eqs = [
                0 ~ p.Φ + n.Φ   # conservation of energy
                0 ~ p.ṁ + n.ṁ
                p.ṁ ~ Ṁ
                p.T ~ n.T 
                n.P ~ p.P
            ]
        ODESystem(eqs, t, [], ps; name = name, systems = [n,p], defaults = [Ṁ => 1.0])
    end

    @component function GasFlowValve(;name)
        @named n = ThermoPin()
        @named p = ThermoPin()

        ps = @variables Ṁ(t)
        eqs = [
                0 ~ p.Φ + n.Φ   # conservation of energy
                0 ~ p.ṁ + n.ṁ
                p.ṁ ~ Ṁ
                p.T ~ n.T 
                n.P ~ p.P
            ]
        ODESystem(eqs, t, ps, []; name = name, systems = [n,p], defaults = [Ṁ => 1.0])
    end
    
    
    @component function ThermoHeatSource(; name, Qin = 1e6)
        @named p = ThermoPin()
        @named n = ThermoPin()
        @named q = FixedHeatFlowPin(Qin = Qin)
        #   0 variables 3 equations for pin

        eqs = [
            q.Q̇ ~ p.Φ + n.Φ             # conservation of energy
            0  ~ p.ṁ + n.ṁ  # p.ṁ ~ n.ṁ                               # conservation of mass
            n.T ~ p.T + q.Q̇/(p.ṁ*p.cp)
            n.P ~ p.P
        ]
        ODESystem(eqs,t,[],[]; name = name, systems = [p,n,q])
    end
    
    @component function ThermoHeatTransfer(; name)
        @named p = ThermoPin()
        @named n = ThermoPin()
        @named q = HeatTransferPin()
        st = @variables Q̇(t)=0.0 C(t)=5192
        eqs = [
            0 ~ p.ṁ + n.ṁ           # p.ṁ ~ n.ṁ                               # conservation of mass
            q.Q̇ ~ p.Φ + n.Φ             # conservation of energy            
            C ~ p.ṁ * p.cp          # duty
            0 ~ q.Q̇ - Q̇ 
            n.T ~ p.T + q.Q̇/(p.ṁ*p.cp)
            n.P ~ p.P
        ]
        ODESystem(eqs,t,[Q̇,C],[]; name = name, systems = [p,n,q], defaults = [Q̇ => 0.0, C => 5192])
    end

    @component function FixedThermoHeatTransfer(; name, Qin = 1e6)
        @named p = ThermoPin()
        @named n = ThermoPin()
        @named q = HeatTransferPin()
        st = @variables C(t)=5192
        ps = @parameters Q̇=Qin
        eqs = [
            0 ~ p.ṁ + n.ṁ           # p.ṁ ~ n.ṁ                               # conservation of mass
            q.Q̇ ~ p.Φ + n.Φ             # conservation of energy            
            C ~ p.ṁ * p.cp          # duty
            0 ~ q.Q̇ - Q̇ 
            n.T ~ p.T + Q̇/C
            n.P ~ p.P
        ]
        ODESystem(eqs,t,[C],[Q̇]; name = name, systems = [p,n,q], defaults = [C => 5192])
    end
    
    @component function ActiveThermoCompressor(; name, η = 1.0, rp = 3.5)
        @named p = ThermoPin()
        @named n = ThermoPin()
        @named w = WorkPin()
    
        ps = @parameters rp = rp η = η
        eqs = [
            0 ~ p.ṁ + n.ṁ #p.ṁ ~ n.ṁ                               # conservation of mass
            n.P ~ p.P * rp
            n.T ~ p.T * ((1.0 - η - (rp^((p.k-1)/p.k))) / (-η))
            w.Ẇ ~ p.ṁ * p.cp * (n.T - p.T)
            w.Ẇ ~ p.Φ + n.Φ
        ]
        ODESystem(eqs,t,[],ps; name = name, systems = [p,n,w], defaults = [η => η, rp => rp])
    end

    @component function PassiveThermoCompressor(; name, η = 1.0)
        @named p = ThermoPin()
        @named n = ThermoPin()
        @named w = WorkPin()
    
        ps = @parameters η = η
        sts = @variables rp(t) = 3.5

        eqs = [
            0 ~ p.ṁ + n.ṁ #p.ṁ ~ n.ṁ                               # conservation of mass
            n.P ~ p.P * rp
            n.T ~ p.T * ((1.0 - η - (rp^((p.k-1)/p.k))) / (-η))
            w.Ẇ ~ p.ṁ * p.cp * (n.T - p.T)
            w.Ẇ ~ p.Φ + n.Φ
        ]
        ODESystem(eqs,t,sts,ps; name = name, systems = [p,n,w], defaults = [rp => 3.5])
    end
    
    @component function ActiveThermoTurbine(; name, η = 1.0, rp = 3.5)
        @named p = ThermoPin(Pdef = 80, Tdef = 800)
        @named n = ThermoPin(Pdef = 80/3.5, Tdef = 500)
        @named w = WorkPin()
        ps = @parameters rp = rp η = η
        eqs = [
            0 ~ p.ṁ + n.ṁ                               # conservation of mass
            p.P ~ n.P * rp
            n.T ~ p.T *(η + rp^((p.k-1)/p.k)- η*(rp^((p.k-1)/p.k))) / (rp^((p.k-1)/p.k))
            w.Ẇ ~ p.ṁ * p.cp * (n.T - p.T)
            w.Ẇ ~ p.Φ + n.Φ
        ]
        ODESystem(eqs,t,[],ps; name = name, systems =  [p,n,w])
    end

    @component function PassiveThermoTurbine(; name, η = 1.0)
        @named p = ThermoPin(Pdef = 80, Tdef = 800)
        @named n = ThermoPin(Pdef = 80/3.5, Tdef = 500)
        @named w = WorkPin()
        ps = @parameters η = η
        sts = @variables rp(t) = 3.5
        eqs = [
            0 ~ p.ṁ + n.ṁ                               # conservation of mass
            p.P ~ n.P * rp
            n.T ~ p.T *(η + rp^((p.k-1)/p.k)- η*(rp^((p.k-1)/p.k))) / (rp^((p.k-1)/p.k))
            w.Ẇ ~ p.ṁ * p.cp * (n.T - p.T)
            w.Ẇ ~ p.Φ + n.Φ
        ]
        ODESystem(eqs,t,sts,ps; name = name, systems =  [p,n,w],  defaults = [rp => 3.5])
    end
    
    @component function IdealCooler(;name)
        @named p = ThermoPin()
        @named n = ThermoPin()
        @named q = HeatTransferPin()
        sts = @variables Q̇(t) = 0.0
        # 0 variables, 3 equations for pin
        eqs = [
            q.Q̇ ~ p.Φ + n.Φ             # conservation of energy
            0 ~ p.ṁ + n.ṁ
            n.P  ~  p.P                 # no pressure   
            q.Q̇ ~ p.ṁ * p.cp * (n.T - p.T)
            Q̇ ~ q.Q̇
            ]
        ODESystem(eqs, t, sts, []; name = name, systems = [p,n,q], defaults = [Q̇ => 0])
    end

    @component function PassiveElement(;name)
        @named p = ThermoPin()
        @named n = ThermoPin()
        sts = @variables Q̇(t) = 0.0 ΔP(t) ΔT(t) Δṁ(t)
        # 0 variables, 3 equations for pin
        eqs = [
            n.ṁ + p.ṁ ~ Δṁ    #has to be negative
            ΔT ~ p.T - n.T
            Q̇ ~ -(p.ṁ * p.cp * ΔT)
            ΔP ~ p.P - n.P
            Q̇ ~ p.Φ + n.Φ             # conservation of energy
            ]
        ODESystem(eqs, t, [Q̇,ΔP,ΔT,Δṁ], []; name = name, systems = [p,n], defaults = [Q̇ => 0, ΔT =>0, ΔP => 0, Δṁ => 0])
    end
    
    @component function Regenerator(;name, ϵ = 0.95)
        @named A = ThermoHeatTransfer()
        @named B = ThermoHeatTransfer()
        ps = @parameters ϵ = ϵ
        @variables Q̇(t)=0.0 Cmin(t)=0.0
        eqs = [
            Q̇ ~ ϵ * (A.p.T - B.p.T) * ((A.C < B.C) *  A.C + (A.C >= B.C) * B.C)      #   heat transfer out of A -> B , if A/T > B/T 
            0 ~ A.Q̇ + B.Q̇
            A.Q̇ ~ -Q̇
        ]
        ODESystem(eqs, t, [Q̇], ps; name = name, systems = [A,B], defaults = [ϵ => 0.95] )
    end
end

## compound components
begin
    @component function  Intercooler(;name, Tout = 300)
        @named p = ThermoPin()
        @named n = ThermoPin()

        # have to switch signs for internal circiuts
        sts = @variables Q̇(t) = 0.0

        @named st = SetTemperature(T = Tout) 
        @named ic = IdealCooler()

        eqs = vcat(gas_CopyConnect(p,ic.p),
                    gas_connect(ic.n,st.p),
                    gas_CopyConnect(st.n,n),
                    Q̇ ~ ic.Q̇)
        ODESystem(eqs, t, [Q̇], []; name = name, systems = [p,n,st,ic], defaults = [Q̇ =>0.0] )
    end
end

@named ICc = Intercooler()
@named reg = Regenerator()
##


##############################################################################
using Printf
function showsol(c,sol)
    for cel in c
        cvec = []
        if hasproperty(cel,:p)
            push!(cvec,cel.p)
        end
        if hasproperty(cel,:n)
            push!(cvec,cel.n)
        end
        for n in cvec

            @printf "\t %s ṁ = %.2f \t P = % .2f \t T = %.2f " n.name sol[n.ṁ][end] sol[n.P][end] sol[n.T][end]

            if hasproperty(cel,:Q̇)
                @printf " Q̇ = %.2e " sol[cel.Q̇][end] 
            elseif hasproperty(cel,:q)
                @printf " Q̇ = %.2e " sol[cel.q.Q̇][end] 
            end
            if hasproperty(cel,:C)
                @printf " C = %.2e " sol[cel.C][end] 
            end
            if hasproperty(cel,:w)
                @printf " Ẇ = %.2e " sol[cel.w.Ẇ][end] 
            end
            @printf "\n"
        end
    end

end
# testing Regenerator
begin
    @named hotres = SinglePortReservoir(P=100, T=500)
    @named hmsrc   = GasFlowSource(Ṁ = 1.0)

    @named coolres = SinglePortReservoir(P=100, T=400)
    @named cmsrc   = GasFlowSource(Ṁ = 2.0)

    @named reg    = Regenerator()
    connections = vcat(gas_connect(hotres.n,hmsrc.p),
            gas_connect(hmsrc.n,reg.A.p),
            gas_connect(coolres.n,cmsrc.p),
            gas_connect(cmsrc.n,reg.B.p))

    systemNames =[hotres,hmsrc,coolres,cmsrc,reg]

    @named odemodel = ODESystem(connections,t; systems = systemNames)
    
    smodel = substitute(odemodel,propDict)
    simple_mod = structural_simplify(smodel)
    prob = ODAEProblem(simple_mod,Pair[],(0.0,1.0))
    sol=solve(prob)
    showsol([hotres,hmsrc,coolres,cmsrc,reg.A,reg.B],sol)
end

##
# testing intercooler
begin
    @named hotres = SinglePortReservoir(P=100, T=400)
    @named hmsrc  = GasFlowSource(Ṁ = 1.0)
    @named IC       = Intercooler(Tout = 320)

    connections = vcat(gas_connect(hotres.n,hmsrc.p),
                gas_connect(hmsrc.n,IC.p))

    systemNames =[hotres,hmsrc,IC]

    @named odemodel = ODESystem(connections,t; systems = systemNames)
    
    smodel = substitute(odemodel,propDict)
    simple_mod = structural_simplify(smodel)
    prob = ODAEProblem(simple_mod,Pair[],(0.0,1.0))
    sol=solve(prob)
    showsol(systemNames,sol)
end

# TURBINE
begin
    @named hotres   = SinglePortReservoir(P=100, T=1000)
    @named hmsrc    = GasFlowSource(Ṁ = 1.0)
    @named turbine  = PassiveThermoTurbine(η = 0.9)
    @named ep       = SetPressure(P=33)

    connections = vcat(gas_connect(hotres.n,hmsrc.p),
                    gas_connect(hmsrc.n,turbine.p), gas_connect(turbine.n,ep.p))

    systemNames =[hotres,hmsrc,turbine,ep]

    @named odemodel = ODESystem(connections,t; systems = systemNames)
    
    smodel = substitute(odemodel,propDict)

    simple_mod = structural_simplify(smodel)

    prob = ODAEProblem(simple_mod,Pair[],(0.0,1.0))
    sol=solve(prob)
    showsol(systemNames,sol)
end

# HEATER
begin
    @named cres = SinglePortReservoir(P=100, T=10)
    @named msrc  = GasFlowSource(Ṁ = 75)
    @named heat  = FixedThermoHeatTransfer(Qin = 150e6)
    @named ep       = SetPressure(P=100)

    connections = vcat(gas_connect(cres.n,msrc.p),
                    gas_connect(msrc.n,heat.p),
                    gas_connect(heat.n,ep.p))

    systemNames =[cres,msrc,heat]

    @named odemodel = ODESystem(connections,t; systems = systemNames)
    
    smodel = substitute(odemodel,propDict)
    simple_mod = structural_simplify(smodel)
    prob = ODAEProblem(simple_mod,Pair[],(0.0,1.0))
    sol=solve(prob)
    showsol(systemNames,sol)
end

# SIMPLE BRAYTON
begin
    @named cmsrc   = GasFlowSource(Ṁ = 100)
    @named comp    = ThermoCompressor()
    @named heat    = ThermoHeatSource(Qin = 100e6)
    @named turbine = PassiveThermoTurbine()
    @named cool    = IdealCooler()
    @named res     = TwoPortReservoir(P = 10, T = 300)

    connections = vcat(gas_connect(res.n,cmsrc.p),
            gas_connect(cmsrc.n,comp.p),
            gas_connect(comp.n,heat.p),
            gas_connect(heat.n,turbine.p),
            gas_connect(turbine.n,cool.p),
            gas_connect(cool.n,res.p))

    # systemNames =[coolres,cmsrc,comp,heat,turbine,cool,outres]
    systemNames = [res,cmsrc,comp,heat,turbine,cool]

    @named odemodel = ODESystem(connections,t; systems = systemNames)
    
    smodel = substitute(odemodel,propDict)
    simple_mod = structural_simplify(smodel)
    prob = ODAEProblem(simple_mod,Pair[],(0.0,1.0))
    sol=solve(prob)
    @show Qnet = sol[cool.q.Q̇] + sol[heat.q.Q̇]
    @show Wnet = sol[turbine.w.Ẇ] + sol[comp.w.Ẇ]
    @show Wnet + Qnet
    @show Wnet/sol[heat.q.Q̇]
end

# Brayton with regen
begin
    TminCycle = 300; 
    PminCycle = 15;
    Qload = 1e6;
    Mflow = 75;

    @named res     = TwoPortReservoir(P = PminCycle, T = TminCycle)
    @named valve   = GasFlowSource(Ṁ = Mflow)
    @named comp1   = ActiveThermoCompressor(rp = 3.8, η = 0.9)
    @named ic1     = Intercooler(Tout = TminCycle)
    @named comp2   = ActiveThermoCompressor(rp = 3.5, η = 0.95)
    @named ic2     = Intercooler(Tout = TminCycle)
    @named comp3   = ActiveThermoCompressor(rp = 3.5, η = 0.95)
    @named reg     = Regenerator()
    @named heat    = FixedThermoHeatTransfer(Qin = 1000e6)
    @named turbine = PassiveThermoTurbine()
    @named cool    = IdealCooler()
    
    connections = vcat(gas_connect(res.n,valve.p),
    gas_connect(valve.n,comp1.p),
    gas_connect(comp1.n,ic1.p),
    gas_connect(ic1.n,comp2.p),
    gas_connect(comp2.n,ic2.p),
    gas_connect(ic2.n,comp3.p),
    hx_connect(reg,comp3,heat,turbine,cool),
    gas_connect(heat.n,turbine.p),
    gas_connect(cool.n,res.p));

    # systemNames =[coolres,valve,comp,heat,turbine,cool,outres]
    systemNames = [res,valve,comp1,ic1,comp2,ic2,comp3,reg,heat,turbine,cool];
    @named odemodel = ODESystem(connections,t; systems = systemNames);
    smodel = substitute(odemodel,propDict);
    simple_sys = structural_simplify(smodel);
end


## dynamic model
begin
    TminCycle = 300; 
    PminCycle = 15;
    Qload = 1e6;
    Mflow = 75;

    @named res     = TwoPortReservoir(P = PminCycle, T = TminCycle)
    @named valve   = GasFlowValve()
    @named comp1   = ActiveThermoCompressor(rp = 1.7, η = 0.9)
    @named ic1     = Intercooler(Tout = TminCycle)
    @named comp2   = ActiveThermoCompressor(rp = 1.5, η = 0.95)
    @named ic2     = Intercooler(Tout = TminCycle)
    @named comp3   = ActiveThermoCompressor(rp = 1.5, η = 0.95)
    @named reg     = Regenerator()
    # @named heat    = FixedThermoHeatTransfer(Qin = 300e6)
    @named heat = ThermoHeatTransfer()
    @named turbine = PassiveThermoTurbine()
    @named cool    = IdealCooler()
    
    connections = vcat(gas_connect(res.n,valve.p),
    gas_connect(valve.n,comp1.p),
    gas_connect(comp1.n,ic1.p),
    gas_connect(ic1.n,comp2.p),
    gas_connect(comp2.n,ic2.p),
    gas_connect(ic2.n,comp3.p),
    hx_connect(reg,comp3,heat,turbine,cool),
    gas_connect(heat.n,turbine.p),
    gas_connect(cool.n,res.p));

    mass_flow_fcn(t) = t < 5 ? 100 : 75
    Q̇in_fcn(t) = 100e6 + 35e6 *sin(t)  

    @register_symbolic mass_flow_fcn(t)

    controlled_eqs = [valve.p.ṁ ~ mass_flow_fcn(t),
                        heat.Q̇ ~ Q̇in_fcn(t)]

    @named aux_sys = ODESystem(controlled_eqs,t)


    # systemNames =[coolres,valve,comp,heat,turbine,cool,outres]
    systemNames = [res,valve,comp1,ic1,comp2,ic2,comp3,reg,heat,turbine,cool];
    @named odemodel = ODESystem(connections,t; systems = systemNames);

    @named odemodel = extend(aux_sys,odemodel)
    smodel = substitute(odemodel,propDict);
    simple_sys = structural_simplify(smodel)
    tspan = (0.0,10.0)

    u0 = [valve.p.ṁ => 100
        turbine.p.T=> 1000
        heat.p.cp => 5192
        turbine.rp=>3.5]
        
    prob = ODEProblem(simple_sys,u0,tspan,)
    sol   = solve(prob, Rodas4(); kwargs...);
    c   = [res,valve,comp1,reg.A,heat,turbine,reg.B,cool];
    p1 = plot(sol,vars = [valve.p.ṁ], title = "Cycle Flow Rate")
    p2 = plot(sol, vars = [ heat.Q̇], title = "Heat Loading")
    p3 = plot(sol,vars = [heat.n.T, reg.B.n.T, reg.B.p.T, heat.p.T,heat.n.T])

    p=plot(p1,p2,p3,layout = (3,1))

    display(p)
    showsol(c,sol)
end


# tspan = (0.0,10.0)


# u0 = [valve.p.ṁ => 50
#         turbine.p.T=> 1000
#         heat.p.cp => 5192
#         turbine.rp=>3.5]

# prob = ODEProblem(simple_sys,u0,tspan,)
# sol   = solve(prob, Rodas4(); kwargs...);
# c   = [res,valve,comp1,reg.A,reg.B,heat,turbine,reg.B,cool];

# p1 = plot(sol,vars = [valve.p.ṁ, heat.Q̇])
# p2 = plot(sol,vars = [heat.n.T, reg.B.n.T, reg.B.p.T, heat.p.T,heat.n.T])

# showsol(c,sol)

# # p=[res.P => 30
#     # res.T => 300
#     # comp1.rp =>1.25
#     # comp1.η => 1.0
#     # ic1.st.T => 300
#     # comp2.rp => 1.5
#     # comp2.η => 1.0
#     # heat.Q̇ => 800e6
#     # turbine.η => 0.9
#     # reg.ϵ   => 0.9];


# # u0=[    aux_sys₊valve.p.
# #     turbine₊p₊T(t) [defaults to 800]
# #     heat₊p₊cp(t) [defaults to 5192]
# #     turbine₊rp(t)]
# kwargs = (abstol=1e-10, reltol=1e-10)
# prob = ODEProblem(simple_sys,u0,tspan,)
# sol1e   = solve(prob, Rodas4(); kwargs...);
# showsol(c,sol1e)


# sol = solve(prob)
# plot(sol,vars = [turbine.p.ṁ])
# # probae = ODAEProblem(simple_mod, u0,tspan,p)

# begin
#     simple_mod = structural_simplify(smodel);

#     tspan = (0.0,10.0);

#     # u0 = [reg.A.n.T => 500 
#     #     reg.A.p.T => 600
#     #     reg.B.n.T => 550
#     #     cmsrc.p.ṁ => 100
#     #     heat.p.T => 550
#     #     heat.n.T => 1000
#     #     turbine.p.T => 1000];

#     # p=[res.P => 30
#     #     res.T => 300
#     #     cmsrc.Ṁ => 10
#     #     comp1.rp => 3.5
#     #     comp.η => 1.0
#     #     heat.q.Qin => 10e6
#     #     turbine.η => 0.9
#     #     reg.ϵ   => 0.9];


#     # prob = ODEProblem(simple_mod, u0,tspan)
#     # parameters(simple_mod)


#     probde = ODEProblem(simple_mod,u0,tspan,p)
#     probae = ODAEProblem(simple_mod, u0,tspan,p)

#     sol1=solve(probae, Rodas4())
#     sol2=solve(probae,  Tsit5(), abstol = 1e-8, reltol = 1e-8)
#     sol3=solve(probn, Rodas4());
#     sol4=solve(probn, Rodas4());

#     sol1de=solve(probde, Rodas4());
# end


#     # pf(xin) = prob.f.observed(xin,prob.u0,prob.p,1.0).

# using Printf
# # nodes = [res.p, cmsrc.p, comp.p, reg.A.p, reg.A.n, heat.p, turbine.p, reg.B.p, reg.B.n];
# c   = [res,cmsrc,comp,reg.A,reg.B,heat,turbine,reg.B,cool];
# ht  = [reg.A,reg.B,cool,heat];

# sol=sol2;

# sol[heat.n.T]

# begin
#     for cel in c
#         for n in [cel.p,cel.n]
#             @printf "\t %s ṁ = %.2f \t P = % .2f \t T = %.2f\n" n.name sol[n.ṁ][end] sol[n.P][end] sol[n.T][end]
#         end
#     end

#     for hel in ht
#         n = hel.q
#         @printf "\t %s Q̇ = %.2f ΔT = %.2f \n" n.name sol[n.Q̇][end] (sol[hel.n.T][end]-sol[hel.p.T][end])
#     end
# end









# nodes = [reg.A.n, reg.B.n, heat.p, heat.n]
# using Printf
# for n in nodes
#     @printf "\t %s ṁ = %.2f \t P = % .2f \t T = %.2f\n" n.name pf(n.ṁ) pf.(n.P) pf(n.T)
# end


# @named traced_sys = modelingtoolkitize(prob)

# mtksys = structural_simplify(traced_sys)


# states(traced_sys)
# observed(mtksys)

# mtkprob = ODEProblem(mtksys,[],tspan)
    
# # odeprob = ODEProblem(mtksys)
# sol=solve(mtkprob, Tsit5(), abstol = 1e-8,reltol=1e-8)
# # # @show Qnet = sol[cool.q.Q̇] + sol[heat.q.Q̇]
# # @show Wnet = sol[turbine.w.Ẇ] + sol[comp.w.Ẇ]
# # @show Wnet + Qnet
# # @show Wnet/sol[heat.q.Q̇]
# # end