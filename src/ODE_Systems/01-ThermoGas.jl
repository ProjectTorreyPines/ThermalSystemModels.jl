# module ThermoGas
#=============================================================#
#                      THERMOGAS
# DEFAULT FLUID = HELIUM
#   units
#       pressure, bar
#       temperature, Kelvin
#       mass flow, kg/s
#
#   Sign convention
#       Heat in (Q) is positive
#       Heat out is negative
#       Work out is negative
#       Work in is positive
#=============================================================#
using ModelingToolkit, Plots, DifferentialEquations, Revise,  Logging, Printf
# using NonlinearSolve, Printf
# Unitful, CoolProp,
Revise.retry()
# Logging.disable_logging(Logging.Warn)
@variables t
# include("03-MTK_UTILS.jl")
# begin
#     hpt_he(x,y) = PropsSI("HMASS","P",x*u"bar","T",y*u"K","He").val;
#     tph_he(x,y) = PropsSI("T","P",x*u"bar","h",y*u"J/kg","He").val;
# end
#   PROPERTIES DICT
cphe(x) = 5192.0    #J/kg/K
cvhe(x) = 3115.8    #J/kg/K
khe(x) = 1.667

gcpfunc(x) = cphe(x)
gcvfunc(x) = cvhe(x)
gkfunc(x) = khe(x)

@variables x
@register_symbolic gcpfunc(x)
@register_symbolic gcvfunc(x)
@register_symbolic gkfunc(x)
propDict = Dict(gcpfunc => cphe, gkfunc => khe)

# PINS AND CONNECTORS
# display(methods(cpfunc))
function ThermoPin(; name, Pdef = 10.0, Tdef = 300, ṁdef = 0.0)
    @variables P(t) = Pdef              #[unit = u"bar"] 
    @variables ṁ(t) = ṁdef              #[unit = u"kg/s"] 
    @variables T(t) = Tdef             #[unit = u"K"] 
    @variables Φ(t) = 1.0
    @variables cp(t) = 5192
    @variables k(t) = 1.667
    eq = [
        cp ~ gcpfunc(T)
        k ~ gkfunc(T)
    ]
    sts = [T, P, ṁ, cp, k, Φ]
    ODESystem(
        eq,
        t,
        sts,
        [];
        name = name,
        defaults = [P => Pdef, T => Tdef, ṁ => ṁdef, Φ => 0.0],
    )
end

function WorkPin(; name)
    sts = @variables Ẇ(t) = 0.0
    ODESystem(Equation[], t, sts, []; name = name)
end

# convention = +Q = Qin
function FixedHeatFlowPin(; name, Qin = 1e6)  #input in W
    sts = @variables Q̇(t) = Qin
    ps = @parameters Qin = Qin
    eqs = [Q̇ ~ Qin]
    ODESystem(eqs, t, sts, ps; name = name, defaults = [Q̇ => Qin])
end

function HeatTransferPin(; name)  #input in W
    @variables Q̇(t) = 0.0
    ODESystem(Equation[], t, [Q̇], []; name = name, defaults = [Q̇ => 0.0])
end

## Connector functions
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
    eqs = basic_connect(pins[1].n, pins[2].p)   # n -> p

    if length(pins) > 2
        for i = 2:length(pins)-1
            eqs = vcat(eqs, basic_connect(pins[i].n, pins[i+1].p))
        end
    end
    return eqs
end

function hx_connect(hx, compAin, compAout, compBin, compBout)

    eqs = vcat(
        gas_connect(compAin.n, hx.A.p),
        gas_connect(hx.A.n, compAout.p),
        gas_connect(compBin.n, hx.B.p),
        gas_connect(hx.B.n, compBout.p),
    )

    return eqs
end

## RESERVOIRs - TEMPERATURE SETTERS
"""
    function SetPressure(; name, P = 0.1)

"""
function SetPressure(; name, P = 0.1)
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
    ODESystem(eqs, t, [], ps; name = name, systems = [p, n])
end

function SetTemperature(; name, T = 300)
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
    ODESystem(eqs, t, [], ps; name = name, systems = [p, n], defaults = [T => 300])
end

#   HELIUM COMPONENTS
function SinglePortReservoir(; name, P = 0.1, T = 300)
    @named n = ThermoPin(Pdef = P, Tdef = T)
    ps = @parameters P = P T = T
    eqs = [
        n.P ~ P
        n.T ~ T
        n.Φ ~ 0
    ]
    ODESystem(eqs, t, [], ps; name = name, systems = [n])
end

function TwoPortReservoir(; name, P = 0.1, T = 300)
    @named p = ThermoPin(Pdef = P, Tdef = T)
    @named n = ThermoPin(Pdef = P, Tdef = T)
    ps = @parameters P = P T = T
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
    ODESystem(eqs, t, sts, ps; name = name, systems = [p, n], defaults = [ΔΦ => 0.0])
end

function GasFlowSource(; name, Ṁ = 1.0)
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
    ODESystem(eqs, t, [], ps; name = name, systems = [n, p])
end

function GasFlowValve(; name)
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
    ODESystem(eqs, t, ps, []; name = name, systems = [n, p], defaults = [Ṁ => 1.0])
end

function ThermoHeatSource(; name, Qin = 1e6)
    @named p = ThermoPin()
    @named n = ThermoPin()
    @named q = FixedHeatFlowPin(Qin = Qin)
    #   0 variables 3 equations for pin

    eqs = [
        q.Q̇ ~ p.Φ + n.Φ             # conservation of energy
        0 ~ p.ṁ + n.ṁ  # p.ṁ ~ n.ṁ                               # conservation of mass
        n.T ~ p.T + q.Q̇ / (p.ṁ * p.cp)
        n.P ~ p.P
    ]
    ODESystem(eqs, t, [], []; name = name, systems = [p, n, q])
end

function ThermoHeatTransfer(; name, ΔP = 0.0)
    @named p = ThermoPin()
    @named n = ThermoPin()
    @named q = HeatTransferPin()
    st = @variables Q̇(t) = 0.0 C(t) = 5192
    ps = @parameters ΔP = ΔP
    eqs = [
        0 ~ p.ṁ + n.ṁ           # p.ṁ ~ n.ṁ                               # conservation of mass
        q.Q̇ ~ p.Φ + n.Φ         # conservation of energy            
        C ~ p.ṁ * p.cp          # duty
        0 ~ q.Q̇ - Q̇
        n.T ~ p.T + q.Q̇ / (p.ṁ * p.cp)
        n.P ~ p.P - ΔP
    ]
    ODESystem(
        eqs,
        t,
        [Q̇, C],
        ps;
        name = name,
        systems = [p, n, q],
        defaults = [Q̇ => 0.0, C => 5192],
    )
end

"""
    FlowControlThermoHeatTransfer(; name, ΔP = 0.0, Tout)
    Ability to change mass flow rate to achieve desired outlet temperature
"""
function FlowControlThermoHeatTransfer(; name, ΔP = 0.0, Tout = 773.0)
    @named p = ThermoPin()
    @named n = ThermoPin()
    @named q = HeatTransferPin()
    st = @variables Q̇(t) = 0.0 C(t) = 5192 ṁ(t) = 1.0
    ps = @parameters ΔP = ΔP Tout = Tout
    eqs = [
        n.T ~ Tout
        p.ṁ ~ q.Q̇ / ((n.T - p.T) * p.cp)           # q = mcp*ΔT
        0 ~ p.ṁ + n.ṁ           # p.ṁ ~ n.ṁ                               # conservation of mass
        q.Q̇ ~ p.Φ + n.Φ         # conservation of energy            
        C ~ p.ṁ * p.cp          # duty
        0 ~ q.Q̇ - Q̇
        n.P ~ p.P - ΔP
    ]
    ODESystem(
        eqs,
        t,
        [Q̇, C, ṁ],
        ps;
        name = name,
        systems = [p, n, q],
        defaults = [Q̇ => 0.0, C => 5192.0, ṁ => 1.0],
    )
end

function FixedThermoHeatTransfer(; name, Qin = 1e6)
    @named p = ThermoPin()
    @named n = ThermoPin()
    @named q = HeatTransferPin()
    st = @variables C(t) = 5192
    ps = @parameters Q̇ = Qin
    eqs = [
        0 ~ p.ṁ + n.ṁ           # p.ṁ ~ n.ṁ                               # conservation of mass
        q.Q̇ ~ p.Φ + n.Φ             # conservation of energy            
        C ~ p.ṁ * p.cp          # duty
        0 ~ q.Q̇ - Q̇
        n.T ~ p.T + Q̇ / C
        n.P ~ p.P
    ]
    ODESystem(eqs, t, [C], [Q̇]; name = name, systems = [p, n, q], defaults = [C => 5192])
end

function ActiveThermoCompressor(; name, η = 1.0, rp = 3.5)
    @named p = ThermoPin()
    @named n = ThermoPin()
    @named w = WorkPin()

    ps = @parameters rp = rp η = η
    eqs = [
        0 ~ p.ṁ + n.ṁ #p.ṁ ~ n.ṁ                               # conservation of mass
        n.P ~ p.P * rp
        n.T ~ p.T * ((1.0 - η - (rp^((p.k - 1) / p.k))) / (-η))
        w.Ẇ ~ p.ṁ * p.cp * (n.T - p.T)
        w.Ẇ ~ p.Φ + n.Φ
    ]
    ODESystem(eqs, t, [], ps; name = name, systems = [p, n, w])
end

function PassiveThermoCompressor(; name, η = 1.0)
    @named p = ThermoPin()
    @named n = ThermoPin()
    @named w = WorkPin()

    ps = @parameters η = η

    eqs = [
        0 ~ p.ṁ + n.ṁ #p.ṁ ~ n.ṁ                               # conservation of mass
        n.T ~ p.T * ((1.0 - η - ((n.P / p.P)^((p.k - 1) / p.k))) / (-η))
        w.Ẇ ~ p.ṁ * p.cp * (n.T - p.T)
        w.Ẇ ~ p.Φ + n.Φ
    ]
    ODESystem(eqs, t, [], ps; name = name, systems = [p, n, w])
end

function ActiveThermoTurbine(; name, η = 1.0, rp = 3.5)
    @named p = ThermoPin(Pdef = 80, Tdef = 800)
    @named n = ThermoPin(Pdef = 80 / 3.5, Tdef = 500)
    @named w = WorkPin()
    ps = @parameters rp = rp η = η
    eqs = [
        0 ~ p.ṁ + n.ṁ                               # conservation of mass
        p.P ~ n.P * rp
        n.T ~
            p.T * (η + rp^((p.k - 1) / p.k) - η * (rp^((p.k - 1) / p.k))) /
            (rp^((p.k - 1) / p.k))
        w.Ẇ ~ p.ṁ * p.cp * (n.T - p.T)
        w.Ẇ ~ p.Φ + n.Φ
    ]
    ODESystem(eqs, t, [], ps; name = name, systems = [p, n, w])
end

function PassiveThermoTurbine(; name, η = 1.0)
    @named p = ThermoPin(Pdef = 80, Tdef = 800)
    @named n = ThermoPin(Pdef = 80 / 3.5, Tdef = 500)
    @named w = WorkPin()
    ps = @parameters η = η
    sts = @variables rp(t) = 3.5
    eqs = [
        0 ~ p.ṁ + n.ṁ                               # conservation of mass
        p.P ~ n.P * rp
        n.T ~
            p.T * (η + rp^((p.k - 1) / p.k) - η * (rp^((p.k - 1) / p.k))) /
            (rp^((p.k - 1) / p.k))
        w.Ẇ ~ p.ṁ * p.cp * (n.T - p.T)
        w.Ẇ ~ p.Φ + n.Φ
    ]
    ODESystem(eqs, t, sts, ps; name = name, systems = [p, n, w], defaults = [rp => 3.5])
end

function IdealCooler(; name)
    @named p = ThermoPin()
    @named n = ThermoPin()
    @named q = HeatTransferPin()
    sts = @variables Q̇(t) = 0.0
    # 0 variables, 3 equations for pin
    eqs = [
        q.Q̇ ~ p.Φ + n.Φ             # conservation of energy
        0 ~ p.ṁ + n.ṁ
        n.P ~ p.P                 # no pressure   
        q.Q̇ ~ p.ṁ * p.cp * (n.T - p.T)
        Q̇ ~ q.Q̇
    ]
    ODESystem(eqs, t, sts, []; name = name, systems = [p, n, q], defaults = [Q̇ => 0])
end

function PassiveElement(; name)
    @named p = ThermoPin()
    @named n = ThermoPin()
    sts = @variables Q̇(t) = 0.0 ΔP(t) = 0.0 ΔT(t) = 0.0 Δṁ(t) = 0.0
    # 0 variables, 3 equations for pin
    eqs = [
        n.ṁ + p.ṁ ~ Δṁ    #has to be negative
        ΔT ~ p.T - n.T
        Q̇ ~ -(p.ṁ * p.cp * ΔT)
        ΔP ~ p.P - n.P
        Q̇ ~ p.Φ + n.Φ             # conservation of energy
    ]
    ODESystem(
        eqs,
        t,
        [Q̇, ΔP, ΔT, Δṁ],
        [];
        name = name,
        systems = [p, n],
        defaults = [Q̇ => 0, ΔT => 0, ΔP => 0, Δṁ => 0],
    )
end

function ReliefElement(; name)
    @named p = ThermoPin()
    @named n = ThermoPin()
    @named q = HeatTransferPin()
    # 0 variables, 3 equations for pin
    eqs = [
        n.ṁ + p.ṁ ~ 0           # has to be negative
        q.Q̇ ~ p.Φ + n.Φ         # conservation of energy
    ]
    ODESystem(eqs, t, [], []; name = name, systems = [p, n, q])
end

function throttle(; name)
    @named p = ThermoPin()
    @named n = ThermoPin()
    sts = @variables ΔP(t) = 0.0 ΔΦ(t) = 0.0
    eqs = [
        n.ṁ + p.ṁ ~ 0
        p.T ~ n.T
        ΔP ~ p.P - n.P
        n.Φ ~ 0
        ΔΦ ~ p.Φ + n.Φ
    ]              # conservation of energy
    ODESystem(eqs, t, sts, []; name = name, systems = [p, n], defaults = [ΔP => 0, ΔΦ => 0])
end

function Regenerator(; name, ϵ = 0.95)
    @named A = ThermoHeatTransfer()
    @named B = ThermoHeatTransfer()
    ps = @parameters ϵ = ϵ
    @variables Q̇(t) = 0.0 Cmin(t) = 0.0
    eqs = [
        Q̇ ~ ϵ * (A.p.T - B.p.T) * ((A.C < B.C) * A.C + (A.C >= B.C) * B.C)      #   heat transfer out of A -> B , if A/T > B/T 
        0 ~ A.Q̇ + B.Q̇
        A.Q̇ ~ -Q̇
    ]
    ODESystem(eqs, t, [Q̇], ps; name = name, systems = [A, B], defaults = [ϵ => 0.95])
end

## compound components
function Intercooler(; name, Tout = 300)
    @named p = ThermoPin()
    @named n = ThermoPin()

    # have to switch signs for internal circiuts
    sts = @variables Q̇(t) = 0.0

    @named st = SetTemperature(T = Tout)
    @named ic = IdealCooler()

    eqs = vcat(
        gas_CopyConnect(p, ic.p),
        gas_connect(ic.n, st.p),
        gas_CopyConnect(st.n, n),
        Q̇ ~ ic.Q̇,
    )
    ODESystem(
        eqs,
        t,
        [Q̇],
        [];
        name = name,
        systems = [p, n, st, ic],
        defaults = [Q̇ => 0.0],
    )
end

function showsol(c, sol)
    for cel in c
        cvec = []
        if hasproperty(cel, :p)
            push!(cvec, cel.p)
        end
        if hasproperty(cel, :n)
            push!(cvec, cel.n)
        end
        for n in cvec

            @printf "\t %s ṁ = %.2f \t P = % .2f \t T = %.2f " n.name sol[n.ṁ][end] sol[n.P][end] sol[n.T][end]

            if hasproperty(cel, :Q̇)
                @printf " Q̇ = %.2e " sol[cel.Q̇][end]
            elseif hasproperty(cel, :q)
                @printf " Q̇ = %.2e " sol[cel.q.Q̇][end]
            end
            if hasproperty(cel, :C)
                @printf " C = %.2e " sol[cel.C][end]
            end
            if hasproperty(cel, :w)
                @printf " Ẇ = %.2e " sol[cel.w.Ẇ][end]
            end
            @printf "\n"
        end
    end
end

###################################################
#                   TEST FUNCTIONS                
###################################################

function testRegen()
    @named hotres = SinglePortReservoir(P = 100, T = 500)
    @named hmsrc = GasFlowSource(Ṁ = 1.0)

    @named coolres = SinglePortReservoir(P = 100, T = 400)
    @named cmsrc = GasFlowSource(Ṁ = 2.0)

    @named reg = Regenerator()
    connections = vcat(
        gas_connect(hotres.n, hmsrc.p),
        gas_connect(hmsrc.n, reg.A.p),
        gas_connect(coolres.n, cmsrc.p),
        gas_connect(cmsrc.n, reg.B.p),
    )

    systemNames = [hotres, hmsrc, coolres, cmsrc, reg]

    @named odemodel = ODESystem(connections, t; systems = systemNames)

    smodel = substitute(odemodel, propDict)
    simple_mod = structural_simplify(smodel)
    prob = ODAEProblem(simple_mod, Pair[], (0.0, 1.0))
    sol = solve(prob)
    showsol([hotres, hmsrc, coolres, cmsrc, reg.A, reg.B], sol)
end

function testIntercool()
    @named hotres = SinglePortReservoir(P = 100, T = 400)
    @named hmsrc = GasFlowSource(Ṁ = 1.0)
    @named IC = Intercooler(Tout = 320)

    connections = vcat(gas_connect(hotres.n, hmsrc.p), gas_connect(hmsrc.n, IC.p))

    systemNames = [hotres, hmsrc, IC]

    @named odemodel = ODESystem(connections, t; systems = systemNames)

    smodel = substitute(odemodel, propDict)
    simple_mod = structural_simplify(smodel)
    prob = ODAEProblem(simple_mod, Pair[], (0.0, 1.0))
    sol = solve(prob)
    showsol(systemNames, sol)
end

function testTurbine(; Tin = 1000, Pin = 50, rp = 3.0, returnmode = :sys)
    @named hotres = SinglePortReservoir(P = Pin, T = Tin)
    @named hmsrc = GasFlowSource(Ṁ = 1.0)
    @named turbine = PassiveThermoTurbine(η = 0.9)
    @named ep = SetPressure(P = Pin / rp)

    connections = vcat(
        gas_connect(hotres.n, hmsrc.p),
        gas_connect(hmsrc.n, turbine.p),
        gas_connect(turbine.n, ep.p),
    )

    systemNames = [hotres, hmsrc, turbine, ep]

    @named odemodel = ODESystem(connections, t; systems = systemNames)

    if returnmode == :sys
        return odemodel
    else
        smodel = substitute(odemodel, propDict)

        simple_mod = structural_simplify(smodel)

        prob = ODAEProblem(simple_mod, Pair[], (0.0, 1.0))

        sol = solve(prob)

        showsol(systemNames, sol)
        return sol
    end
end

function testHeater(Tin = 300, Pin = 100, Qin = 15e6)
    @named cres = SinglePortReservoir(P = Pin, T = Tin)
    @named msrc = GasFlowSource(Ṁ = 75)
    @named heat = FixedThermoHeatTransfer(Qin = Qin)
    @named ep = SetPressure(P = 100)

    connections = vcat(
        gas_connect(cres.n, msrc.p),
        gas_connect(msrc.n, heat.p),
        gas_connect(heat.n, ep.p),
    )

    systemNames = [cres, msrc, heat]

    @named odemodel = ODESystem(connections, t; systems = systemNames)

    smodel = substitute(odemodel, propDict)
    simple_mod = structural_simplify(smodel)
    prob = ODAEProblem(simple_mod, Pair[], (0.0, 1.0))
    sol = solve(prob)
    showsol(systemNames, sol)
end

###################################################
#                   Cases            
###################################################
function simpleBrayton()
    @named cmsrc = GasFlowSource(Ṁ = 100)
    @named comp = ThermoCompressor()
    @named heat = ThermoHeatSource(Qin = 100e6)
    @named turbine = PassiveThermoTurbine()
    @named cool = IdealCooler()
    @named res = TwoPortReservoir(P = 10, T = 300)

    connections = vcat(
        gas_connect(res.n, cmsrc.p),
        gas_connect(cmsrc.n, comp.p),
        gas_connect(comp.n, heat.p),
        gas_connect(heat.n, turbine.p),
        gas_connect(turbine.n, cool.p),
        gas_connect(cool.n, res.p),
    )

    # systemNames =[coolres,cmsrc,comp,heat,turbine,cool,outres]
    systemNames = [res, cmsrc, comp, heat, turbine, cool]

    @named odemodel = ODESystem(connections, t; systems = systemNames)

    smodel = substitute(odemodel, propDict)
    simple_mod = structural_simplify(smodel)
    prob = ODAEProblem(simple_mod, Pair[], (0.0, 1.0))
    sol = solve(prob)
    @show Qnet = sol[cool.q.Q̇] + sol[heat.q.Q̇]
    @show Wnet = sol[turbine.w.Ẇ] + sol[comp.w.Ẇ]
    @show Wnet + Qnet
    @show Wnet / sol[heat.q.Q̇]

end

function simpleBraytonRegen()
    TminCycle = 300
    PminCycle = 15
    Qload = 1e6
    Mflow = 75

    @named res = TwoPortReservoir(P = PminCycle, T = TminCycle)
    @named valve = GasFlowSource(Ṁ = Mflow)
    @named comp1 = ActiveThermoCompressor(rp = 3.8, η = 0.9)
    @named ic1 = Intercooler(Tout = TminCycle)
    @named comp2 = ActiveThermoCompressor(rp = 3.5, η = 0.95)
    @named ic2 = Intercooler(Tout = TminCycle)
    @named comp3 = ActiveThermoCompressor(rp = 3.5, η = 0.95)
    @named reg = Regenerator()
    @named heat = FixedThermoHeatTransfer(Qin = 1000e6)
    @named turbine = PassiveThermoTurbine()
    @named cool = IdealCooler()

    connections = vcat(
        gas_connect(res.n, valve.p),
        gas_connect(valve.n, comp1.p),
        gas_connect(comp1.n, ic1.p),
        gas_connect(ic1.n, comp2.p),
        gas_connect(comp2.n, ic2.p),
        gas_connect(ic2.n, comp3.p),
        hx_connect(reg, comp3, heat, turbine, cool),
        gas_connect(heat.n, turbine.p),
        gas_connect(cool.n, res.p),
    )

    # systemNames =[coolres,valve,comp,heat,turbine,cool,outres]
    systemNames = [res, valve, comp1, ic1, comp2, ic2, comp3, reg, heat, turbine, cool]
    @named odemodel = ODESystem(connections, t; systems = systemNames)
    smodel = substitute(odemodel, propDict)
    simple_sys = structural_simplify(smodel)

end

function GasReference(; name, Pref = 10, T_ref = 500)
    # ground pin, Pref in Bar, Tref in K
    p = ThermoPin()
    eqs = [
        P.Φ ~ 0
        p.ṁ ~ 0
        p.P ~ Pref
        p.T ~ Tref
    ]

end
