using ModelingToolkit, Plots, DifferentialEquations, Revise, ModelingToolkitStandardLibrary, Unitful
using Symbolics, Logging

@variables t


cphe(x) = 5.192
cvhe(x) = 5.192 - 2.0769
@variables x
@register_symbolic cpfunc(x)
@register_symbolic kfunc(x)


# display(methods(cpfunc))
@connector function ThermoPin(; name, Pdef = 10.0, Tdef = 300, ṁdef = 0.0)
    @variables P(t) = Pdef             #[unit = u"bar"] 
    @variables ṁ(t) = ṁdef             #[unit = u"kg/s"] 
    @variables T(t)  = Tdef          #[unit = u"K"]
    @variables cp(t) = 5192
    @variables k(t)  = 1.667
    eq = [
        cp ~ cpfunc(T)
        k ~ kfunc(T)
    ]
    sts = [T, P, ṁ, cp, k]
    ODESystem(eq, t, sts, []; name = name, defaults = [P => Pdef, T => Tdef, ṁ => ṁdef] )
end

@connector function WorkPin(; name)
    sts = @variables Ẇ(t) = 0.0 #[output = true]
    ODESystem(Equation[], t, sts, []; name = name)
end

# convention = +Q = Qin
@connector function HeatFlowPin(; name, Qin = 1e6)  #input in W
    sts = @variables Q̇(t)=Qin #[input = true]
    ps = @parameters Qin = Qin
    eqs=[
        Q̇ ~ Qin
    ]
    ODESystem(eqs,t, sts, ps; name = name, defaults = [Q̇ => Qin])
end

@component function FixedFlowSource(;name, mflow = 1.0, P = 10, T = 300)
    @named src = ThermoPin(Pdef = P, Tdef = T, ṁdef = mflow)
    ps = @parameters ṁ=mflow P=P T=T
    # 0 variables, 3 equations for pin
    eqs = [
        src.T ~ T
        src.P ~ P
        ṁ + src.ṁ ~ 0   #has to be negative
    ]
    ODESystem(eqs, t, [], ps; name = name, systems = [src], defaults = [P => P, T =>T, ṁ => mflow] )
end

@component function ThermoHeat(; name)
    @named p = ThermoPin()
    @named n = ThermoPin()
    @named h = HeatFlowPin()
    #   0 variables 3 equations for pin
    eqs = [
        0 ~ p.ṁ + n.ṁ# p.ṁ ~ n.ṁ                               # conservation of mass
        n.T ~ p.T + h.Q̇/(p.ṁ*p.cp)
        n.P ~ p.P
    ]
    ODESystem(eqs,t,[],[]; name = name, systems = [p,n,h])
end

@component function ThermoCompressor(; name, η = 1.0, rp = 3.5)
    @named p = ThermoPin()
    @named n = ThermoPin()
    @named w = WorkPin()

    ps = @parameters rp = rp η = η
    eqs = [
        0 ~ p.ṁ + n.ṁ #p.ṁ ~ n.ṁ                               # conservation of mass
        n.P ~ p.P * rp
        n.T ~ p.T * ((1.0 - η - (rp^((p.k-1)/p.k))) / (-η))
        w.Ẇ ~ p.ṁ * p.cp * (n.T - p.T)
    ]
    ODESystem(eqs,t,[],ps; name = name, systems = [p,n,w], defaults = [η => η, rp => rp])
end

@component function ThermoTurbine(; name, η = 1.0, rp = 3.5)
    @named p = ThermoPin(Pdef = 80, Tdef = 800)
    @named n = ThermoPin(Pdef = 80/3.5, Tdef = 500)
    @named w = WorkPin()
    ps = @parameters rp = rp η = η
    eqs = [
        0 ~ p.ṁ + n.ṁ                               # conservation of mass
        n.P ~ p.P * rp
        n.T ~ p.T *(η + rp^((p.k-1)/p.k)- η*(rp^((p.k-1)/p.k))) / (rp^((p.k-1)/p.k))
        w.Ẇ ~ p.ṁ * p.cp * (n.T - p.T)
    ]
    ODESystem(eqs,t,[],ps; name = name, systems =  [p,n,w],  defaults = [η => η, rp => rp])
end

@component function WorkSensor(;name)
    @named src = WorkPin()
    var = @variables Ẇ(t)
    eqs= [0 ~ Ẇ - src.Ẇ]
    ODESystem(eqs,t,var,[]; name = name, systems =  [src])
end

@component function PassiveElement(;name, Pdef = 10, Tdef = 300)
    @named p = ThermoPin()
    @named n = ThermoPin(Pdef = Pdef, Tdef = Tdef)
    sts = @variables Q̇(t) = 0.0
    ps = @parameters P = Pdef T = Tdef
    # 0 variables, 3 equations for pin
    eqs = [
        n.ṁ + p.ṁ ~ 0   #has to be negative
        Q̇ ~ -(p.ṁ * p.cp * (p.T - n.T))
        ]
    ODESystem(eqs, t, [Q̇], ps; name = name, systems = [p,n], defaults = [P => Pdef, T =>Tdef] )
end

function connect_pins(pins...)
    # mass balance
    eqs = [
        sum(pin -> pin.ṁ, pins) ~ 0.0, # sum(mdot) = 0
    ]

    # pressure balance
    for pin in pins[2:end]
        push!(eqs, pins[1].P ~ pin.P) # p_1 = p_2, p_1 = p_3,...
    end

    return eqs
end


"""
    basic_connect(pinA,pinB) 
    Stuff
"""
function basic_connect(pinA,pinB)   # creates thermoPin connections  P = P, T = T, mflow = mflow
    eqs = [
        sum(pin -> pin.ṁ, [pinA,pinB]) ~ 0.0, # sum(mdot) = 0
    ]
    push!(eqs, pinA.P ~ pinB.P)
    push!(eqs, pinA.T ~ pinB.T)
    return eqs
end


"""
    Dont use this for flowsource
"""
function series_connect(pins)
    @assert length(pins)>2

    eqs =  basic_connect(pins[1].n,pins[2].p)

    for i = 2:length(pins)-1
        eqs =     vcat(eqs,  basic_connect(pins[i].n,pins[i+1].p))
    end
    return eqs
end


# mflow(t) = 100*t
@register_symbolic mf(t)

flow    = FixedFlowSource(name = :flow, mflow = mf(t))
cmp     = ThermoCompressor(name = :cmp)
heat    = ThermoHeat(name = :heat)
trb     = ThermoTurbine(name = :trb)
psv     = PassiveElement(name = :psv, Tdef = 300, Pdef = 300)

src_connections = basic_connect(flow.src,cmp.p)
all_connections = vcat(src_connections,series_connect([cmp,heat,trb]))

cphe(x) = 5192
khe(x) = 1.667
propDict = Dict(cpfunc => cphe, kfunc => khe)


@named model = ODESystem(all_connections,t; systems = [flow,cmp,heat,trb])
model = substitute(model,propDict)
sys = alias_elimination(model)
# equations(sys)
# sys.defaults
# sys.states


params = [
        flow.T    => 300
        flow.P    => 10.0
        flow.ṁ    => 10
        cmp.rp  => 3.5
        cmp.η   => 0.8
        heat.h.Qin => 1e6
        trb.rp    => 3.5
        trb.η     => 0.9
        ]


mod = structural_simplify(model)
tspan = (0.0,10.0)
prob = ODEProblem(mod,mod.defaults,tspan,params)
@show win = prob.f.observed(cmp.w.Ẇ, prob.u0, prob.p, 9.0)
# @show prob.f.observed(th.p.T, prob.u0, prob.p, 9.0)
# @show prob.f.observed(th.n.T, prob.u0, prob.p, 9.0)
@show prob.f.observed(cmp.p.T, prob.u0, prob.p, 9.0)
@show prob.f.observed(cmp.n.T, prob.u0, prob.p, 9.0)


@show qTi = prob.f.observed(heat.p.T, prob.u0, prob.p, 9.0)
@show qTo = prob.f.observed(heat.n.T, prob.u0, prob.p, 9.0)
@show qin = prob.f.observed(heat.h.Q̇, prob.u0, prob.p, 9.0)
# @show prob.f.observed(th.p.T, prob.u0, prob.p, 9.0)
@show   twork = prob.f.observed(trb.w.Ẇ, prob.u0, prob.p, 9.0)
@show   prob.f.observed(trb.n.T, prob.u0, prob.p, 9.0)
# @show  qwst = prob.f.observed(psv.Q̇, prob.u0, prob.p, 9.0)
@show qin + win 
@show twork