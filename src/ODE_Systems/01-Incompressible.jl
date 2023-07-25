# module ThermoIncompressible
using ModelingToolkit, Plots, DifferentialEquations, Revise, CoolProp, Symbolics, Logging
using NonlinearSolve, Printf
using Plots
@variables t
Logging.disable_logging(Logging.Warn)
# include("03-MTK_UTILS.jl")
# data from [2]
cppb(x)     = (0.195 - 9.116e-6 .* x) .* 1000  # J/kgK
vpb(x)      = 1 / (10520.35 - 1.19051 .* x)      # m³/kg
spb(x)      = 0.195*1e3*log(x) - 9.11e-6*1000*x - 1210.316002277804
hpb(x,y)      = cppb(x)*x +  vpb(x) * y * 1e5 

Lcpfunc(x) = cppb(x)
Lvfunc(x) = vpb(x)
Lsfunc(x) = spb(x)
Lhfunc(x,y) = hpb(x,y)

@variables x
@register_symbolic Lcpfunc(x)
@register_symbolic Lvfunc(x)
@register_symbolic Lsfunc(x)
@register_symbolic Lhfunc(x,y)
propDict = Dict(Lcpfunc => cppb, 
            Lvfunc => vpb,
            Lsfunc => spb,
            Lhfunc => hpb)


@connector function IncompressiblePin(; name, Pdef = 50, Tdef = 555, ṁdef = 0.0)
    across_var  = @variables  P(t)=Pdef T(t)=Tdef s(t)=1.0 cp(t)=187 v(t)=.001
    thru_var    = @variables  ṁ(t)=ṁdef Φ(t)=1.0                     # mass flow and energy flow

    eq = [
        cp ~ Lcpfunc(T)
        v ~ Lvfunc(T)
        s ~ Lsfunc(T)
    ]
    sts = [T, P, ṁ, cp,s,v, Φ]
    ODESystem(eq, t, sts, []; name = name, defaults = [P => Pdef,T=>Tdef, ṁ =>ṁdef, Φ => 0.0])
end

@connector function WorkPin(; name)
    sts = @variables Ẇ(t) = 0.0 
    ODESystem(Equation[], t, sts, []; name = name)
end

@connector function HeatTransferPin(; name)  #input in W
    sts = @variables Q̇(t) = 0.0
    ODESystem(Equation[], t, sts, []; name = name, defaults = [Q̇ => 0.0])
end

@connector function FixedHeatFlowPin(; name, Qin = 1e3)  #input in W
    sts = @variables Q̇(t)=Qin 
    ps = @parameters Qin = Qin
    eqs=[Q̇ ~ Qin]
    ODESystem(eqs,t, sts, ps; name = name, defaults = [Q̇ => Qin])
end

@component function IncompressibleOnePort(; name)
    @named p = IncompressiblePin()
    @named n = IncompressiblePin()

    sts=@variables LHS(t) ΔP(t)=0.0 ΔT(t)=0.0

    eqs =[
        LHS ~ p.Φ + n.Φ             # conservation of energy
        0 ~ p.ṁ + n.ṁ               # mass flow 
        ΔP ~ n.P - p.P
        ΔT ~ n.T - p.T
    ]
    ODESystem(eqs,t, sts, []; name = name, systems = [p,n])
end

@component function SetPressure(;name, P = 0.1)
    @named oneport = IncompressibleOnePort()
    @unpack ΔT, ΔP, LHS, p = oneport
    ps = @parameters P = P 
    eqs = [
        LHS ~ 0
        ΔT ~ 0
        ΔP ~ 0
        p.P ~  P
        ]
        extend(ODESystem(eqs,t,[],ps; name = name), oneport)
end

@component function SetPressure2(;name, P = 0.1)
    @named p = IncompressiblePin(Pdef = P)
    @named n = IncompressiblePin(Pdef = P)
    ps = @parameters P = P 
    eqs = [
        p.P ~ P
        n.P ~ p.P
        n.T ~ p.T
        0 ~ p.Φ + n.Φ             # conservation of energy
        0 ~ p.ṁ + n.ṁ
    ]
    ODESystem(eqs, t,[], ps; name = name, systems = [p,n])
end

@component function SetTemperature(;name, T = 300)
    @named oneport = IncompressibleOnePort()
    @unpack ΔT, ΔP, LHS, p = oneport
    ps = @parameters T = T 
    eqs = [
        p.T ~ T
        ΔT ~ 0
        ΔP ~ 0
        LHS ~ 0
    ]
    extend(ODESystem(eqs,t,[],ps; name = name), oneport)
end

@component function IncompressibleFlowValve(; name)
    @named p = IncompressiblePin()
    @named n = IncompressiblePin()

    sts = @variables ṁ(t) = 1.0
    eqs = [
        0 ~ p.Φ + n.Φ   # conservation of energy
        0 ~ p.ṁ + n.ṁ
        p.ṁ ~ ṁ
        p.T ~ n.T 
        n.P ~ p.P
    ] 
    ODESystem(eqs,t,sts,[]; name = name,systems = [p,n],defaults = [ṁ => 1.0])
end

@component function IncompressibleFixedFlowSource(; name, Ṁ = 1.0)
    @named p = IncompressiblePin()
    @named n = IncompressiblePin()

    ps = @parameters Ṁ = Ṁ
    eqs = [
        0 ~ p.Φ + n.Φ   # conservation of energy
        0 ~ p.ṁ + n.ṁ
        p.ṁ ~ Ṁ
        p.T ~ n.T 
        n.P ~ p.P
    ] 
    ODESystem(eqs,t,[],ps; name = name,systems = [p,n])
end

@component function PassiveIncompressiblePump(; name, η = 0.9)
    @named oneport = IncompressibleOnePort()
    @named w = WorkPin()
    @unpack ΔP, LHS, p,n = oneport
    ps = @parameters η = η

    eqs = [LHS ~ w.Ẇ
            w.Ẇ ~ p.ṁ * p.v * ΔP * 1e5 / η
            n.T  ~ p.T + (w.Ẇ - p.v*ΔP)/p.cp # work, multiply by 100 to get to kPa
            ] 
    extend(ODESystem(eqs,t,[],ps;systems = [w], name = name), oneport)
end

@component function PassiveIncompressiblePump2(; name, η = 1.0)
    @named p = IncompressiblePin()
    @named n = IncompressiblePin()
    @named w = WorkPin()

    ps = @parameters η = η

    eqs = [
        0 ~ p.ṁ + n.ṁ #p.ṁ ~ n.ṁ                               # conservation of mass
        w.Ẇ ~ p.ṁ * p.v * (n.P - p.P) * 1e5 / η
        w.Ẇ ~ p.Φ + n.Φ
        n.T  ~ p.T + (p.v * (n.P - p.P) * 1e5 / η - p.v*(n.P - p.P)*1e5)/p.cp
    ]
    ODESystem(eqs,t,[],ps; name = name, systems = [p,n,w])
end



@component function IncompressibleHeatTransfer(; name)
    @named p = IncompressiblePin()
    @named n = IncompressiblePin()
    @named q = HeatTransferPin()

    st = @variables Q̇(t)=0.0 C(t)=187

    eqs = [
        0 ~ p.ṁ + n.ṁ           # p.ṁ ~ n.ṁ                               # conservation of mass
        q.Q̇ ~ p.Φ + n.Φ         # conservation of energy            
        C ~ p.ṁ * p.cp          # duty
        0 ~ q.Q̇ - Q̇ 
        n.T ~ p.T + q.Q̇/C
        n.P ~ p.P
    ]
    ODESystem(eqs,t,[Q̇,C],[]; name = name, systems = [p,n,q], defaults = [Q̇ => 0.0, C => 187])
end

@component function FlowControlIncompressibleHeatTransfer(; name, ΔP = 0.0, Tout = 1000.0)
    @named p = IncompressiblePin()
    @named n = IncompressiblePin()
    @named q = HeatTransferPin()

    st = @variables Q̇(t)=0.0 C(t)=187 ṁ(t)=1.0
    ps = @parameters ΔP = ΔP Tout = Tout
    eqs = [
        n.T ~ Tout
        ṁ ~ q.Q̇ /((n.T - p.T)*p.cp)           # q = mcp*ΔT
        p.ṁ ~ ṁ
        0 ~ p.ṁ + n.ṁ           # p.ṁ ~ n.ṁ                               # conservation of mass
        q.Q̇ ~ p.Φ + n.Φ         # conservation of energy            
        C ~ ṁ * p.cp          # duty
        0 ~ q.Q̇ - Q̇ 
        n.P ~ p.P - ΔP
    ]
    ODESystem(eqs,t,[Q̇,C,ṁ],ps; name = name, systems = [p,n,q], defaults = [Q̇ => 0.0, C => 187])
end

@component function SinglePortReservoir(;name, P = 0.1, T = 300)
    @named n = IncompressiblePin(Pdef = P, Tdef = T)
    ps = @parameters P=P T=T
    eqs = [
        n.P ~ P
        n.T ~ T
        n.Φ ~ 0
    ]
    ODESystem(eqs, t,[], ps; name = name, systems = [n])
end

@component function TwoPortReservoir(;name, P = 0.1, T = 300)
    @named oneport = IncompressibleOnePort()
    @unpack ΔT, ΔP, LHS, n,p = oneport
    ps=@parameters T=T P=P
    eqs = [LHS ~ 0
        n.Φ ~ 0
        p.T ~ T
        p.P ~ P
        ΔP ~ 0
        ΔT ~ 0]

    extend(ODESystem(eqs,t,[],ps; name = name, defaults = [T => 300.0, P => 10]), oneport)
end

@component function throttle(; name)
    @named oneport = IncompressibleOnePort()
    @unpack LHS,p,ΔT,n = oneport
    sts=@variables Ė_loss(t)=0.0
    eqs = [LHS ~ p.Φ
            ΔT ~ 0
            Ė_loss ~ p.Φ + n.Φ]
    extend(ODESystem(eqs,t,sts,[]; name = name, defaults = [Ė_loss => 0.0]), oneport)
end

@component function IdealCooler(;name)
    @named p = IncompressiblePin()
    @named n = IncompressiblePin()
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

@component function ReliefElement(;name)
    @named p = IncompressiblePin()
    @named n = IncompressiblePin()
    @named q = HeatTransferPin()
    # 0 variables, 3 equations for pin
    eqs = [
        n.ṁ + p.ṁ ~ 0           # has to be negative
        q.Q̇ ~ p.Φ + n.Φ         # conservation of energy
        ]
    ODESystem(eqs, t,[], []; name = name, systems = [p,n,q])
end

"""
    incompressible_connect(pins...)
    Adds mass and energy balance eqs to all nodes
"""
function incompressible_connect(pins...)

    eqs = [
        sum(pin -> pin.ṁ, pins) ~ 0.0, # Mass
        sum(pin -> pin.Φ, pins) ~ 0.0, # Energy
    ]

    for pin in pins[2:end]
        push!(eqs, pins[1].P ~ pin.P)   # p_1 = p_2, p_1 = p_3,...
        push!(eqs, pins[1].T~ pin.T)   # p_1 = p_2, p_1 = p_3,...
    end
    return eqs
end

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


# function breeder_loop()
#     @named res = TwoPortReservoir(P=10,T=500)
#     @named valve = IncompressibleFlowSource()
#     @named pump = PassiveIncompressiblePump()
#     @named pset = SetPressure(P=40)
#     @named qin  = IncompressibleHeaatTransfer()
#     @named qpsv = IdealCooler()
#     @named throt = throttle()
#     # @named res2 = SinglePortReservoir(P=40,T=500)

#     connections = vcat(incompressible_connect(res.n,valve.p),
#                         incompressible_connect(valve.n, pump.p),
#                         incompressible_connect(pump.n,pset.p),
#                         incompressible_connect(pset.n,qin.p),
#                         incompressible_connect(qin.n,qpsv.p),
#                         incompressible_connect(qpsv.n,throt.p),
#                         incompressible_connect(throt.n,res.p))

#     mass_flow_fcn(t) = 100 + 15 * sin(t)
#     Q̇in_fcn(t) = 100e6 + 35e6 * cos(t)  

#     @register_symbolic mass_flow_fcn(t)
#     @register_symbolic Q̇in_fcn(t)

#     controlled_eqs = [valve.p.ṁ ~ mass_flow_fcn(t),
#                             qin.Q̇ ~ Q̇in_fcn(t)]

#     systemNames = [res,valve,pump,pset,qin,qpsv,throt];#,throt];

#     @named aux_sys = ODESystem(controlled_eqs,t)
#     @named odemodel = ODESystem(connections,t; systems =  systemNames)
#     @named odemodel = extend(aux_sys,odemodel)

#     smodel = substitute(odemodel,propDict)
#     simple_sys = structural_simplify(smodel)

#     tspan = (0.0,5.0)

#     saveat = LinRange(1.0,4.0,10)
#     u0 = [valve.ṁ => 100]

#     kwargs = (abstol=1e-10, reltol=1e-2, saveat = saveat)
#     p = [500,10,1.0,40]

#     probae = ODAEProblem(simple_sys,[],tspan, p)
#     probde = ODEProblem(simple_sys,[],tspan, p)

#     solde   = solve(probde, Rodas4(); kwargs...);
#     solae   = solve(probae, Rodas4(); kwargs...);

#     p1=plot(solae, vars = [valve.p.ṁ])
#     p2=plot(solae, vars = [qin.Q̇, qpsv.Q̇])

#     p = plot(p1,p2,layout = (2,1))


#     display(p)
#     showsol(systemNames,sol)
#     sol.t
# end

# end