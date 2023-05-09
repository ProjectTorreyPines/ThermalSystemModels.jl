# module ThermoSteam
using ModelingToolkit, Plots, Revise, Unitful, Symbolics, Logging, XSteam, Printf
using Graphs, Plots, GraphRecipes
@variables t
Logging.disable_logging(Logging.Warn)

##
#       Need todo: Coolprops, validation
#               Connect work nodes - master energy balance?
##
#   PROPERTIES DICT
@register_symbolic stm_cpfunc(x,y)
@register_symbolic stm_cphfunc(x,y)
@register_symbolic stm_vphfunc(x,y)

@register_symbolic stm_Tpsfunc(x,y)
@register_symbolic stm_Tphfunc(x,y)

@register_symbolic stm_xphfunc(x,y)
@register_symbolic stm_xptfunc(x,y)

@register_symbolic stm_sptfunc(x,y)
@register_symbolic stm_sphfunc(x,y)

@register_symbolic stm_hptfunc(x,y)
@register_symbolic stm_hpsfunc(x,y)
@register_symbolic stm_hsatfunc(x)

hydro_prop_dict = Dict(stm_Tpsfunc => T_ps,
                        stm_Tphfunc => T_ph,
                        stm_xphfunc => x_ph,
                        stm_sptfunc => s_pT,
                        stm_sphfunc => s_ph,
                        stm_hptfunc => h_pT,
                        stm_hpsfunc => h_ps,
                        stm_hsatfunc => hL_p,
                        stm_vphfunc => v_ph,
                        stm_cphfunc => Cp_ph);

#   STEAM PINS
"""
    BasicSteamPin()
    Self computes T,s,x,V
    Must have methods for ṁ,Φ,P,h
"""
@connector function BasicSteamPin(; name, Pdef=0.1)
    across_var  = @variables  P(t)=Pdef T(t)=190 s(t)=0.0 h(t)=191 x(t)=0.0 v(t)=.001
    thru_var    = @variables  ṁ(t)=0.0 Φ(t)=0.0                     # mass flow and energy flow
    
    ps = []
    
    eqs = Equation[T ~ stm_Tphfunc(P,h)
                    s ~ stm_sphfunc(P,h)
                    x ~ stm_xphfunc(P,h)
                    v ~ stm_vphfunc(P,h)]

    ODESystem(eqs, t, [thru_var..., across_var...], ps; name = name, defaults = [P => Pdef, ṁ => 0, Φ => 0, h=>191])
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

"""
    hydro_connect(pins...)
    Adds mass and energy balance eqs to all nodes

"""
function hydro_connect(pins...)

    eqs = [
        sum(pin -> pin.ṁ, pins) ~ 0.0, # Mass
        sum(pin -> pin.Φ, pins) ~ 0.0, # Energy
    ]

    for pin in pins[2:end]
        push!(eqs, pins[1].P ~ pin.P)   # p_1 = p_2, p_1 = p_3,...
        push!(eqs, pins[1].h ~ pin.h)   # p_1 = p_2, p_1 = p_3,...
    end
    return eqs
end

function hydro_basic_connect(n,p)
    eqs = [
        0 ~ p.ṁ + n.ṁ
        0 ~ p.Φ + n.Φ
        p.h ~ n.h
        p.P ~ n.P
    ]
    return eqs
end

function hydro_series_connect(comp ;  returnmode = :eq)
    eqs =  hydro_basic_connect(comp[1].n,comp[2].p)   # n -> p

    if length(comp)>2
        for i = 2:length(comp)-1
            eqs = vcat(eqs,  hydro_basic_connect(comp[i].n,comp[i+1].p))
        end
    end

    if returnmode == :ode
        return ODESystem(eqs,t; name = :connections)
    end
    
    return eqs
end

function extenda(odevec::Vector{ODESystem})
    if length(odevec)==1
        return odevec[1]
    else
        return extend(odevec[1],extenda(odevec[2:end]))
    end
end

@component function HeatTransferPort(; name)
    @named p = HeatTransferPin()
    @named n = HeatTransferPin()
    eqs = [0 ~p.Q̇ + n.Q̇]
    @variables Q̇(t) 
    ODESystem(Equation[],t, [Q̇], []; name = name, systems = [p,n])
end

@component function Reservoir(;name, P = 0.1)
    @named n = BasicSteamPin(Pdef = P)
    ps = @parameters P=P 
    eqs = [
        n.P ~ P
        n.h  ~ stm_hsatfunc(P)
        n.Φ ~ 0
    ]
    ODESystem(eqs, t,[], ps; name = name, systems = [n], defaults = [P => 0.1])
end

@component function SuperHeatedReservoir(;name, P = 150, T = 600)
    @named n = BasicSteamPin(Pdef = P)
    ps = @parameters P=P T=T 
    eqs = [
        n.P ~ P
        0 ~ stm_hptfunc(P,T) + n.h
        n.Φ ~ -n.ṁ * n.h
    ]
    ODESystem(eqs, t,[], ps; name = name, systems = [n], defaults = [P => 0.1, T => T])
end

@component function ioReservoir(;name, P = 0.1, fixboth = false)
    @named p = BasicSteamPin()
    @named n = BasicSteamPin(Pdef = P)
    ps = @parameters P=P
    eqs = [
        n.P ~ P
        n.h  ~ stm_hsatfunc(P)
        n.Φ ~ 0
        0 ~ n.ṁ + p.ṁ
    ]
    if fixboth
        push!(eqs,p.P ~ P)
        # push!(eqs,p.h ~ n.h)
    end
    ODESystem(eqs, t,[], ps; name = name, systems = [n,p], defaults = [P => 0.1 ])
end

@component function SetPressure(;name, P = 0.1)
    @named p = BasicSteamPin(Pdef = P)
    @named n = BasicSteamPin(Pdef = P)
    ps = @parameters P=P 
    eqs = [
        p.P ~ P
        n.P ~ p.P
        n.h ~ p.h
        0 ~ p.Φ + n.Φ             # conservation of energy
        0 ~ p.ṁ + n.ṁ
    ]
    ODESystem(eqs, t,[], ps; name = name, systems = [p,n], defaults = [P => 0.1])
end

@component function SteamFlowSource(;name, ṁ = 1.0)
    @named p = BasicSteamPin()
    @named n = BasicSteamPin()
    ps = @parameters Ṁ = ṁ
    eqs = [
            0 ~ p.Φ + n.Φ   # conservation of energy
            0 ~ p.ṁ + n.ṁ
            p.ṁ ~ Ṁ
            p.h ~ n.h 
            n.P ~ p.P
        ]
    ODESystem(eqs, t, [], ps; name = name, systems = [n,p], defaults = [Ṁ => ṁ])
    # compose(ODESystem(eqs, t, [], ps; name = name, defaults = [Ṁ => 1.0]),p,n)
end

@component function SteamFlowValve(;name)
    @named p = BasicSteamPin()
    @named n = BasicSteamPin()
    sts = @variables Ṁ(t)
    eqs = [
            0 ~ p.Φ + n.Φ   # conservation of energy
            0 ~ p.ṁ + n.ṁ
            p.ṁ ~ Ṁ
            p.h ~ n.h 
            n.P ~ p.P
        ]
    ODESystem(eqs, t,sts, []; name = name, systems = [n,p], defaults = [Ṁ => 1.0])
    # compose(ODESystem(eqs, t, [], ps; name = name, defaults = [Ṁ => 1.0]),p,n)
end

@component function AdiabaticPump(;name, η = 0.6, setpressure = false ,Pout = 10)
    #Pout in bar
    @named p = BasicSteamPin()
    @named n = BasicSteamPin()
    @named w = WorkPin()
    
    ps = @parameters η = η P=Pout

    eqs = Equation[
        w.Ẇ ~ p.Φ + n.Φ                         # conservation of energy
        0 ~ p.ṁ + n.ṁ
        n.h  ~ p.h + p.v *100* (n.P - p.P) / η  # work, multiply by 100 to get to kPa
        w.Ẇ  ~ p.ṁ * (n.h - p.h)
    ]

    if setpressure
        eqs = vcat(eqs,n.P  ~ P)
    end
    # extenda([ODESystem(eqs, t,[], ps; name = name, defaults = [η => 0.6, P => Pout]), p,n,w])
    ODESystem(eqs, t,[], ps; name = name, systems = [p,n,w], defaults = [η => .6, P => 10])
end

@component function Splitter(;name)
    @named p = BasicSteamPin()
    @named y = BasicSteamPin()
    @named z = BasicSteamPin()

    #  0 ~ pMax₊P - trbsimo₊split₊z₊P(t)
    eqs = [
        p.h ~ y.h               # enthalpy
        0 ~ y.h - z.h
        p.P ~ y.P               # pressure
        z.P ~ y.P 
        0 ~ p.ṁ + y.ṁ + z.ṁ
        0 ~ p.Φ + y.Φ + z.Φ
    ]

    ODESystem(eqs, t,[], []; name = name, systems = [p,y,z])
end

@component function AdiabaticTurbine(;name, η = 1.0, setpressure = false, Pout = 0.1)
    #Pout in bar
    @named p = BasicSteamPin()
    @named n = BasicSteamPin()
    @named w = WorkPin()
    ps = @parameters P=Pout η = η

    eqs = Equation[
        0 ~ p.ṁ + n.ṁ      
        n.h  ~ p.h - (p.h-stm_hpsfunc(n.P,p.s))*η
        w.Ẇ  ~ p.ṁ * (n.h - p.h)
        w.Ẇ ~ p.Φ + n.Φ                         # conservation of energy
    ]

    if setpressure
        eqs = vcat(eqs, n.P ~ P)
    end

    # extenda([ODESystem(eqs, t,[], ps; name = name, defaults = [η => 0.6, P => Pout]), p,n,w])
    ODESystem(eqs, t,[], ps; name = name, systems = [p,n,w], defaults = [η => 1.0, P => 10])
end   

@component function SIMOAdiabaticTurbine(;name, ηin = 1.0, setpressure = false, Pyin = 10, Pzin  = 0.1)
    #Pout in bar
    ps = @parameters Py=Pyin Pz=Pzin η = ηin
    sp = setpressure

    # EXTERNAL NODES FOR INTERFACING
    @named p = BasicSteamPin()  # inlet node
    @named w = WorkPin()
    @named hp = AdiabaticTurbine(η = ηin, Pout = Pyin, setpressure = sp)
    @named lp = AdiabaticTurbine(η = ηin, Pout = Pzin, setpressure = sp)

    #
    #                yn(-) -------> y
    # IN --> p (+) --|
    #                zn(-) -------> z
    #
    # INTERNAL NODES
    @named yn = BasicSteamPin()
    @named zn = BasicSteamPin()

    split_connect = hydro_connect(p, yn, zn)    # mflow, p.ṁ => positive
    hp_connect = hydro_connect(yn,hp.p)
    lp_connect = hydro_connect(zn,lp.p)


    eqs = vcat(split_connect,hp_connect,lp_connect)

    compose(ODESystem(eqs, t,[], ps; name = name, systems = [yn,zn,p], defaults = [η => 1.0, Py => 10, Pz => 0.1]), hp,lp)
    # ODESystem(eqs, t,[], ps; name = name, systems = [p,hp,lp], defaults = [η => 1.0, Py => 10, Pz => 0.1])
    # extend(ODESystem(eqs, t,sts, ps; name = name, systems = [hp,lp], defaults = [η => 1.0 Py => 10 Pz => 0.1]),split)
end

@component function SteamHeatTransfer(;name)
    @named p = BasicSteamPin()
    @named n = BasicSteamPin()
    @named q = HeatTransferPin()
    st = @variables Q̇(t)=0.0 C(t)=4000
    eqs = [
        0 ~ p.ṁ + n.ṁ           # p.ṁ ~ n.ṁ                               # conservation of mass
        q.Q̇ ~ p.Φ + n.Φ         # conservation of energy            
        C ~ p.ṁ * stm_cphfunc(p.P, p.h) * 1000   # duty T/h = cp
        0 ~ q.Q̇ - Q̇ 
        n.h ~ p.h + q.Q̇/(p.ṁ)/1000
        n.P ~ p.P
    ]
    ODESystem(eqs,t,[Q̇,C],[]; name = name, systems = [p,n,q], defaults = [Q̇ => 0.0, C => 400])
end

@component function IdealBoiler(;name, Tout = 350)
    #Pout in bar
    @named p = BasicSteamPin()
    @named n = BasicSteamPin()
    @named q = HeatTransferPin()
    
    ps = @parameters T = Tout

    eqs = Equation[
        q.Q̇ ~ p.Φ + n.Φ             # conservation of energy
        0 ~ p.ṁ + n.ṁ
        n.P  ~  p.P                 # no pressure
        n.h  ~ stm_hptfunc(p.P,T)       # work, multiply by 100 to get to kPa
        q.Q̇ ~ p.ṁ * (n.h - p.h)
    ]
    ODESystem(eqs, t,[], ps; name = name, systems = [p,n,q], defaults = [ T => Tout])
end

@component function IdealCondensor(;name)
    #Pout in bar
    @named p = BasicSteamPin()
    @named n = BasicSteamPin()
    @named q = HeatTransferPin()

    sts = @variables Q̇(t)=0.0
    
    eqs = Equation[
        q.Q̇ ~ p.Φ + n.Φ             # conservation of energy
        0 ~ p.ṁ + n.ṁ
        n.P  ~  p.P                 # no pressure
        n.h  ~ stm_hsatfunc(p.P)        # work, multiply by 100 to get to kPa
        q.Q̇ ~ p.ṁ * (n.h - p.h)
        q.Q̇ ~ Q̇
    ]
    ODESystem(eqs, t,sts, []; name = name, systems = [p,n,q], defaults = [Q̇ => 0.0])
end

@component function OpenFeedwaterHeater(;name)
    # flows x and y are the inlets
    @named y = BasicSteamPin()
    @named z = BasicSteamPin()
    @named n = BasicSteamPin()  

    sts=@variables yfrac(t) = 0.5


    continuous_events = [[yfrac ~ 1] => [yfrac ~ 1]
                        [yfrac ~ 0] => [yfrac ~ 0.0]]
    eqs  =[
        n.P ~ z.P
        n.h ~ stm_hsatfunc(n.P) # sat liquid
        0 ~ n.h - (yfrac * y.h + (1-yfrac) * z.h)
        0 ~ n.ṁ*yfrac + y.ṁ
        0 ~ n.ṁ + y.ṁ + z.ṁ
        0 ~ n.Φ + y.Φ + z.Φ
    ]
    ODESystem(eqs, t,sts, []; name = name, systems = [y,z,n], defaults = [yfrac => 0.5], continuous_events)
end
# end
@component function throttle(;name)
    # flows x and y are the inlets
    @named p = BasicSteamPin()
    @named n = BasicSteamPin()  

    sts=@variables Ė_loss(t) = 0.5 ΔP(t) = 0.0
    eqs  =[
        n.ṁ + p.ṁ ~ 0            
        p.h ~ n.h
        ΔP ~ p.P - n.P
        Ė_loss ~ p.Φ + n.Φ               # conservation of energy
        ]
    ODESystem(eqs, t,sts, []; name = name, systems = [p,n], defaults = [Ė_loss=> 0.0, ΔP => 0.0])
end

function plotmod(smodel)
    # structural_simplify(model)

    var2idx = Dict(s => i for (i, s) in enumerate(states(smodel)))
    idx2var = Dict(i => s for (i, s) in enumerate(states(smodel)))

    varvar = varvar_dependencies(asgraph(smodel), variable_dependencies(smodel))
        
    eqgraph = asgraph(equation_dependencies(smodel),var2idx)
    equation_dep_digraph = asdigraph(eqgraph,smodel)
    equation_dep_digraph.fadjlist
    varvar = varvar_dependencies(asgraph(smodel), variable_dependencies(smodel))
    varvar.badjlist
    # # asd = asdigraph(digr,odesys)
    nm = [idx2var[idx] for idx in 1:nv(varvar)]
    p=plot(varvar,
            size = (1500,1500),
            names = nm,
            method = :sfpd)
    

    return p,var2idx,idx2var,varvar,eqgraph
        
end

    # function FeedwaterRankine()
    #     @named iores        = ioReservoir(P=0.1)
    #     @named valve        = SteamFlowSource(ṁ = 1)
    #     @named pumpB        = AdiabaticPump(Pout = 150, setpressure = true)
    #     @named boil         = IdealBoiler(Tout = 600)
    #     @named turbine      = SIMOAdiabaticTurbine(setpressure = true, Pyin = 10, Pzin = 0.1,ηin = 1.0) 
    #     @named pumpA        = AdiabaticPump(Pout = 10, setpressure = true)
    #     @named cond         =  IdealCondensor()
    #     @named openfw       =  OpenFeedwaterHeater()

    #     connections = vcat(hydro_connect(pumpB.n,valve.p),
    #                         hydro_connect(valve.n, boil.p),
    #                         hydro_connect(boil.n, turbine.p),
    #                         hydro_connect(turbine.hp.n, openfw.y),
    #                         hydro_connect(turbine.lp.n, cond.p),
    #                         hydro_connect(cond.n,iores.p),
    #                         hydro_connect(iores.n,pumpA.p),
    #                         hydro_connect(pumpA.n,openfw.z),
    #                         hydro_connect(openfw.n,pumpB.p))

    #     systemNames =[valve,
    #                     boil,
    #                     turbine,
    #                     pumpB,
    #                     iores,
    #                     pumpA,
    #                     cond,
    #                     openfw]
                        
    #     @named odemodel = ODESystem(connections,t; systems = systemNames)
    #     smodel = substitute(odemodel,hydro_prop_dict)
    #     simple_mod = structural_simplify(smodel)
    #     prob = ODAEProblem(simple_mod,Pair[],(0.0,1.0))

    #     plotmod(smodel)
    #     prob.f.observed(openfw.n.P,prob.u0,prob.p,1.0)

    #     pf(xin) = prob.f.observed(xin,prob.u0,prob.p,1.0)
    #     nodes = [turbine.p, turbine.yn, turbine.zn, turbine.hp.p, turbine.lp.p,turbine.hp.n, turbine.lp.n,  cond.p, cond.n, pumpA.n,pumpB.p, openfw.z, openfw.y]
    #     for n in nodes
    #         @printf "\t %s ṁ = %.2f \t P = % .2f \t h = %.2f\n" n.name pf(n.ṁ) pf.(n.P) pf(n.h)
    #     end
    # end
    
    # export FeedwaterRankine
    # begin
    #     @named res          = Reservoir(P = 10)
    #     @named iores        = ioReservoir(P=10)
    #     @named freenodeA    = BasicSteamPin(Pdef = 150)
    #     @named valve        = SteamFlowSource(ṁ = 1)
    #     @named pumpB        = AdiabaticPump(Pout = 150, setpressure = true)
    #     @named boil         = IdealBoiler(Tout = 600)
    #     @named turbine      = SIMOAdiabaticTurbine(setpressure = true, Pyin = 10, Pzin = 0.1,ηin = 1.0) 
    #     @named pump         = AdiabaticPump(Pout = 10, setpressure = true)
    #     @named cond         = IdealCondensor()
    #     @named openfw       = OpenFeedwaterHeater()

    #     connections = vcat(hydro_connect(res.n, valve.p),
    #                         hydro_connect(valve.n,boil.p),
    #                         hydro_connect(boil.n, turbine.p),
    #                         hydro_connect(turbine.hp.n, openfw.y),
    #                         hydro_connect(turbine.lp.n, cond.p),
    #                         hydro_connect(cond.n,pump.p),
    #                         hydro_connect(pump.n,openfw.z),
    #                         hydro_connect(openfw.n,freenodeB),
    #                         hydro_connect(freenodeB,pumpB.p),
    #                         hydro_connect(pumpB.n,iores.p))

    #     systemNames =[res,valve,boil,turbine,pump,cond,openfw,freenodeB]

    #     @named odemodel = ODESystem(connections,t; systems = systemNames)

    #     smodel = substitute(odemodel,hydro_prop_dict)
    #     simple_mod = structural_simplify(smodel)
    #     prob = ODAEProblem(simple_mod,Pair[],(0.0,1.0))
    # end

    # begin
    #     println("Turbine")
    #     nodes = [turbine.p, turbine.yn, turbine.zn, turbine.hp.p, turbine.lp.p, cond.p, cond.n, openfw.z, openfw.y]
    #     for n in nodes
    #         @printf "\t %s ṁ = %.2f \t P = % .2f\n" n.name pf(n.ṁ) pf.(n.P)
    #     end


# end