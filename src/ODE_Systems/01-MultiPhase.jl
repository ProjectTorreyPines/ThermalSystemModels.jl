# module ThermoSteam
using ModelingToolkit, Plots, Revise, Unitful, Symbolics, Logging, XSteam, Printf
using Graphs, Plots, GraphRecipes
using DifferentialEquations
@variables t
Logging.disable_logging(Logging.Warn)
include("03-MTK_UTILS.jl")
##
#       Need todo: Coolprops, validation
#               Connect work nodes - master energy balance?
##
#   PROPERTIES DICT
# @register_symbolic stm_cpfunc(x,y)
# @register_symbolic stm_cphfunc(x,y)
# @register_symbolic stm_vphfunc(x,y)

# @register_symbolic stm_Tpsfunc(x,y)
# @register_symbolic stm_Tphfunc(x,y)

# @register_symbolic stm_xphfunc(x,y)
# @register_symbolic stm_xptfunc(x,y)

# @register_symbolic stm_sptfunc(x,y)
# @register_symbolic stm_sphfunc(x,y)

# @register_symbolic stm_hptfunc(x,y)
# @register_symbolic stm_hpsfunc(x,y)
# @register_symbolic stm_hsatfunc(x)
#   PROPERTIES DICT
# converting XSteam to SI units
# W_cpfunc(x,y)  = Cp_ph(x,y/1e3)*1e3
# W_cphfunc(x,y) = Cp_ph(x,y/1e3)*1e3
# W_vphfunc(x,y) = v_ph(x, y/1e3)
# W_Tpsfunc(x,y) = T_ps(x, y/1e3)+273.15
# W_Tphfunc(x,y) = T_ph(x, y/1e3)+273.15
# W_xphfunc(x,y) = x_ph(x, y/1e3)
# W_xptfunc(x,y) = x_pT(x, y-273.15)
# W_sptfunc(x,y) = s_pT(x, y-273.15)*1e3
# W_sphfunc(x,y) = s_ph(x, y/1e3)*1e3
# W_hptfunc(x,y) = h_pT(x, y-273.15)*1e3
# W_hpsfunc(x,y) = h_ps(x,y/1e3)*1e3
# W_hsatfunc(x)  = hL_p(x)*1e3

stm_cpfunc(x,y) = Cp_ph(x,y/1e3)*1e3
stm_cphfunc(x,y)= Cp_ph(x,y/1e3)*1e3
stm_vphfunc(x,y)= v_ph(x, y/1e3)
stm_Tpsfunc(x,y)= T_ps(x, y/1e3)+273.15
stm_Tphfunc(x,y)= T_ph(x, y/1e3)+273.15
stm_xphfunc(x,y)= x_ph(x, y/1e3)
stm_xptfunc(x,y)= x_pT(x, y-273.15)
stm_sptfunc(x,y)= s_pT(x, y-273.15)*1e3
stm_sphfunc(x,y)= s_ph(x, y/1e3)*1e3
stm_hptfunc(x,y)= h_pT(x, y-273.15)*1e3
stm_hpsfunc(x,y)= h_ps(x,y/1e3)*1e3
stm_hsatfunc(x) = hL_p(x)*1e3


@variables x , y
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


# hydro_prop_dict=Dict(stm_cphfunc => W_cphfunc ,
#                     stm_vphfunc  => W_vphfunc ,
#                     stm_Tpsfunc  => W_Tpsfunc ,
#                     stm_Tphfunc  => W_Tphfunc ,
#                     stm_xphfunc  => W_xphfunc ,
#                     stm_xptfunc  => W_xptfunc ,
#                     stm_sptfunc  => W_sptfunc ,
#                     stm_sphfunc  => W_sphfunc ,
#                     stm_hptfunc  => W_hptfunc ,
#                     stm_hpsfunc  => W_hpsfunc ,
#                     stm_hsatfunc  => W_hsatfunc)


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
                        


# @variables xx, yy
# hst(x,y)=hydro_prop_dict(stm_hptfunc(x,y))
# hst(30,500)
# # hydro_prop_arr = [stm_Tpsfunc => T_ps,
# #                         stm_Tphfunc => T_ph,
#                         stm_xphfunc => x_ph,
#                         stm_sptfunc => s_pT,
#                         stm_sphfunc => s_ph,
#                         stm_hptfunc => h_pT,
#                         stm_hpsfunc => h_ps,
#                         stm_hsatfunc => hL_p,
#                         stm_vphfunc => v_ph,
#                         stm_cphfunc => Cp_ph];

    # hydro_prop_dict = Dict(stm_cpfunc  => PropsSI("C","P",x*u"bar","T",y*u"K","Water").val,
    #                             stm_cphfunc => PropsSI("C","P",x*u"bar","H",y*u"J/kg/K","Water").val,
    #                             stm_vphfunc => 1/PropsSI("D","P",x*u"bar","H",y*u"J/kg","Water").val,
    #                             stm_Tpsfunc => PropsSI("T","P",x*u"bar","S",y*u"J/kg","Water").val,
    #                             stm_Tphfunc => PropsSI("T","P",x*u"bar","H",y*u"J/kg","Water").val,
    #                             stm_xphfunc => x_ph(x,y/1e3),
    #                             stm_xptfunc => x_pT(x,y-273.15),
    #                             stm_sptfunc => PropsSI("SMASS","P",x*u"bar","T",y*u"K","Water").val,
    #                             stm_sphfunc => PropsSI("SMASS","P",x*u"bar","H",y*u"J/kg","Water").val,
    #                             stm_hptfunc => PropsSI("HMASS","P",x*u"bar","T",y*u"K","Water").val,
    #                             stm_hpsfunc => PropsSI("HMASS","P",x*u"bar","S",y*u"J/kg","Water").val,
    #                             stm_hsatfunc => hL_p(x))
    # hydro_prop_dict = Dict(stm_cphfunc => PropsSI("CPMASS","P",x*u"bar","H",y*u"J/kg","Water").val,
    #                             stm_vphfunc => 1/PropsSI("D","P",x*u"bar","H",y*u"J/kg","Water").val,
    #                             stm_Tpsfunc => PropsSI("T","P",x*u"bar","S",y*u"J/kg","Water").val,
    #                             stm_Tphfunc => PropsSI("T","P",x*u"bar","H",y*u"J/kg","Water").val,
    #                             stm_xphfunc => x_ph(x,y/1e3),
    #                             stm_xptfunc => x_pT(x,y-273.15),
    #                             stm_sptfunc => PropsSI("SMASS","P",x*u"bar","T",y*u"K","Water").val,
    #                             stm_sphfunc => PropsSI("SMASS","P",x*u"bar","H",y*u"J/kg","Water").val,
    #                             stm_hptfunc => PropsSI("HMASS","P",x*u"bar","T",y*u"K","Water").val,
    #                             stm_hpsfunc => PropsSI("HMASS","P",x*u"bar","S",y*u"J/kg","Water").val,
    #                             stm_hsatfunc => hL_p(x))
    # hydro_prop_dict = Dict(stm_Tpsfunc => PropsSI("T","P",x*u"bar","S",y*u"J/kg/K","Water").val,,
    #                         stm_Tphfunc => PropsSI("T","P",x*u"bar","H",y*u"J/kg/K","Water").val,
    #                         stm_xphfunc => x_ph,
    #                         stm_sptfunc => PropsSI("SMASS","P",30u"bar","T",350u"K","Water").val,
    #                         stm_sphfunc => PropsSI("SMASS","P",x*u"bar","H",y*u"J/kg/K","Water").val,
    #                         stm_hptfunc => PropsSI("HMASS","P",x*u"bar","T",y*u"K","Water").val,
    #                         stm_hpsfunc => PropsSI("HMASS","P",x*u"bar","S",y*u"J/kg/K","Water").val,
    #                         stm_hsatfunc => hL_p,
    #                         stm_vphfunc => v_ph)

# #   STEAM PINS
"""
    BasicSteamPin()
    Self computes T,s,x,V
    Must have methods for ṁ,Φ,P,h
""" 
@connector function BasicSteamPin(; name, Pdef=0.1)
    across_var  = @variables  P(t)=Pdef T(t)=300 s(t)=0.0 h(t)=191e3 x(t)=0.0 v(t)=.001
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
    eqs = [Q̇ ~ Qin]
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

@component function MultiPhaseGnd(;name, P = 0.1)
    @named n = BasicSteamPin(Pdef = P)
    ps = @parameters P=P 
    eqs = [
        n.P ~ P
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

@component function TwoPortReservoir(;name, P = 0.1)
    @named p = BasicSteamPin()
    @named n = BasicSteamPin(Pdef = P)
    ps = @parameters P=P
    eqs = [
        n.P ~ P
        p.P ~ P
        n.h  ~ stm_hsatfunc(P)
        p.h ~ n.h
        p.Φ ~ 0 #p.ṁ*p.h         # E = Pressure * Volumetric flow rate = PRessures * mass flow rate / density
        0 ~ n.ṁ + p.ṁ                                                           # 1/density = specific volume
        0 ~ n.Φ + p.Φ
    ]
    ODESystem(eqs, t,[], ps; name = name, systems = [n,p], defaults = [P => 0.1 ])
end

@component function ContinuityReservoir2(;name)
    @named p = BasicSteamPin()
    @named n = BasicSteamPin()
    eqs = [
        n.P ~ p.P
        n.h ~ p.h
        p.Φ ~ 0
        0 ~ n.Φ + p.Φ
        0 ~ n.ṁ + p.ṁ
    ]
    ODESystem(eqs, t,[], []; name = name, systems = [n,p])
end

@component function ContinuityReservoir(;name)
    @named p = BasicSteamPin()
    @named n = BasicSteamPin()
    sts = @variables ΔΦ(t)=0.0
    eqs = [
        n.P ~ p.P
        n.h ~ p.h
        n.Φ ~ 0
        ΔΦ ~ p.Φ
        0 ~ n.ṁ + p.ṁ
    ]
    ODESystem(eqs, t,sts, []; name = name, systems = [n,p], defaults = [ΔΦ => 0.1 ])
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

@component function TunableSteamFlowValve(;name, ṁ = 1.0)
    @named p = BasicSteamPin()
    @named n = BasicSteamPin()
    ps = @parameters ṁ = ṁ [tunable = true]
    eqs = [
            0 ~ p.Φ + n.Φ   # conservation of energy
            0 ~ p.ṁ + n.ṁ
            p.ṁ ~ ṁ
            p.h ~ n.h 
            n.P ~ p.P
        ]
    ODESystem(eqs, t, [], ps; name = name, systems = [n,p])
    # compose(ODESystem(eqs, t, [], ps; name = name, defaults = [Ṁ => 1.0]),p,n)
end

@component function AdiabaticPump(;name, η = 1.0, setpressure = false, Pout = 10, controlinlet = false)
    #Pout in bar
    @named p = BasicSteamPin()
    @named n = BasicSteamPin()
    @named w = WorkPin()
    
    ps = @parameters η = η P=Pout

    eqs = Equation[
        w.Ẇ ~ p.Φ + n.Φ                         # conservation of energy
        0 ~ p.ṁ + n.ṁ
        n.h  ~ p.h + p.v *1e5* (n.P - p.P) / η  # work, multiply by 100 to get to kPa
        0  ~ p.ṁ * (n.h - p.h)  + w.Ẇ 
    ]

    if setpressure
        eqs = vcat(eqs,n.P  ~ P)
    end

    if controlinlet
        eqs = vcat(eqs,p.h  ~ stm_hsatfunc(p.P))
    end
    # extenda([ODESystem(eqs, t,[], ps; name = name, defaults = [η => 0.6, P => Pout]), p,n,w])
    ODESystem(eqs, t,[], ps; name = name, systems = [p,n,w], defaults = [η => .6, P => 10])
end

@component function Splitter(;name)
    @named p = BasicSteamPin()
    @named n1 = BasicSteamPin()
    @named n2 = BasicSteamPin()

    #  0 ~ pMax₊P - trbsimo₊split₊z₊P(t)
    eqs = [
        p.h ~ n1.h               # enthalpy
        0 ~ n1.h - n2.h
        p.P ~ n1.P               # pressure
        n2.P ~ n1.P 
        0 ~ p.ṁ + n1.ṁ + n2.ṁ
        0 ~ p.Φ + n1.Φ + n2.Φ
    ]

    ODESystem(eqs, t,[], []; name = name, systems = [p,n1,n2])
end

@component function AdiabaticTurbine(;name, η = 1.0, setpressure = false, Pout = 0.1)
    #Pout in bar
    @named p = BasicSteamPin()
    @named n = BasicSteamPin()
    @named w = WorkPin()
    ps = @parameters P=Pout η = η

    sts = @variables Ẇ(t)=0

    eqs = Equation[
        Ẇ ~ p.Φ + n.Φ 
        w.Ẇ ~ Ẇ
        0 ~ p.ṁ + n.ṁ      
        n.h  ~ p.h - (p.h-stm_hpsfunc(n.P,p.s))*η
        0 ~ p.ṁ * (n.h - p.h) + w.Ẇ
    ]

    if setpressure
        eqs = vcat(eqs, n.P ~ P)
    end

    # extenda([ODESystem(eqs, t,[], ps; name = name, defaults = [η => 0.6, P => Pout]), p,n,w])
    ODESystem(eqs, t,sts, ps; name = name, systems = [p,n,w], defaults = [η => 1.0, P => Pout, Ẇ => 0])
end   

@component function SIMOAdiabaticTurbine(;name, ηin = 1.0, setpressure = false, Pyin = 10, Pzin  = 0.1)
    #Pout in bar
    ps = @parameters Py=Pyin Pz=Pzin η = ηin
    sp = setpressure

    # EXTERNAL NODES FOR INTERFACING
    @named p = BasicSteamPin()  # inlet node

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
        C ~ p.ṁ * stm_cphfunc(p.P, p.h) # * 1000   # duty T/h = cp
        0 ~ q.Q̇ - Q̇ 
        n.h ~ p.h + q.Q̇/(p.ṁ)
        n.P ~ p.P
    ]
    ODESystem(eqs,t,[Q̇,C],[]; name = name, systems = [p,n,q], defaults = [Q̇ => 0.0, C => 400])
end

@component function TunableSteamHeatTransfer(;name, Q̇in = 150e6)
    @named p = BasicSteamPin()
    @named n = BasicSteamPin()
    @named q = HeatTransferPin()
    ps = @parameters Q̇in = Q̇in [tunable = true]
    st = @variables C(t)=4000
    eqs = [
        0 ~ p.ṁ + n.ṁ           # p.ṁ ~ n.ṁ                               # conservation of mass
        q.Q̇ ~ p.Φ + n.Φ         # conservation of energy            
        C ~ p.ṁ * stm_cphfunc(p.P, p.h) # * 1000   # duty T/h = cp
        0 ~ q.Q̇ - Q̇in 
        n.h ~ p.h + q.Q̇/(p.ṁ)
        n.P ~ p.P
    ]
    ODESystem(eqs,t,[C],ps; name = name, systems = [p,n,q], defaults = [C => 400])
end

@component function IdealBoiler(;name, Tout = 350)
    #Pout in bar
    @named p = BasicSteamPin()
    @named n = BasicSteamPin()
    @named q = HeatTransferPin()
    
    ps = @parameters T = Tout

    eqs = Equation[
        q.Q̇ ~ p.Φ + n.Φ                 # conservation of energy
        0 ~ p.ṁ + n.ṁ
        n.P  ~  p.P                     # no pressure
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

@component function PassiveCondensor(;name)
    #Pout in bar
    @named p = BasicSteamPin()
    @named n = BasicSteamPin()
    @named q = HeatTransferPin()

    sts = @variables Q̇(t)=0.0
    
    eqs = Equation[
        q.Q̇ ~ p.ṁ * (n.h - p.h)
        0 ~ p.ṁ + n.ṁ
        q.Q̇ ~ Q̇
    ]
    ODESystem(eqs, t,sts, []; name = name, systems = [p,n,q], defaults = [Q̇ => 0.0])
end

@component function OpenFeedwaterHeater(;name)
    # flows x and y are the inlets
    @named p1 = BasicSteamPin()
    @named p2 = BasicSteamPin()
    @named n = BasicSteamPin()  

    sts=@variables yfrac(t) = 0.5 [bounds = (0.0,1.0)]


    continuous_events = [[yfrac ~ 1] => [yfrac ~ 1]
                        [yfrac ~ 0] => [yfrac ~ 0.0]]
    eqs  =[
        n.P ~ p2.P
        p1.P ~ n.P
        n.h ~ stm_hsatfunc(n.P) # sat liquid
        0 ~ n.h - (yfrac * p1.h + (1-yfrac) * p2.h)
        0 ~ n.ṁ*yfrac + p1.ṁ
        0 ~ n.ṁ + p1.ṁ + p2.ṁ
        0 ~ n.Φ + (p1.Φ + p2.Φ)
    ]
    ODESystem(eqs, t,sts, []; name = name, systems = [p1,p2,n], defaults = [yfrac => 0.5], continuous_events)
end

@component function MixingChamber(;name)
    # flows x and y are the inlets
    @named p1 = BasicSteamPin()
    @named p2 = BasicSteamPin()
    @named n = BasicSteamPin()  

    eqs  =[
        n.P ~ p1.P
        p1.P ~ p2.P
        n.h ~ 1/(p1.ṁ + p2.ṁ) * (p1.ṁ*p1.h + p2.ṁ*p2.h)
        0 ~ n.ṁ + p1.ṁ + p2.ṁ
        0 ~ n.Φ + p1.Φ + p2.Φ
    ]
    ODESystem(eqs, t,[], []; name = name, systems = [n,p1,p2])
end

@component function ClosedFeedwaterHeater(;name)

    # flows x and y are the inlets
    @named p1 = BasicSteamPin()
    @named n1 = BasicSteamPin()
    @named p2 = BasicSteamPin()  
    @named n2 = BasicSteamPin()
    # sts=@variables yfrac(t)=0.5 Q̇(t)=0.0
    #     #   Steam_OFW CLOSED FEEDWATER HEATER
    #     #                _______________
    #     #               |               |
    #     #           -> ___/\/\/\/\/\/\______ Steam_known = outlet
    #     #               |               |
    #     #               |_______________|
    #     # yfrac*p1.h + (1-yfrac)*p2.h ~ yfrac*n1.h + (1-yfrac)*n2.h
    #     # p1.ṁ ~ yfrac*(n1.ṁ + p1.ṁ)
    #     # p2.ṁ ~ (1-yfrac)*(n1.ṁ + p1.ṁ)
    sts = @variables  yfrac(t)=0.5 [bounds = (0.0, 1.0)]
    # continuous_events = [[yfrac ~ 1] => [yfrac ~ 1]
    #                     [yfrac ~ 0] => [yfrac ~ 0.0]]
    eqs  =[
        n1.P ~ p1.P # pressure
        n2.P ~ p2.P
        n2.h ~ n1.h # enthalpies
        0 ~ n1.ṁ + p1.ṁ
        0 ~ n2.ṁ + p2.ṁ
        yfrac ~ (n2.h - p2.h)/((p1.h-n1.h)+(n2.h-p2.h))
        p1.ṁ ~ yfrac * (p1.ṁ + p2.ṁ)
        p1.Φ ~ yfrac * (p1.Φ + p2.Φ)
        0 ~ n1.Φ + p1.Φ + p1.ṁ * p1.h + n1.ṁ  * n1.h
        0 ~ n2.Φ + p2.Φ + p2.ṁ * p2.h + n2.ṁ  * n2.h
    ]
    # 0 ~ n1.ṁ*n1.h+n2.ṁ*n2.h+p1.ṁ*p1.h+p2.ṁ*p2.h
    # n1.Φ + p1.Φ ~  n2.Φ  + p2.Φ
    # n1.Φ + p1.Φ ~  n2.Φ  + p2.Φ
    # p2.ṁ * (n2.h - p2.h) ~ p1.ṁ * (n1.h-p1.h)
    # 0 ~ n1.ṁ*n1.h+n2.ṁ*n2.h+p1.ṁ*p1.h+p2.ṁ*p2.h
    ODESystem(eqs, t,sts, []; name = name, systems = [p1,n1,p2,n2], defaults = [ yfrac => 0.5])
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
        n.Φ  ~ 0]         #conservation of energy
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

@component function FeedwaterRankine(; name,  Pmin = 0.1, Pmid = 10, Pmax = 150)
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
                            hydro_connect(valve.n, pumpB.p),          # pump -> boilder
                            hydro_connect(pumpB.n,gnd.p),
                            hydro_connect(gnd.n,boil.p),
                            hydro_connect(boil.n, turbn.p),             # boiler -> turbine
                            hydro_connect(turbn.hp.n, openfw.y),        # turbine -> openfw y
                            hydro_connect(turbn.lp.n, cndnsr.p),        #   '------> condensor
                            hydro_connect(cndnsr.n, pumpA.p),
                            hydro_connect(pumpA.n,openfw.z),                                
                            work_connect(WorkRes, turbn.lp.w,turbn.hp.w, pumpA.w, pumpB.w),
                            heat_connect(ColdUtil, cndnsr.q),
                            heat_connect(HotUtil, boil.q));
            
        systems = [valve, pumpA, pumpB, boil, turbn, cndnsr, openfw, gnd, WorkRes, ColdUtil, HotUtil]
        ODESystem(connections,t; name = name, systems = systems)
end

@component function TunableFeedwaterRankine(; name,  Pmin = 0.1, Pmid = 10, Pmax = 150)
    # Control elements
    @named gnd          = ContinuityReservoir()
    @named valve        = TunableSteamFlowValve(ṁ = 1.0 )
    @named pumpA        = AdiabaticPump(Pout = Pmid, setpressure = false)
    @named pumpB        = AdiabaticPump(Pout = Pmax, setpressure = true)

    #Boiler
    # @named boil         = IdealBoiler(Tout = 600+273)
    @named boil         = TunableSteamHeatTransfer(Q̇in = 150e6)
    @named turbn        = SIMOAdiabaticTurbine(setpressure = true,Pyin = Pmid, Pzin = Pmin, ηin = 1.0) 

    @named cndnsr       = IdealCondensor()
    @named openfw       = OpenFeedwaterHeater()

    @named WorkRes     = WorkPin()
    @named ColdUtil    = HeatTransferPin()
    @named HotUtil     = HeatTransferPin()

    connections = vcat(Steam.hydro_connect(openfw.n,valve.p),
                    hydro_connect(valve.n, pumpB.p),          # pump -> boilder
                    hydro_connect(pumpB.n,gnd.p),
                    hydro_connect(gnd.n,boil.p),
                    hydro_connect(boil.n, turbn.p),             # boiler -> turbine
                    hydro_connect(turbn.hp.n, openfw.p1),        # turbine -> openfw y
                    hydro_connect(turbn.lp.n, cndnsr.p),        #   '------> condensor
                    hydro_connect(cndnsr.n, pumpA.p),
                    hydro_connect(pumpA.n,openfw.p2),                                
                    work_connect(WorkRes, turbn.lp.w,turbn.hp.w, pumpA.w, pumpB.w),
                    heat_connect(ColdUtil, cndnsr.q),
                    heat_connect(HotUtil, boil.q));
    
    systems = [valve, pumpA, pumpB, boil, turbn, cndnsr, openfw, gnd, WorkRes, ColdUtil, HotUtil]
    ODESystem(connections,t; name = name, systems = systems)
end

@component function Water_loop(; name, Pmin=32, Pmax=40)
    @named WorkRes     = WorkPin()
    @named ColdUtil    = HeatTransferPin()
    @named res          = ioReservoir(P = Pmin, fixboth = true)
    @named valve        = SteamFlowValve()
    @named pump         = AdiabaticPump(η = 0.9, setpressure = true, Pout = Pmax)
    @named pset         = SetPressure(P=Pmin)
    @named HeatIn       = SteamHeatTransfer()
    @named HeatTx       = SteamHeatTransfer()
    @named HeatRej      = IdealCondensor()
    @named throttle     = throttle()

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

#     begin 
#     Plow = 0.1;
#     Pmid = 12;
#     Pmax = 150;
#     # Control elements
#     @named gnd          = ContinuityReservoir()
#     @named valve        = SteamFlowValve()
#     @named pumpA        = AdiabaticPump(Pout = Pmid, setpressure = false)
#     @named pumpB        = AdiabaticPump(Pout = Pmax, setpressure = true)
#     #Boiler
#     # @named boil         = IdealBoiler(Tout = 600+273)
#     @named boil         = SteamHeatTransfer()
#     @named turbn        = SIMOAdiabaticTurbine(setpressure = true,ηin = 1.0) 
#     @named cndnsr       = IdealCondensor()
#     @named openfw       = OpenFeedwaterHeater()
#     @named WorkRes     = WorkPin()
#     @named ColdUtil    = HeatTransferPin()
#     @named HotUtil     = HeatTransferPin()
#     connections = vcat(hydro_connect(openfw.n,valve.p),
#                         hydro_connect(valve.n, pumpB.p),          # pump -> boilder
#                         hydro_connect(pumpB.n,gnd.p),
#                         hydro_connect(gnd.n,boil.p),
#                         hydro_connect(boil.n, turbn.p),             # boiler -> turbine
#                         hydro_connect(turbn.hp.n, openfw.y),        # turbine -> openfw y
#                         hydro_connect(turbn.lp.n, cndnsr.p),        #   '------> condensor
#                         hydro_connect(cndnsr.n, pumpA.p),
#                         hydro_connect(pumpA.n,openfw.z),                                
#                         work_connect(WorkRes, turbn.lp.w,turbn.hp.w, pumpA.w, pumpB.w),
#                         heat_connect(ColdUtil, cndnsr.q),
#                         heat_connect(HotUtil, boil.q));
#     systems = [valve, pumpA, pumpB, boil, turbn, cndnsr, openfw, gnd, WorkRes, ColdUtil, HotUtil]
#     @named odesys = ODESystem(connections,t; systems = systems)
#     masfcn(t) = 100 + 15 * sin(t)
#     Qf(t) = 100e6 + 35e6 * cos(t)  
#     # @unpack valve, HeatIn, HeatTx, HeatRej,circulator = odesys
#     @unpack valve, boil,gnd, cndnsr,openfw, pumpA, WorkRes, turbn = odesys
#     controlled_eqs = [valve.p.ṁ ~ masfcn(t),
#                         boil.Q̇ ~ Qf(t)]
#     @named aux_sys = ODESystem(controlled_eqs,t)
#     @named odesys = extend(odesys, aux_sys)
#     odesys = substitute(odesys,hydro_prop_dict)
#     simple_sys = structural_simplify(odesys)
#     saveat = LinRange(1.0,9.0,25)
#     kwargs = (abstol=1e-10, reltol=1e-2, saveat = saveat)
#     tspan = (0.0,10.0)
#     states(simple_sys)
#     prob    = ODEProblem(simple_sys,[],tspan,)
#     sol   = solve(prob, Rodas4(); kwargs...);
#     # c   = [res,valve,comp1,reg.A,heat,turbine,reg.B,cool];
#     p1 = plot(sol, vars = [valve.p.ṁ], title = "Cycle Flow Rate")
#     p2 = plot(sol, vars = [boil.Q̇, cndnsr.Q̇], title = "Heat Loading")
#     p3 = plot(sol,vars = [WorkRes.Ẇ], title = "Net Work")
#     #     # p3 = plot(sol,vars = [pump.w.Ẇ])
#     p=plot(p1,p2,p3,layout = (3,1))
# end
# include("03-MTK_UTILS.jl")
# begin
# systems =  [valve, pumpA, pumpB, boil, turbn.hp, turbn.lp, cndnsr, openfw, gnd]
# showsol(systems,sol)
# end
# begin
#     Plow = 0.1;
#     Pmid = 12;
#     Pmax = 150;
#     @named gnd          = EGnd(P = Pmid)
#     @named pset         = SetPressure(P = Plow)
#     @named iores        = ioReservoir(P=Pmid, fixboth = false)
#     @named valve        = SteamFlowSource(ṁ= 1.0)
#     @named pumpB        = AdiabaticPump(Pout = Pmax, setpressure = true)
#     @named boil         = IdealBoiler(Tout = 600+273)
#     @named turbn        = SIMOAdiabaticTurbine(setpressure = false, Pyin = Pmid, Pzin = Plow, ηin = 1.0) 
#     @named pumpA        = AdiabaticPump(Pout = Pmid, setpressure = true)
#     @named cndnsr       = IdealCondensor()
#     @named openfw       = OpenFeedwaterHeater()
#     @named WorkRes      = WorkPin()
#     @named ColdUtil     = HeatTransferPin()
#     @named HotUtil      = HeatTransferPin()
#     # connections = vcat(hydro_connect(iores.n,valve.p),
#     #                     hydro_connect(valve.n, pumpB.p),
#     #                     hydro_connect(pumpB.n,boil.p),
#     #                     hydro_connect(boil.n, turbn.p),
#     #                     hydro_connect(turbn.hp.n, openfw.y),
#     #                     hydro_connect(turbn.lp.n, cndnsr.p),
#     #                     hydro_connect(cndnsr.p,pset.n),
#     #                     hydro_connect(pset.p,pumpA.p),
#     #                     hydro_connect(pumpA.n,openfw.z),
#     #                     hydro_connect(openfw.n,iores.p),
#     #                     work_connect(WorkRes, turbn.lp.w,turbn.hp.w, pumpA.w, pumpB.w),
#     #                     heat_connect(ColdUtil, cndnsr.q),
#     #                     heat_connect(HotUtil, boil.q));
#     # pumpB -> Valve -> boil
#     connections = vcat(hydro_connect(pumpB.n, boil.p),          # pump -> boilder
#                     hydro_connect(valve.n, pumpB.p, gnd.n),       # valve -> Pump (0 point)
#                     hydro_connect(boil.n, turbn.p),             # boiler -> turbine
#                     hydro_connect(turbn.hp.n, openfw.y),        # turbine -> openfw y
#                     hydro_connect(turbn.lp.n, cndnsr.p),        #   '------> condensor
#                     hydro_connect(cndnsr.n, pset.p),
#                     hydro_connect(pset.n, pumpA.p),
#                     hydro_connect(pumpA.n,openfw.z),
#                     hydro_connect(openfw.n,valve.p),
#                     work_connect(WorkRes, turbn.lp.w,turbn.hp.w, pumpA.w, pumpB.w),
#                     heat_connect(ColdUtil, cndnsr.q),
#                     heat_connect(HotUtil, boil.q));
#     systems = [gnd,pset,valve, pumpA, pumpB, boil, turbn, openfw,cndnsr, WorkRes, ColdUtil, HotUtil]
#     @named odesys = ODESystem(connections,t;systems = systems)
#     # @unpack valve, boil,iores,cndnsr, cndnsr,openfw, pumpA, WorkRes, turbn = odesys
#     @unpack valve, boil,cndnsr, cndnsr,openfw, pumpA, WorkRes, turbn = odesys
#     # masfcn(t) = 10.0
#     # controlled_eqs = [valve.p.ṁ ~ masfcn(t)]
#     #                     # boil.Q̇ ~ Qf(t)]
#     # @named aux_sys = ODESystem(controlled_eqs,t)
#     # @named odesys = extend(odesys, aux_sys)
#     odesys = substitute(odesys,hydro_prop_dict)
#     odesys = alias_elimination(odesys)
#     simple_sys = structural_simplify(alias_elimination(odesys))
#     saveat = LinRange(1.0,10,25)
#     kwargs = (abstol=1e-10, reltol=1e-2, saveat = saveat)
#     tspan = (0.0,10.0)
#     prob    = ODAEProblem(odesys,[],tspan,)
#     # sol     = solve(prob, Rodas4(); kwargs...);

#     sys = modelingtoolkitize(prob)

#     # # # c   = [res,valve,comp1,reg.A,heat,turbn,reg.B,cool];
#     # p1 = plot(sol, vars = [valve.p.ṁ], title = "Cycle Flow Rate")
#     # p2 = plot(sol, vars = [boil.q.Q̇], title = "Heat Loading")
#     # p3 = plot(sol,vars = [WorkRes.Ẇ], title = "Net Work")
#     # p4 = plot(sol, vars = [cndnsr.Q̇], title = "condensor")
#     # #     # p3 = plot(sol,vars = [pump.w.Ẇ])
#     # p=plot(p1,p2,p3,p4,layout = (4,1))end
# end
# include("03-MTK_UTILS.jl")
# begin
# systems = [valve,pumpB,boil,turbn,turbn.hp,turbn.lp,openfw,cndnsr,pumpA, WorkRes, ColdUtil, HotUtil];
# showsol(systems,sol)
# end
# x_ph(150,1.67e3)
# sol[cndnsr.p.h]
# sol[pumpA.n.Φ]
# sol[pumpA.p.Φ]
# sol[pumpA.w.Ẇ]
# sol[turbn.hp.w.Ẇ]
# sol[openfw.yfrac]
# sol[turbn.p.Φ]
# sol[turbn.hp.n.Φ]
# sol[turbn.lp.n.Φ]
# sol[turbn.hp.n.Φ][end]
# sol[turbn.lp.n.Φ][end]
# sol[openfw.z.Φ]
# sol[boil.n.h]
# sol[boil.p.h]
# sol[boil.Q̇]
# sol[boil.p.h]
# sol[boil.Q̇] .- sol[boil.p.ṁ] .* (sol[boil.n.h] .- sol[boil.p.h])
# x_ph(150,2.09e3)
# boil.p.h
# sol[cndnsr.Q̇]
# sol[cndnsr.n.T]
# sol[boil.n.x]
# sol[turbn.lp.n.P]
# h_ps(.1,2.105)
# iores.n,pumpA.p,pumpA.n
# # connections = vcat(hydro_connect(pumpB.n,valve.p),
# #                     hydro_connect(valve.n, boil.p),
#                     hydro_connect(boil.n, turbine.p),
#                     hydro_connect(turbine.hp.n, openfw.y),
#                     hydro_connect(turbine.lp.n, condensor.p),
#                     hydro_connect(condensor.n,iores.p),
#                     hydro_connect(iores.n,pumpA.p),
#                     hydro_connect(pumpA.n,openfw.z),
#                     hydro_connect(openfw.n,pumpB.p),
#                     work_connect(WorkRes, turbine.lp.w,turbine.hp.w, pumpA.w, pumpB.w),
#                     heat_connect(ColdUtil, condensor.q),
#                     heat_connect(HotUtil, boil.q));