# module ThermoSteam
using ModelingToolkit, Revise, Unitful, Symbolics, Logging, Printf, XSteamTP # Unitful
# using DifferentialEquations

@variables t
# Logging.disable_logging(Logging.Warn)
# include("03-MTK_UTILS.jl")
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

stm_cpfunc(x, y) = Cp_ph(x, y / 1e3) * 1e3
stm_cphfunc(x, y) = Cp_ph(x, y / 1e3) * 1e3
stm_vphfunc(x, y) = v_ph(x, y / 1e3)
stm_Tpsfunc(x, y) = T_ps(x, y / 1e3) + 273.15
stm_Tphfunc(x, y) = T_ph(x, y / 1e3) + 273.15
stm_xphfunc(x, y) = x_ph(x, y / 1e3)
stm_xptfunc(x, y) = x_pT(x, y - 273.15)
stm_sptfunc(x, y) = s_pT(x, y - 273.15) * 1e3
stm_sphfunc(x, y) = s_ph(x, y / 1e3) * 1e3
stm_hptfunc(x, y) = h_pT(x, y - 273.15) * 1e3
stm_hpsfunc(x, y) = h_ps(x, y / 1e3) * 1e3
stm_hsatfunc(x)   = hL_p(x) * 1e3


@variables x, y
@register_symbolic stm_cpfunc(x, y)
@register_symbolic stm_cphfunc(x, y)
@register_symbolic stm_vphfunc(x, y)
@register_symbolic stm_Tpsfunc(x, y)
@register_symbolic stm_Tphfunc(x, y)
@register_symbolic stm_xphfunc(x, y)
@register_symbolic stm_xptfunc(x, y)
@register_symbolic stm_sptfunc(x, y)
@register_symbolic stm_sphfunc(x, y)
@register_symbolic stm_hptfunc(x, y)
@register_symbolic stm_hpsfunc(x, y)
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


hydro_prop_dict = Dict(
    stm_Tpsfunc => T_ps,
    stm_Tphfunc => T_ph,
    stm_xphfunc => x_ph,
    stm_sptfunc => s_pT,
    stm_sphfunc => s_ph,
    stm_hptfunc => h_pT,
    stm_hpsfunc => h_ps,
    stm_hsatfunc => hL_p,
    stm_vphfunc => v_ph,
    stm_cphfunc => Cp_ph,
);



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
    BasicSteamPin(; name, Pdef = 0.1)

DOCSTRING
    Basic Steam Pin
    Self computes T,s,x,V
    Must have methods for ṁ,Φ,P,h
    across_var =
    @variables P(t) = Pdef T(t) = 300 s(t) = 0.0 h(t) = 191e3 x(t) = 0.0 v(t) = 0.001
    thru_var = @variables ṁ(t) = 0.0 Φ(t) = 0.0                     # mass flow and energy flow
ELEM TYPE: PIN
EQUATIONS:
    T ~ stm_Tphfunc(P, h)
    s ~ stm_sphfunc(P, h)
    x ~ stm_xphfunc(P, h)
    v ~ stm_vphfunc(P, h)
INPUTS
- name: Name of the system, symbol. Or use the @named macro when initializing.
- Pdef = 0.1: DESCRIPTION
"""
function BasicSteamPin(; name, Pdef = 0.1)
    across_var =
        @variables P(t) = Pdef T(t) = 300 s(t) = 0.0 h(t) = 191e3 x(t) = 0.0 v(t) = 0.001
    thru_var = @variables ṁ(t) = 0.0 Φ(t) = 0.0                     # mass flow and energy flow
    ps = []
    eqs = Equation[
        T ~ stm_Tphfunc(P, h)
        s ~ stm_sphfunc(P, h)
        x ~ stm_xphfunc(P, h)
        v ~ stm_vphfunc(P, h)
    ]
    ODESystem(
        eqs,
        t,
        [thru_var..., across_var...],
        ps;
        name = name,
        defaults = [P => Pdef, ṁ => 0, Φ => 0, h => 191],
    )
end

function WorkPin(; name)
    sts = @variables Ẇ(t) = 0.0
    ODESystem(Equation[], t, sts, []; name = name)
end

function HeatTransferPin(; name)  #input in W
    sts = @variables Q̇(t) = 0.0
    ODESystem(Equation[], t, sts, []; name = name, defaults = [Q̇ => 0.0])
end

function FixedHeatFlowPin(; name, Qin = 1e3)  #input in W
    sts = @variables Q̇(t) = Qin
    ps = @parameters Qin = Qin
    eqs = [Q̇ ~ Qin]
    ODESystem(eqs, t, sts, ps; name = name, defaults = [Q̇ => Qin])
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

function hydro_basic_connect(n, p)
    eqs = [
        0 ~ p.ṁ + n.ṁ
        0 ~ p.Φ + n.Φ
        p.h ~ n.h
        p.P ~ n.P
    ]
    return eqs
end

function hydro_series_connect(comp; returnmode = :eq)
    eqs = hydro_basic_connect(comp[1].n, comp[2].p)   # n -> p

    if length(comp) > 2
        for i = 2:length(comp)-1
            eqs = vcat(eqs, hydro_basic_connect(comp[i].n, comp[i+1].p))
        end
    end

    if returnmode == :ode
        return ODESystem(eqs, t; name = :connections)
    end

    return eqs
end

function extenda(odevec::Vector{ODESystem})
    if length(odevec) == 1
        return odevec[1]
    else
        return extend(odevec[1], extenda(odevec[2:end]))
    end
end

function HeatTransferPort(; name)
    @named p = HeatTransferPin()
    @named n = HeatTransferPin()
    eqs = [0 ~p.Q̇ + n.Q̇]
    @variables Q̇(t)
    ODESystem(Equation[], t, [Q̇], []; name = name, systems = [p, n])
end

"""
    Reservoir(; name, P = 0.1)

DOCSTRING


ELEM TYPE: UTILITY
EQUATIONS:
n.P ~ P
n.h ~ stm_hsatfunc(P)
n.Φ ~ 0
INPUTS
- name: Name of the system, symbol. Or use the @named macro when initializing.
-  = 0.1`: DESCRIPTION
"""
function Reservoir(; name, P = 0.1)
    @named n = BasicSteamPin(Pdef = P)
    ps = @parameters P = P
    eqs = [
        n.P ~ P
        n.h ~ stm_hsatfunc(P)
        n.Φ ~ 0
    ]
    ODESystem(eqs, t, [], ps; name = name, systems = [n])
end

"""
    MultiPhaseGnd(; name, P = 0.1)

DOCSTRING


ELEM TYPE: UTILITY
EQUATIONS:
n.P ~ P
n.Φ ~ 0
INPUTS
- name: Name of the system, symbol. Or use the @named macro when initializing.
-  = 0.1`: DESCRIPTION
"""
function MultiPhaseGnd(; name, P = 0.1)
    @named n = BasicSteamPin(Pdef = P)
    ps = @parameters P = P
    eqs = [
        n.P ~ P
        n.Φ ~ 0
    ]
    ODESystem(eqs, t, [], ps; name = name, systems = [n])
end

"""
    SuperHeatedReservoir(; name, P = 150, T = 600)

DOCSTRING


ELEM TYPE: UTILITY
EQUATIONS:
n.P ~ P
0 ~ stm_hptfunc(P, T) + n.h
n.Φ ~ -n.ṁ * n.h
INPUTS
- name: Name of the system, symbol. Or use the @named macro when initializing.
- `P = 150`: DESCRIPTION
- `T = 600`: DESCRIPTION
"""
function SuperHeatedReservoir(; name, P = 150, T = 600)
    @named n = BasicSteamPin(Pdef = P)
    ps = @parameters P = P T = T
    eqs = [
        n.P ~ P
        0 ~ stm_hptfunc(P, T) + n.h
        n.Φ ~ -n.ṁ * n.h
    ]
    ODESystem(eqs, t, [], ps; name = name, systems = [n])
end

"""
    ioReservoir(; name, P = 0.1, fixboth = false)

DOCSTRING
Fluid reference element with the option to set both the inlet and outlet pressures

ELEM TYPE: UTILITY
EQUATIONS:
n.P ~ P
n.h ~ stm_hsatfunc(P)
n.Φ ~ 0
0 ~ n.ṁ + p.ṁ
INPUTS
- name: Name of the system, symbol. Or use the @named macro when initializing.
- P = 0.1: Pressure (bar)
- `fixboth = false`: Option to fix both the inlet an d the outlet
"""
function ioReservoir(; name, P = 0.1, fixboth = false)
    @named p = BasicSteamPin()
    @named n = BasicSteamPin(Pdef = P)
    ps = @parameters P = P
    eqs = [
        n.P ~ P
        n.h ~ stm_hsatfunc(P)
        n.Φ ~ 0
        0 ~ n.ṁ + p.ṁ
    ]
    if fixboth
        push!(eqs, p.P ~ P)
        # push!(eqs,p.h ~ n.h)
    end
    ODESystem(eqs, t, [], ps; name = name, systems = [n, p])
end

"""
    TwoPortReservoir(; name, P = 0.1)

DOCSTRING
Basic
ELEM TYPE: UTILITY
EQUATIONS:
n.P ~ P
p.P ~ P
n.h ~ stm_hsatfunc(P)
p.h ~ n.h
p.Φ ~ 0             
0 ~ n.ṁ + p.ṁ                                                           
0 ~ n.Φ + p.Φ
INPUTS
- name: Name of the system, symbol. Or use the @named macro when initializing.
- P = 0.1: Pressure for both the p and n ports
"""
function TwoPortReservoir(; name, P = 0.1)
    @named p = BasicSteamPin()
    @named n = BasicSteamPin(Pdef = P)
    ps = @parameters P = P
    eqs = [
        n.P ~ P
        p.P ~ P
        n.h ~ stm_hsatfunc(P)
        p.h ~ n.h
        p.Φ ~ 0             #p.ṁ*p.h         # E = Pressure * Volumetric flow rate = PRessures * mass flow rate / density
        0 ~ n.ṁ + p.ṁ                                                           # 1/density = specific volume
        0 ~ n.Φ + p.Φ
    ]
    ODESystem(eqs, t, [], ps; name = name, systems = [n, p])
end

"""
    ContinuityReservoir2(; name)

DOCSTRING


ELEM TYPE: UTILITY
EQUATIONS:
n.P ~ p.P
n.h ~ p.h
p.Φ ~ 0
0 ~ n.Φ + p.Φ
0 ~ n.ṁ + p.ṁ
INPUTS
- name: Name of the system, symbol. Or use the @named macro when initializing.
"""
function ContinuityReservoir2(; name)
    @named p = BasicSteamPin()
    @named n = BasicSteamPin()
    eqs = [
        n.P ~ p.P
        n.h ~ p.h
        p.Φ ~ 0
        0 ~ n.Φ + p.Φ
        0 ~ n.ṁ + p.ṁ
    ]
    ODESystem(eqs, t, [], []; name = name, systems = [n, p])
end

"""
    ContinuityReservoir(; name)

DOCSTRING


ELEM TYPE: UTILITY
EQUATIONS:
n.P ~ p.P
n.h ~ p.h
n.Φ ~ 0
ΔΦ ~ p.Φ
0 ~ n.ṁ + p.ṁ
INPUTS
- name: Name of the system, symbol. Or use the @named macro when initializing.
"""
function ContinuityReservoir(; name)
    @named p = BasicSteamPin()
    @named n = BasicSteamPin()
    sts = @variables ΔΦ(t) = 0.0
    eqs = [
        n.P ~ p.P
        n.h ~ p.h
        n.Φ ~ 0
        ΔΦ ~ p.Φ
        0 ~ n.ṁ + p.ṁ
    ]
    ODESystem(eqs, t, sts, []; name = name, systems = [n, p], defaults = [ΔΦ => 0.1])
end

"""
    SetPressure(; name, P = 0.1)

DOCSTRING
Ideal pressure source

ELEM TYPE: SOURCE
EQUATIONS:
p.P ~ P
n.P ~ p.P
n.h ~ p.h
0 ~ p.Φ + n.Φ           
0 ~ p.ṁ + n.ṁ
INPUTS
- name: Name of the system, symbol. Or use the @named macro when initializing.
- P = 0.1: DESCRIPTION
"""
function SetPressure(; name, P = 0.1)
    @named p = BasicSteamPin(Pdef = P)
    @named n = BasicSteamPin(Pdef = P)
    ps = @parameters P = P
    eqs = [
        p.P ~ P
        n.P ~ p.P
        n.h ~ p.h
        0 ~ p.Φ + n.Φ             # conservation of energy
        0 ~ p.ṁ + n.ṁ
    ]
    ODESystem(eqs, t, [], ps; name = name, systems = [p, n], defaults = [P => 0.1])
end

"""
    SteamFlowSource(; name, ṁ = 1.0)

DOCSTRING


ELEM TYPE: SOURCE
EQUATIONS:
0 ~ p.Φ + n.Φ 
0 ~ p.ṁ + n.ṁ
p.ṁ ~ Ṁ
p.h ~ n.h
n.P ~ p.P
INPUTS
- name: Name of the system, symbol. Or use the @named macro when initializing.
- `ṁ = 1.0`: Mass flow rate (kg/s)
"""
function SteamFlowSource(; name, ṁ = 1.0)
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
    ODESystem(eqs, t, [], ps; name = name, systems = [n, p], defaults = [Ṁ => ṁ])
    # compose(ODESystem(eqs, t, [], ps; name = name, defaults = [Ṁ => 1.0]),p,n)
end

"""
    SteamFlowValve(; name)

DOCSTRING


ELEM TYPE: UTILITY
EQUATIONS:
    0 ~ p.Φ + n.Φ  
    0 ~ p.ṁ + n.ṁ
    p.ṁ ~ Ṁ
    p.h ~ n.h
    n.P ~ p.P
INPUTS
- name: Name of the system, symbol. Or use the @named macro when initializing.
"""
function SteamFlowValve(; name)
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
    ODESystem(eqs, t, sts, []; name = name, systems = [n, p], defaults = [Ṁ => 1.0])
    # compose(ODESystem(eqs, t, [], ps; name = name, defaults = [Ṁ => 1.0]),p,n)
end

"""
    TunableSteamFlowValve(; name, ṁ = 1.0)

DOCSTRING


ELEM TYPE: UTILITY
EQUATIONS:
    0 ~ p.Φ + n.Φ  
    0 ~ p.ṁ + n.ṁ
    p.ṁ ~ ṁ
    p.h ~ n.h
    n.P ~ p.P
INPUTS
- name: Name of the system, symbol. Or use the @named macro when initializing.
- `ṁ = 1.0`: Mass flow rate (kg/s)
"""
function TunableSteamFlowValve(; name, ṁ = 1.0)
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
    ODESystem(eqs, t, [], ps; name = name, systems = [n, p])
    # compose(ODESystem(eqs, t, [], ps; name = name, defaults = [Ṁ => 1.0]),p,n)
end

"""
    AdiabaticPump(; name, η = 1.0, setpressure = false, Pout = 10, controlinlet = false)

DOCSTRING
Adiabatic pump
# work, multiply by 100 to get to kPa
@named p = BasicSteamPin()
@named n = BasicSteamPin()
@named w = WorkPin()
ps = @parameters η = η P = Pout

ELEM TYPE: COMPONENT
EQUATIONS:
w.Ẇ ~ p.Φ + n.Φ                       
0 ~ p.ṁ + n.ṁ
n.h ~ p.h + p.v * 1e5 * (n.P - p.P) / η 
w.Ẇ ~ p.ṁ * (n.h - p.h)
INPUTS
- name: Name of the system, symbol. Or use the @named macro when initializing.
- `η = 1.0`: Isentropic Effciency
- `setpressure = false`: Option to constrain outlet pressure
- `Pout = 10`: Pressure (bar)
- `controlinlet = false`: DESCRIPTION
"""
function AdiabaticPump(;
    name,
    η = 1.0,
    setpressure = false,
    Pout = 10,
    controlinlet = false,
)
    #Pout in bar
    @named p = BasicSteamPin()
    @named n = BasicSteamPin()
    @named w = WorkPin()
    ps = @parameters η = η P = Pout

    eqs = Equation[
        w.Ẇ ~ p.Φ + n.Φ                         # conservation of energy
        0 ~ p.ṁ + n.ṁ
        n.h ~ p.h + p.v * 1e5 * (n.P - p.P) / η  # work, multiply by 100 to get to kPa
        w.Ẇ ~ p.ṁ * (n.h - p.h)
    ]

    if setpressure
        eqs = vcat(eqs, n.P ~ P)
    end

    if controlinlet
        eqs = vcat(eqs, p.h ~ stm_hsatfunc(p.P))
    end
    # extenda([ODESystem(eqs, t,[], ps; name = name, defaults = [η => 0.6, P => Pout]), p,n,w])
    ODESystem(eqs, t, [], ps; name = name, systems = [p, n, w])
end

"""
    Splitter(; name)

DOCSTRING


ELEM TYPE: UTILITY
EQUATIONS:

INPUTS
- name: Name of the system, symbol. Or use the @named macro when initializing.
"""
function Splitter(; name)
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

    ODESystem(eqs, t, [], []; name = name, systems = [p, n1, n2])
end

"""
    AdiabaticTurbine(; name, η = 1.0, setpressure = false, Pout = 0.1)

DOCSTRING


ELEM TYPE: COMPONENT
EQUATIONS:
Ẇ ~ p.Φ + n.Φ
w.Ẇ ~ Ẇ
0 ~ p.ṁ + n.ṁ
n.h ~ p.h - (p.h - stm_hpsfunc(n.P, p.s)) * η
w.Ẇ ~ p.ṁ * (n.h - p.h)
INPUTS
- name: Name of the system, symbol. Or use the @named macro when initializing.
- `η = 1.0`: Isentropic Effciency
- `setpressure = false`: Option to constrain outlet pressure
- `Pout = 0.1`: DESCRIPTION
"""
function AdiabaticTurbine(; name, η = 1.0, setpressure = false, Pout = 0.1)
    #Pout in bar
    @named p = BasicSteamPin()
    @named n = BasicSteamPin()
    @named w = WorkPin()
    ps = @parameters P = Pout η = η

    sts = @variables Ẇ(t) = 0

    eqs = Equation[
        Ẇ ~ p.Φ + n.Φ
        w.Ẇ ~ Ẇ
        0 ~ p.ṁ + n.ṁ
        n.h ~ p.h - (p.h - stm_hpsfunc(n.P, p.s)) * η
        w.Ẇ ~ p.ṁ * (n.h - p.h)
    ]

    if setpressure
        eqs = vcat(eqs, n.P ~ P)
    end

    # extenda([ODESystem(eqs, t,[], ps; name = name, defaults = [η => 0.6, P => Pout]), p,n,w])
    ODESystem(
        eqs,
        t,
        sts,
        ps;
        name = name,
        systems = [p, n, w],
        defaults = [P => Pout, Ẇ => 0],
    )
end

"""
    SIMOAdiabaticTurbine(; name, ηin = 1.0, setpressure = false, Pyin = 10, Pzin = 0.1)

DOCSTRING
Single input multi output turbinbe
EXTERNAL NODES FOR INTERFACING
@named p = BasicSteamPin()  # inlet node
@named hp = AdiabaticTurbine(η = ηin, Pout = Pyin, setpressure = sp)
@named lp = AdiabaticTurbine(η = ηin, Pout = Pzin, setpressure = sp)
@named yn = BasicSteamPin()
@named zn = BasicSteamPin()
split_connect = hydro_connect(p, yn, zn)    # mflow, p.ṁ => positive
hp_connect = hydro_connect(yn, hp.p)
lp_connect = hydro_connect(zn, lp.p)
#                yn(-) -------> y
# IN --> p (+) --|
#                zn(-) -------> z

ELEM TYPE: COMPONENT
EQUATIONS:

INPUTS
- name: Name of the system, symbol. Or use the @named macro when initializing.
- `ηin = 1.0`: DESCRIPTION
- `setpressure = false`: Option to constrain outlet pressure
- `Pyin = 10`: DESCRIPTION
- `Pzin = 0.1`: DESCRIPTION
"""
function SIMOAdiabaticTurbine(;
    name,
    ηin = 1.0,
    setpressure = false,
    Pyin = 10,
    Pzin = 0.1,
)
    #Pout in bar
    ps = @parameters Py = Pyin Pz = Pzin η = ηin
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
    hp_connect = hydro_connect(yn, hp.p)
    lp_connect = hydro_connect(zn, lp.p)


    eqs = vcat(split_connect, hp_connect, lp_connect)

    compose(ODESystem(eqs, t, [], ps; name = name, systems = [yn, zn, p]), hp, lp)
    # ODESystem(eqs, t,[], ps; name = name, systems = [p,hp,lp], defaults = [η => 1.0, Py => 10, Pz => 0.1])
    # extend(ODESystem(eqs, t,sts, ps; name = name, systems = [hp,lp], defaults = [η => 1.0 Py => 10 Pz => 0.1]),split)
end

"""
    SteamHeatTransfer(; name)

DOCSTRING


ELEM TYPE: COMPONENT
EQUATIONS:
    0 ~ p.ṁ + n.ṁ           
    q.Q̇ ~ p.Φ + n.Φ             
    C ~ p.ṁ * stm_cphfunc(p.P, p.h) # * 1000
    0 ~ q.Q̇ - Q̇
    n.h ~ p.h + q.Q̇ / (p.ṁ)
    n.P ~ p.P
INPUTS
- name: Name of the system, symbol. Or use the @named macro when initializing.
"""
function SteamHeatTransfer(; name)
    @named p = BasicSteamPin()
    @named n = BasicSteamPin()
    @named q = HeatTransferPin()
    st = @variables Q̇(t) = 1.0 C(t) = 4000
    eqs = [
        0 ~ p.ṁ + n.ṁ           # p.ṁ ~ n.ṁ                               # conservation of mass
        q.Q̇ ~ p.Φ + n.Φ         # conservation of energy            
        C ~ p.ṁ * stm_cphfunc(p.P, p.h) # * 1000   # duty T/h = cp
        0 ~ q.Q̇ - Q̇
        n.h ~ p.h + q.Q̇ / (p.ṁ)
        n.P ~ p.P
    ]
    ODESystem(eqs, t, [Q̇, C], []; name = name, systems = [p, n, q])
end

"""
    TunableSteamHeatTransfer(; name, Q̇in = 1.5e8)

DOCSTRING


ELEM TYPE: COMPONENT
EQUATIONS:
0 ~ p.ṁ + n.ṁ           # p.ṁ ~ n.ṁ                               # conservation of mass
q.Q̇ ~ p.Φ + n.Φ         # conservation of energy            
C ~ p.ṁ * stm_cphfunc(p.P, p.h) # * 1000   # duty T/h = cp
0 ~ q.Q̇ - Q̇in
n.h ~ p.h + q.Q̇ / (p.ṁ)
n.P ~ p.P
INPUTS
- name: Name of the system, symbol. Or use the @named macro when initializing.
- `Q̇in = 1.5e8`: DESCRIPTION
"""
function TunableSteamHeatTransfer(; name, Q̇in = 150e6)
    @named p = BasicSteamPin()
    @named n = BasicSteamPin()
    @named q = HeatTransferPin()
    ps = @parameters Q̇in = Q̇in [tunable = true]
    st = @variables C(t) = 4000
    eqs = [
        0 ~ p.ṁ + n.ṁ           # p.ṁ ~ n.ṁ                               # conservation of mass
        q.Q̇ ~ p.Φ + n.Φ         # conservation of energy            
        C ~ p.ṁ * stm_cphfunc(p.P, p.h) # * 1000   # duty T/h = cp
        0 ~ q.Q̇ - Q̇in
        n.h ~ p.h + q.Q̇ / (p.ṁ)
        n.P ~ p.P
    ]
    ODESystem(eqs, t, [C], ps; name = name, systems = [p, n, q], defaults = [C => 400])
end

"""
    IdealBoiler(; name, Tout = 350)

DOCSTRING


ELEM TYPE: COMPONENT
EQUATIONS:
q.Q̇ ~ p.Φ + n.Φ                 # conservation of energy
0 ~ p.ṁ + n.ṁ
n.P ~ p.P                     # no pressure
n.h ~ stm_hptfunc(p.P, T)       # work, multiply by 100 to get to kPa
q.Q̇ ~ p.ṁ * (n.h - p.h)
INPUTS
- name: Name of the system, symbol. Or use the @named macro when initializing.
- `Tout = 350`: DESCRIPTION
"""
function IdealBoiler(; name, Tout = 350)
    #Pout in bar
    @named p = BasicSteamPin()
    @named n = BasicSteamPin()
    @named q = HeatTransferPin()

    ps = @parameters T = Tout

    eqs = Equation[
        q.Q̇ ~ p.Φ + n.Φ                 # conservation of energy
        0 ~ p.ṁ + n.ṁ
        n.P ~ p.P                     # no pressure
        n.h ~ stm_hptfunc(p.P, T)       # work, multiply by 100 to get to kPa
        q.Q̇ ~ p.ṁ * (n.h - p.h)
    ]
    ODESystem(eqs, t, [], ps; name = name, systems = [p, n, q], defaults = [T => Tout])
end

"""
    IdealCondensor(; name)

DOCSTRING


ELEM TYPE: COMPONENT
EQUATIONS:
q.Q̇ ~ p.Φ + n.Φ             # conservation of energy
0 ~ p.ṁ + n.ṁ
n.P ~ p.P                 # no pressure
n.h ~ stm_hsatfunc(p.P)        # work, multiply by 100 to get to kPa
q.Q̇ ~ p.ṁ * (n.h - p.h)
q.Q̇ ~ Q̇
INPUTS
- name: Name of the system, symbol. Or use the @named macro when initializing.
"""
function IdealCondensor(; name)
    #Pout in bar
    @named p = BasicSteamPin()
    @named n = BasicSteamPin()
    @named q = HeatTransferPin()

    sts = @variables Q̇(t) = 0.0

    eqs = Equation[
        q.Q̇ ~ p.Φ + n.Φ             # conservation of energy
        0 ~ p.ṁ + n.ṁ
        n.P ~ p.P                 # no pressure
        n.h ~ stm_hsatfunc(p.P)        # work, multiply by 100 to get to kPa
        q.Q̇ ~ p.ṁ * (n.h - p.h)
        q.Q̇ ~ Q̇
    ]
    ODESystem(eqs, t, sts, []; name = name, systems = [p, n, q], defaults = [Q̇ => 0.0])
end

"""
    PassiveCondensor(; name)

DOCSTRING


ELEM TYPE: COMPONENT
EQUATIONS:
q.Q̇ ~ p.ṁ * (n.h - p.h)
0 ~ p.ṁ + n.ṁ
q.Q̇ ~ Q̇
INPUTS
- name: Name of the system, symbol. Or use the @named macro when initializing.
"""
function PassiveCondensor(; name)
    #Pout in bar
    @named p = BasicSteamPin()
    @named n = BasicSteamPin()
    @named q = HeatTransferPin()

    sts = @variables Q̇(t) = 0.0

    eqs = Equation[
        q.Q̇ ~ p.ṁ * (n.h - p.h)
        0 ~ p.ṁ + n.ṁ
        q.Q̇ ~ Q̇
    ]
    ODESystem(eqs, t, sts, []; name = name, systems = [p, n, q], defaults = [Q̇ => 0.0])
end

"""
    ReliefElement(; name, pressurecontrol = false)

DOCSTRING


ELEM TYPE: COMPONENT
EQUATIONS:
n.ṁ + p.ṁ ~ 0         
q.Q̇ ~ p.Φ + n.Φ         
INPUTS
- name: Name of the system, symbol. Or use the @named macro when initializing.
- `pressurecontrol = false`: DESCRIPTION
"""
function ReliefElement(; name, pressurecontrol = false)
    @named p = BasicSteamPin()
    @named n = BasicSteamPin()
    @named q = HeatTransferPin()
    # 0 variables, 3 equations for pin
    eqs = [
        n.ṁ + p.ṁ ~ 0           # has to be negative
        q.Q̇ ~ p.Φ + n.Φ         # conservation of energy
    ]
    if pressurecontrol
        eqs = vcat(eqs, n.P ~ p.P)
    end
    ODESystem(eqs, t, [], []; name = name, systems = [p, n, q])
end

"""
    OpenFeedwaterHeater(; name)

DOCSTRING


ELEM TYPE: COMPONENT
EQUATIONS:
[yfrac ~ 1] => [yfrac ~ 1]
[yfrac ~ 0] => [yfrac ~ 0.0]
n.P ~ p2.P
p1.P ~ n.P
n.h ~ stm_hsatfunc(n.P) # sat liquid
0 ~ n.h - (yfrac * p1.h + (1 - yfrac) * p2.h)
0 ~ n.ṁ * yfrac + p1.ṁ
0 ~ n.ṁ + p1.ṁ + p2.ṁ
0 ~ n.Φ + (p1.Φ + p2.Φ)
INPUTS
- name: Name of the system, symbol. Or use the @named macro when initializing.
"""
function OpenFeedwaterHeater(; name)
    # flows x and y are the inlets
    @named p1 = BasicSteamPin()
    @named p2 = BasicSteamPin()
    @named n = BasicSteamPin()

    sts = @variables yfrac(t) = 0.5 #S[bounds = (0.0,1.0)]


    continuous_events = [
        [yfrac ~ 1] => [yfrac ~ 1]
        [yfrac ~ 0] => [yfrac ~ 0.0]
    ]
    eqs = [
        n.P ~ p2.P
        p1.P ~ n.P
        n.h ~ stm_hsatfunc(n.P) # sat liquid
        0 ~ n.h - (yfrac * p1.h + (1 - yfrac) * p2.h)
        0 ~ n.ṁ * yfrac + p1.ṁ
        0 ~ n.ṁ + p1.ṁ + p2.ṁ
        0 ~ n.Φ + (p1.Φ + p2.Φ)
    ]
    ODESystem(eqs, t, sts, []; name = name, systems = [p1, p2, n])
end

"""
    MixingChamber(; name)

DOCSTRING
    @named p1 = BasicSteamPin()
    @named p2 = BasicSteamPin()
    @named n = BasicSteamPin()

ELEM TYPE: COMPONENT
EQUATIONS:
n.P ~ p1.P
p1.P ~ p2.P
n.h ~ 1 / (p1.ṁ + p2.ṁ) * (p1.ṁ * p1.h + p2.ṁ * p2.h)
0 ~ n.ṁ + p1.ṁ + p2.ṁ
0 ~ n.Φ + p1.Φ + p2.Φ
INPUTS
- name: Name of the system, symbol. Or use the @named macro when initializing.
"""
function MixingChamber(; name)
    # flows x and y are the inlets
    @named p1 = BasicSteamPin()
    @named p2 = BasicSteamPin()
    @named n = BasicSteamPin()

    eqs = [
        n.P ~ p1.P
        p1.P ~ p2.P
        n.h ~ 1 / (p1.ṁ + p2.ṁ) * (p1.ṁ * p1.h + p2.ṁ * p2.h)
        0 ~ n.ṁ + p1.ṁ + p2.ṁ
        0 ~ n.Φ + p1.Φ + p2.Φ
    ]
    ODESystem(eqs, t, [], []; name = name, systems = [n, p1, p2])
end

"""
    ClosedFeedwaterHeater(; name)

DOCSTRING
  Steam_OFW CLOSED FEEDWATER HEATER
    @named p1 = BasicSteamPin()
    @named n1 = BasicSteamPin()
    @named p2 = BasicSteamPin()
    @named n2 = BasicSteamPin()

ELEM TYPE: COMPONENT
EQUATIONS:
    n1.P ~ p1.P # pressure
    n2.P ~ p2.P
    n2.h ~ n1.h # enthalpies
    0 ~ n1.ṁ + p1.ṁ
    0 ~ n2.ṁ + p2.ṁ
    yfrac ~ (n2.h - p2.h) / ((p1.h - n1.h) + (n2.h - p2.h))
    p1.ṁ ~ yfrac * (p1.ṁ + p2.ṁ)
    p1.Φ ~ yfrac * (p1.Φ + p2.Φ)
    0 ~ n1.Φ + p1.Φ + p1.ṁ * p1.h + n1.ṁ * n1.h
    0 ~ n2.Φ + p2.Φ + p2.ṁ * p2.h + n2.ṁ * n2.h
INPUTS
- name: Name of the system, symbol. Or use the @named macro when initializing.
"""
function ClosedFeedwaterHeater(; name)

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
    sts = @variables yfrac(t) = 0.5 [bounds = (0.0, 1.0)]
    # continuous_events = [[yfrac ~ 1] => [yfrac ~ 1]
    #                     [yfrac ~ 0] => [yfrac ~ 0.0]]
    eqs = [
        n1.P ~ p1.P # pressure
        n2.P ~ p2.P
        n2.h ~ n1.h # enthalpies
        0 ~ n1.ṁ + p1.ṁ
        0 ~ n2.ṁ + p2.ṁ
        yfrac ~ (n2.h - p2.h) / ((p1.h - n1.h) + (n2.h - p2.h))
        p1.ṁ ~ yfrac * (p1.ṁ + p2.ṁ)
        p1.Φ ~ yfrac * (p1.Φ + p2.Φ)
        0 ~ n1.Φ + p1.Φ + p1.ṁ * p1.h + n1.ṁ * n1.h
        0 ~ n2.Φ + p2.Φ + p2.ṁ * p2.h + n2.ṁ * n2.h
    ]
    # 0 ~ n1.ṁ*n1.h+n2.ṁ*n2.h+p1.ṁ*p1.h+p2.ṁ*p2.h
    # n1.Φ + p1.Φ ~  n2.Φ  + p2.Φ
    # n1.Φ + p1.Φ ~  n2.Φ  + p2.Φ
    # p2.ṁ * (n2.h - p2.h) ~ p1.ṁ * (n1.h-p1.h)
    # 0 ~ n1.ṁ*n1.h+n2.ṁ*n2.h+p1.ṁ*p1.h+p2.ṁ*p2.h
    ODESystem(
        eqs,
        t,
        sts,
        [];
        name = name,
        systems = [p1, n1, p2, n2],
        defaults = [yfrac => 0.5],
    )
end

"""
    throttle(; name)

DOCSTRING


ELEM TYPE: COMPONENT
EQUATIONS:
n.ṁ + p.ṁ ~ 0
p.h ~ n.h
ΔP ~ p.P - n.P
n.Φ ~ 0
INPUTS
- name: Name of the system, symbol. Or use the @named macro when initializing.
"""
function throttle(; name)
    # flows x and y are the inlets
    @named p = BasicSteamPin()
    @named n = BasicSteamPin()

    sts = @variables Ė_loss(t) = 0.5 ΔP(t) = 0.0
    eqs = [
        n.ṁ + p.ṁ ~ 0
        p.h ~ n.h
        ΔP ~ p.P - n.P
        n.Φ ~ 0
    ]         #conservation of energy
    ODESystem(
        eqs,
        t,
        sts,
        [];
        name = name,
        systems = [p, n],
        defaults = [Ė_loss => 0.0, ΔP => 0.0],
    )
end

