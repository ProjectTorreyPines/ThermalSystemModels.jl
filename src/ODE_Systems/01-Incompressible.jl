# module ThermoIncompressible
using ModelingToolkit, Plots, Symbolics, Logging
using Printf
using Plots
# NonlinearSolve,   DifferentialEquations, #CoolProp,
@independent_variables t
# Logging.disable_logging(Logging.Warn)
# include("03-MTK_UTILS.jl")
# data from [2]
cppb(x) = (0.195 - 9.116e-6 .* x) .* 1000  # J/kgK
vpb(x) = 1 / (10520.35 - 1.19051 .* x)      # m³/kg
spb(x) = 0.195 * 1e3 * log(x) - 9.11e-6 * 1000 * x - 1210.316002277804
hpb(x, y) = cppb(x) * x + vpb(x) * y * 1e5

Lcpfunc(x)      = cppb(x)
Lvfunc(x)       = vpb(x)
Lsfunc(x)       = spb(x)
Lhfunc(x, y)    = hpb(x, y)

@variables x
@register_symbolic Lcpfunc(x)
@register_symbolic Lvfunc(x)
@register_symbolic Lsfunc(x)
@register_symbolic Lhfunc(x, y)
propDict = Dict(Lcpfunc => cppb, Lvfunc => vpb, Lsfunc => spb, Lhfunc => hpb)


"""
    IncompressiblePin(; name, Pdef = 50, Tdef = 555, ṁdef = 0.0)

DOCSTRING
Incompresible Fluid pin,
    across_var = @variables P(t) = Pdef T(t) = Tdef s(t) = 1.0 cp(t) = 187 v(t) = 0.001
    thru_var = @variables ṁ(t) = ṁdef Φ(t) = 1.0                     # mass flow and energy flow
    sts = [T, P, ṁ, cp, s, v, Φ]
ELEM TYPE: PIN
EQUATIONS:
    cp ~ Lcpfunc(T)
    v ~ Lvfunc(T)
    s ~ Lsfunc(T)

INPUTS
- `name`: Name of the system, symbol. Or use the @named macro when initializing.
- `Pdef = 50`:  Initial Pressure (bar)
- `Tdef = 555`: Initial Temperature (K)
- `ṁdef = 0.0`: Initial Flow Rate (kg/s)
"""
function IncompressiblePin(; name, Pdef = 50, Tdef = 555, ṁdef = 0.0)
    across_var = @variables P(t) = Pdef T(t) = Tdef s(t) = 1.0 cp(t) = 187 v(t) = 0.001
    thru_var = @variables ṁ(t) = ṁdef Φ(t) = 1.0                     # mass flow and energy flow

    eq = [
        cp ~ Lcpfunc(T)
        v ~ Lvfunc(T)
        s ~ Lsfunc(T)
    ]
    sts = [T, P, ṁ, cp, s, v, Φ]
    ODESystem(
        eq,
        t,
        sts,
        [];
        name = name,
        defaults = [P => Pdef, T => Tdef, ṁ => ṁdef, Φ => 0.0],
    )
end

"""
    WorkPin(; name)

DOCSTRING
    Basic work pin
    sts = @variables Ẇ(t) = 0.0

ELEM TYPE: PIN
EQUATIONS:

INPUTS
- `name`: Name of the system, symbol. Or use the @named macro when initializing.
"""
function WorkPin(; name)
    sts = @variables Ẇ(t) = 0.0
    ODESystem(Equation[], t, sts, []; name = name)
end

"""
    HeatTransferPin(; name)

DOCSTRING
sts = @variables Q̇(t) = 0.0

ELEM TYPE: PIN
EQUATIONS:

INPUTS
- `name`: Name of the system, symbol. Or use the @named macro when initializing.
"""
function HeatTransferPin(; name)  #input in W
    sts = @variables Q̇(t) = 0.0
    ODESystem(Equation[], t, sts, []; name = name, defaults = [Q̇ => 0.0])
end

"""
    FixedHeatFlowPin(; name, Qin = 1000.0)

DOCSTRING
Component with fixed heat flow rate
sts = @variables Q̇(t) = Qin
ps = @parameters Qin = Qin

ELEM TYPE: PIN
EQUATIONS:
eqs = [Q̇ ~ Qin]
INPUTS
- `name`: Name of the system, symbol. Or use the @named macro when initializing.
- `Qin = 1000.0`: Heating rate, (W)
"""
function FixedHeatFlowPin(; name, Qin = 1e3)  #input in W
    sts = @variables Q̇(t) = Qin
    ps = @parameters Qin = Qin
    eqs = [Q̇ ~ Qin]
    ODESystem(eqs, t, sts, ps; name = name, defaults = [Q̇ => Qin])
end

"""
    IncompressibleOnePort(; name)

DOCSTRING
Oneport element for connecting single port components
    p = IncompressiblePin()
    n = IncompressiblePin()
    sts = @variables LHS(t) ΔP(t) = 0.0 ΔT(t) = 0.0

ELEM TYPE: PORT
EQUATIONS:
    LHS ~ p.Φ + n.Φ             # conservation of energy
    0 ~ p.ṁ + n.ṁ               # mass flow
    ΔP ~ n.P - p.P              # pressure drop
    ΔT ~ n.T - p.T              # Temerature change
INPUTS
- `name`: Name of the system, symbol. Or use the @named macro when initializing.
"""
function IncompressibleOnePort(; name)
    @named p = IncompressiblePin()
    @named n = IncompressiblePin()

    sts = @variables LHS(t) ΔP(t) = 0.0 ΔT(t) = 0.0

    eqs = [
        LHS ~ p.Φ + n.Φ             # conservation of energy
        0 ~ p.ṁ + n.ṁ               # mass flow
        ΔP ~ n.P - p.P
        ΔT ~ n.T - p.T
    ]
    ODESystem(eqs, t, sts, []; name = name, systems = [p, n])
end

"""
    SetPressure(; name, P = 0.1)

DOCSTRING
Fixed pressure element

ELEM TYPE: SOURCE
EQUATIONS:

INPUTS
- `name`: Name of the system, symbol. Or use the @named macro when initializing.
- `P = 0.1`: Pressure (bar)
"""
function SetPressure(; name, P = 0.1)
    @named oneport = IncompressibleOnePort()
    @unpack ΔT, ΔP, LHS, p = oneport
    ps = @parameters P = P
    eqs = [
        LHS ~ 0
        ΔT ~ 0
        ΔP ~ 0
        p.P ~ P
    ]
    extend(ODESystem(eqs, t, [], ps; name = name), oneport)
end

"""
    SetPressure2(; name, P = 0.1)

DOCSTRING
    Fixed pressure element that does not use the onport function
    @named p = IncompressiblePin(Pdef = P)
    @named n = IncompressiblePin(Pdef = P)
    ps = @parameters P = P

ELEM TYPE: SOURCE
EQUATIONS:
    p.P ~ P
    n.P ~ p.P
    n.T ~ p.T
    0 ~ p.Φ + n.Φ             # conservation of energy
    0 ~ p.ṁ + n.ṁ
INPUTS
- `name`: Name of the system, symbol. Or use the @named macro when initializing.
- `P = 0.1`: Pressure (bar)
"""
function SetPressure2(; name, P = 0.1)
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
    ODESystem(eqs, t, [], ps; name = name, systems = [p, n])
end

"""
    SetTemperature(; name, T = 300)

DOCSTRING
Fixed temperature element
    @named oneport = IncompressibleOnePort()
    @unpack ΔT, ΔP, LHS, p = oneport
    ps = @parameters T = T

ELEM TYPE: SOURCE
EQUATIONS:
    p.T ~ T
    ΔT ~ 0
    ΔP ~ 0
    LHS ~ 0
INPUTS
- `name`: Name of the system, symbol. Or use the @named macro when initializing.
- `T = 300`: Temperature (K)
"""
function SetTemperature(; name, T = 300)
    @named oneport = IncompressibleOnePort()
    @unpack ΔT, ΔP, LHS, p = oneport
    ps = @parameters T = T
    eqs = [
        p.T ~ T
        ΔT ~ 0
        ΔP ~ 0
        LHS ~ 0
    ]
    extend(ODESystem(eqs, t, [], ps; name = name), oneport)
end

"""
    IncompressibleFlowValve(; name)

DOCSTRING
A flow valve object with a mass flow variable. Can be controlled with an external controller element.
    @named p = IncompressiblePin()
    @named n = IncompressiblePin()
    sts = @variables ṁ(t) = 1.0

ELEM TYPE: COMPONENT
EQUATIONS:
    0 ~ p.Φ + n.Φ   # conservation of energy
    0 ~ p.ṁ + n.ṁ
    p.ṁ ~ ṁ
    p.T ~ n.T
    n.P ~ p.P
INPUTS
- `name`: Name of the system, symbol. Or use the @named macro when initializing.
"""
function IncompressibleFlowValve(; name)
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
    ODESystem(eqs, t, sts, []; name = name, systems = [p, n], defaults = [ṁ => 1.0])
end

"""
    IncompressibleFixedFlowSource(; name, Ṁ = 1.0)

DOCSTRING
    A fixed flow rate source
    @named p = IncompressiblePin()
    @named n = IncompressiblePin()
    ps = @parameters Ṁ = Ṁ

ELEM TYPE: COMPONENT
EQUATIONS:
    0 ~ p.Φ + n.Φ   # conservation of energy
    0 ~ p.ṁ + n.ṁ
    p.ṁ ~ Ṁ
    p.T ~ n.T
    n.P ~ p.P
INPUTS
- `name`: Name of the system, symbol. Or use the @named macro when initializing.
- `Ṁ = 1.0`: mass flow rate of the fixed flow source
"""
function IncompressibleFixedFlowSource(; name, Ṁ = 1.0)
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
    ODESystem(eqs, t, [], ps; name = name, systems = [p, n])
end

"""
    PassiveIncompressiblePump(; name, η = 0.9)

DOCSTRING


ELEM TYPE: COMPONENT
EQUATIONS:

INPUTS
- `name`: Name of the system, symbol. Or use the @named macro when initializing.
- `η = 0.9`: DESCRIPTION
"""
function PassiveIncompressiblePump(; name, η = 0.9)
    @named oneport = IncompressibleOnePort()
    @named w = WorkPin()
    @unpack ΔP, LHS, p, n = oneport
    ps = @parameters η = η

    eqs = [
        LHS ~ w.Ẇ
        w.Ẇ ~ p.ṁ * p.v * ΔP * 1e5 / η
        n.T ~ p.T + (w.Ẇ - p.v * ΔP) / p.cp # work, multiply by 100 to get to kPa
    ]
    extend(ODESystem(eqs, t, [], ps; systems = [w], name = name), oneport)
end

"""
    PassiveIncompressiblePump2(; name, η = 1.0)

DOCSTRING
A passive pump component for incompressible liquids. Calculates work required based off of pressures at the p and n port
    @named p = IncompressiblePin()
    @named n = IncompressiblePin()
    @named w = WorkPin()
    ps = @parameters η = η

ELEM TYPE: COMPONENT
EQUATIONS:
    0 ~ p.ṁ + n.ṁ #p.ṁ ~ n.ṁ                               # conservation of mass
    w.Ẇ ~ p.ṁ * p.v * (n.P - p.P) * 1e5 / η
    w.Ẇ ~ p.Φ + n.Φ
    n.T ~ p.T + (p.v * (n.P - p.P) * 1e5 / η - p.v * (n.P - p.P) * 1e5) / p.cp
INPUTS
- `name`: Name of the system, symbol. Or use the @named macro when initializing.
- `η = 1.0`: Isentropic effeciency of the pump
"""
function PassiveIncompressiblePump2(; name, η = 1.0)
    @named p = IncompressiblePin()
    @named n = IncompressiblePin()
    @named w = WorkPin()
    ps = @parameters η = η

    eqs = [
        0 ~ p.ṁ + n.ṁ #p.ṁ ~ n.ṁ                               # conservation of mass
        w.Ẇ ~ p.ṁ * p.v * (n.P - p.P) * 1e5 / η
        w.Ẇ ~ p.Φ + n.Φ
        n.T ~ p.T + (p.v * (n.P - p.P) * 1e5 / η - p.v * (n.P - p.P) * 1e5) / p.cp
    ]
    ODESystem(eqs, t, [], ps; name = name, systems = [p, n, w])
end

"""
    IncompressibleHeatTransfer(; name)

DOCSTRING
Basic heat transfer port. Includes p,n ports for flow, and q which is a heat transfer port
    @named p = IncompressiblePin()
    @named n = IncompressiblePin()
    @named q = HeatTransferPin()
    st = @variables Q̇(t) = 0.0 C(t) = 187   , Heat rate and duty of the flow

ELEM TYPE: COMPONENT
EQUATIONS:
    0 ~ p.ṁ + n.ṁ           # p.ṁ ~ n.ṁ                               # conservation of mass
    q.Q̇ ~ p.Φ + n.Φ         # conservation of energy
    C ~ p.ṁ * p.cp          # duty
    0 ~ q.Q̇ - Q̇
    n.T ~ p.T + q.Q̇ / C
    n.P ~ p.P
INPUTS
- `name`: Name of the system, symbol. Or use the @named macro when initializing.
"""
function IncompressibleHeatTransfer(; name)
    @named p = IncompressiblePin()
    @named n = IncompressiblePin()
    @named q = HeatTransferPin()

    st = @variables Q̇(t) = 0.0 C(t) = 187

    eqs = [
        0 ~ p.ṁ + n.ṁ           # p.ṁ ~ n.ṁ                               # conservation of mass
        q.Q̇ ~ p.Φ + n.Φ         # conservation of energy
        C ~ p.ṁ * p.cp          # duty
        0 ~ q.Q̇ - Q̇
        n.T ~ p.T + q.Q̇ / C
        n.P ~ p.P
    ]
    ODESystem(
        eqs,
        t,
        [Q̇, C],
        [];
        name = name,
        systems = [p, n, q],
        defaults = [Q̇ => 0.0, C => 187],
    )
end

"""
    FlowControlIncompressibleHeatTransfer(; name, ΔP = 0.0, Tout = 1000.0)

DOCSTRING
A heat transfer element with the ability to change the flow rate to achieve a desired outlet temperature, specified by the inputs.
    @named p = IncompressiblePin()
    @named n = IncompressiblePin()
    @named q = HeatTransferPin()
    st = @variables Q̇(t) = 0.0 C(t) = 187 ṁ(t) = 1.0
    ps = @parameters ΔP = ΔP Tout = Tout

ELEM TYPE: COMPONENT
EQUATIONS:
        n.T ~ Tout
        ṁ ~ q.Q̇ / ((n.T - p.T) * p.cp)           # q = mcp*ΔT
        p.ṁ ~ ṁ
        0 ~ p.ṁ + n.ṁ           # p.ṁ ~ n.ṁ                               # conservation of mass
        q.Q̇ ~ p.Φ + n.Φ         # conservation of energy
        C ~ ṁ * p.cp            # duty
        0 ~ q.Q̇ - Q̇
        n.P ~ p.P - ΔP
INPUTS
- `name`: Name of the system, symbol. Or use the @named macro when initializing.
- `ΔP = 0.0`: Pressure drop across the heat transfer element (bar)
- `Tout = 1000.0`: outlet temperature (K)
"""
function FlowControlIncompressibleHeatTransfer(; name, ΔP = 0.0, Tout = 1000.0)
    @named p = IncompressiblePin()
    @named n = IncompressiblePin()
    @named q = HeatTransferPin()

    st = @variables Q̇(t) = 0.0 C(t) = 187 ṁ(t) = 1.0
    ps = @parameters ΔP = ΔP Tout = Tout
    eqs = [
        n.T ~ Tout
        ṁ ~ q.Q̇ / ((n.T - p.T) * p.cp)           # q = mcp*ΔT
        p.ṁ ~ ṁ
        0 ~ p.ṁ + n.ṁ           # p.ṁ ~ n.ṁ                               # conservation of mass
        q.Q̇ ~ p.Φ + n.Φ         # conservation of energy
        C ~ ṁ * p.cp          # duty
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
        defaults = [Q̇ => 0.0, C => 187],
    )
end

"""
    SinglePortReservoir(; name, P = 0.1, T = 300)

DOCSTRING
A sinlge port flow energy ground to connect with circiuts, only has n pin
    @named n = IncompressiblePin(Pdef = P, Tdef = T)
    ps = @parameters P = P T = T
ELEM TYPE: COMPONENT
EQUATIONS:
    n.P ~ P
    n.T ~ T
    n.Φ ~ 0
INPUTS
- `name`: Name of the system, symbol. Or use the @named macro when initializing.
- `P = 0.1`: Pressure, (bar)
- `T = 300`: Temperature (K)
"""
function SinglePortReservoir(; name, P = 0.1, T = 300)
    @named n = IncompressiblePin(Pdef = P, Tdef = T)
    ps = @parameters P = P T = T
    eqs = [
        n.P ~ P
        n.T ~ T
        n.Φ ~ 0
    ]
    ODESystem(eqs, t, [], ps; name = name, systems = [n])
end

"""
    TwoPortReservoir(; name, P = 0.1, T = 300)

DOCSTRING
A flow conserving ground element with 2 ports.
    @named oneport = IncompressibleOnePort()
    @unpack ΔT, ΔP, LHS, n, p = oneport
    ps = @parameters T = T P = P

ELEM TYPE: COMPONENT
EQUATIONS:
    LHS ~ 0
    n.Φ ~ 0
    p.T ~ T
    p.P ~ P
    ΔP ~ 0
    ΔT ~ 0
INPUTS
- `name`: Name of the system, symbol. Or use the @named macro when initializing.
- `P = 0.1`: Pressure
- `T = 300`: Bar
"""
function TwoPortReservoir(; name, P = 0.1, T = 300)
    @named oneport = IncompressibleOnePort()
    @unpack ΔT, ΔP, LHS, n, p = oneport
    ps = @parameters T = T P = P
    eqs = [
        LHS ~ 0
        n.Φ ~ 0
        p.T ~ T
        p.P ~ P
        ΔP ~ 0
        ΔT ~ 0
    ]

    extend(
        ODESystem(eqs, t, [], ps; name = name, defaults = [T => 300.0, P => 10]),
        oneport,
    )
end

"""
    throttle(; name)

DOCSTRING
An enthalpy conserving throttle valve. Does not control the outlet pressure, but mimics a throttle valve when placed between 2 sources wtih a pressure differential.

ELEM TYPE: COMPONENT
EQUATIONS:

INPUTS
- `name`: Name of the system, symbol. Or use the @named macro when initializing.
"""
function throttle(; name)
    @named oneport = IncompressibleOnePort()
    @unpack LHS, p, ΔT, n = oneport
    sts = @variables Ė_loss(t) = 0.0
    eqs = [
        LHS ~ p.Φ
        ΔT ~ 0
        Ė_loss ~ p.Φ + n.Φ
    ]
    extend(ODESystem(eqs, t, sts, []; name = name, defaults = [Ė_loss => 0.0]), oneport)
end

"""
    IdealCooler(; name)

DOCSTRING


ELEM TYPE: COMPONENT
EQUATIONS:

INPUTS
- `name`: Name of the system, symbol. Or use the @named macro when initializing.
"""
function IdealCooler(; name)
    @named p = IncompressiblePin()
    @named n = IncompressiblePin()
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

"""
    ReliefElement(; name)

DOCSTRING


ELEM TYPE: COMPONENT
EQUATIONS:
    n.ṁ + p.ṁ ~ 0           # has to be negative
    q.Q̇ ~ p.Φ + n.Φ         # conservation of energy
INPUTS
- `name`: Name of the system, symbol. Or use the @named macro when initializing.
"""
function ReliefElement(; name)
    @named p = IncompressiblePin()
    @named n = IncompressiblePin()
    @named q = HeatTransferPin()
    # 0 variables, 3 equations for pin
    eqs = [
        n.ṁ + p.ṁ ~ 0           # has to be negative
        q.Q̇ ~ p.Φ + n.Φ         # conservation of energy
    ]
    ODESystem(eqs, t, [], []; name = name, systems = [p, n, q])
end

"""
    incompressible_connect(pins...)

DOCSTRING


ELEM TYPE: N/A
EQUATIONS:
eqs = [
    sum(pin -> pin.ṁ, pins) ~ 0.0, # Mass
    sum(pin -> pin.Φ, pins) ~ 0.0, # Energy
]

for pin in pins[2:end]
    push!(eqs, pins[1].P ~ pin.P)   # p_1 = p_2, p_1 = p_3,...
    push!(eqs, pins[1].T ~ pin.T)   # p_1 = p_2, p_1 = p_3,...
end

INPUTS
pins..., a set of pins from components in which to apply conservation laws too
"""
function incompressible_connect(pins...)

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
