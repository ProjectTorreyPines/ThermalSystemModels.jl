
#===============================================================#
#                       NODE STRUCTURE                          #
#   node type - calculated    - will change during analysis
#               controlled    - will not change during analysis
#
#===============================================================#
abstract type model_variables end
# abstract type model_component end
abstract type fluidObj end
# abstract type gas_component   <: model_component end
abstract type thru_variable   <: model_variables end
abstract type across_variable <: model_variables end
abstract type component_nodes end
abstract type energy_reservoir end


Base.@kwdef mutable struct gas_node <: component_nodes
    working_fluid::Symbol               = :helium   # fluid [:helium]
    mass_flow::Union{Int64,Float64}     = 1         #   mass flow kg/s
    T::Union{Float64,Int64}             = 300       #   temperature (K)
    P::Union{Float64,Int64}             = 1         #   pressure (Bar)

    Ptype::Symbol                       = :calculated
    Ttype::Symbol                       = :calculated
    fluid_type::Symbol                  = :gas
    node_num::Int                       =   1
end

Base.@kwdef mutable struct liq_node <: component_nodes
    working_fluid::Symbol               = :water        # fluid [:helium]
    mass_flow::Union{Int64,Float64}     = 1             #   mass flow kg/s
    T::Union{Float64,Int64}             = 300           #   temperature (K)
    P::Union{Float64,Int64}             = 1             #   pressure (Bar)
    x::Union{Float64,Int64}             = 1             #   vapour fraction
    h::Union{Float64,Int64}             = 0.0

    Ptype::Symbol                       = :calculated
    Ttype::Symbol                       = :calculated

    fluid_type::Symbol                  =   :liq
    node_num::Int                       =   1
end

mutable struct siso_nodes <: component_nodes
    A::component_nodes #outlet
    B::component_nodes #inlet
end

mutable struct simo_nodes <: component_nodes
    A::component_nodes #outlet
    B::Vector{component_nodes} #inlet
end

mutable struct miso_nodes <: component_nodes
    A::Vector{component_nodes} #outlet
    B::component_nodes  #inlet
end

siso_nodes(AB::Vector{<:component_nodes}) = siso_nodes(AB[1],AB[2])

function inlet_nodes(nodes::siso_nodes)
    return [nodes.A.node_num], [nodes.A.P], [nodes.A.T]
end

function inlet_nodes(nodes::simo_nodes)
    return [nodes.A.node_num], [nodes.A.P], [nodes.A.T]
end

function inlet_nodes(nodes::miso_nodes)
    nd = [a.node_num for a in nodes.A]
    press = [b.P for b in nodes.A]
    Temp = [b.T for b in nodes.A]
    return nd, press, Temp
end

function outlet_nodes(nodes::siso_nodes)
    return [nodes.B.node_num], [nodes.B.P], [nodes.B.T]
end

function outlet_nodes(nodes::simo_nodes)
    nd = [b.node_num for b in nodes.B]
    press = [b.P for b in nodes.B]
    Temp = [b.T for b in nodes.B]
    return nd, press, Temp
end

function outlet_nodes(nodes::miso_nodes)
    return [nodes.B.node_num],[nodes.B.P], [nodes.B.T]
end


# getters and setters
function get_pressure(n::component_nodes)
    return n.P
end

function get_temperature(n::component_nodes)
    return n.T
end

function get_multi_temp(ns::Vector{<:component_nodes})
    ll = length(ns)     #for speed
    TV = zeros(1,ll)    #temperature vecotr
    for i = 1:ll
        TV[i] = get_temperature(ns[i])
    end
    return TV
end

function get_mass_flow(n::component_nodes)
    return n.mass_flow
end

function get_fluid_type(n::component_nodes)
    return n.fluid_type
end

function get_working_fluid(n::component_nodes)
    return n.working_fluid
end

function get_node_fix_type(n::component_nodes)
    return n.Ptype, n.Ttype
end

function get_volume(n::component_nodes)
    return n.v
end

function set_working_fluid(n::component_nodes,fl::Symbol)
    n.working_fluid = fl
end

function fix_temperature(n::component_nodes)
    n.Ttype = :controlled
end

function fix_temperature(n::component_nodes,Tfix::Float64)
    n.T = Tfix
    n.Ttype = :controlled
end

function fix_pressure(n::component_nodes)
    n.Ptype = :controlled
end

function fix_pressure(n::component_nodes,Pfix::Float64)
    n.P = Pfix
    n.Ptype = :controlled
end

function overide_set_temperature(n::component_nodes,Tin::Float64)
    n.T = Tin
end

function overide_set_pressure(n::component_nodes,Pin::Float64)
    n.P = Pin
end

function set_pressure(n::component_nodes,Pin::Union{Float64,Int64}; verbose = false)
    if n.Ptype == :calculated
        n.P = Pin
    elseif verbose
        println("Attempt to set pressure of fixed node")
    end
    return n.P
end

function set_x(ln::liq_node,x::Real)
    ln.x = x
end

function set_h(ln::liq_node,h::Real)
    ln.h = h
end

function get_h(ln::liq_node)
    return ln.h
end


function set_temperature(n::component_nodes,Tin::Union{Float64,Int64})
    if n.Ttype == :calculated
        n.T = Tin
    else
        println("Attempt to set temperature of fixed node")
        println("    Node # $(n.node_num) , Fluid = $(n.working_fluid)")
    end
end

function set_mass_flow(n::component_nodes, mflow::Union{Float64,Int64}; use_existing_mflow_as_reference = false)
    #   if use_existing_mflow_as_reference, then mass flow - 0->1 @ all compoenents
    n.mass_flow = use_existing_mflow_as_reference ? n.mass_flow * mflow : mflow
end

function set_mass_flow(nvec::Vector{<:component_nodes},mflow::Union{Float64,Int64}; use_ref = false)
    for node in nvec
        set_mass_flow(node,mflow; use_existing_mflow_as_reference = use_ref)
    end
end

function set_mass_flow(nds::siso_nodes,mflow::Union{Float64, Int64}; use_ref = false)
    set_mass_flow(nds.A,mflow; use_existing_mflow_as_reference = use_ref)
    set_mass_flow(nds.B,mflow; use_existing_mflow_as_reference = use_ref)
end

function get_inlet_temp(nds::siso_nodes)
    return get_temperature(nds.A)
end

function get_inlet_press(nds::siso_nodes)
    return get_pressure(nds.A)
end

function get_inlet_volume(nds::siso_nodes)
    return get_volume(nds.A)
end

function set_inlet_temp(nds::siso_nodes,Tin::Union{Float64,Int64})
    set_temperature(nds.A,Tin)
end

function set_inlet_press(nds::siso_nodes,Pin::Union{Float64,Int64})
    set_pressure(nds.A,Pin)
end

function get_outlet_temp(nds::siso_nodes)
    return get_temperature(nds.B)
end

function get_outlet_press(nds::siso_nodes)
    return get_pressure(nds.B)
end

function get_outlet_press(nds::simo_nodes)
    Pout = [get_pressure(nd) for nd in nds.B]
    return Pout
end

function get_outlet_volume(nds::siso_nodes)
    return get_volume(nds.B)
end

function set_outlet_temp(nds::siso_nodes,Tin::Union{Float64,Int64})
    set_temperature(nds.B,Tin)
end

function set_outlet_press(nds::siso_nodes,Pin::Union{Float64,Int64})
    set_pressure(nds.B,Pin)
end

function get_TP(nds::siso_nodes)
    T1 = get_inlet_temp(nds)
    P1 = get_inlet_press(nds)
    T2 = get_outlet_temp(nds)
    P2 = get_outlet_press(nds)
    return T1,P1,T2,P2
end

function get_node_fix_type(sn::siso_nodes)
    inP,inT = get_node_fix_type(sn.A)
    oP,oT   = get_node_fix_type(sn.B)
    return inP,inT, oP,oT
end

# constructors
gas_node()                                  = gas_node(working_fluid = :helium, T = 300, P = 1,   mass_flow = 1) 
gas_node_P(pin::Float64; kw...)             = gas_node(P = pin; kw...) 
passive_gas_node(;kw...)                    = gas_node(Ptype = :calculated, Ttype = :calculated, kw...)
init_gas_nodes(numN::Int64)                 = [gas_node(node_num = i) for i=1:numN]
init_gas_nodes(numN::Int64,Pmin::Float64)   = [gas_node_P(Pmin; node_num = i) for i=1:numN]

liq_node()                                      = liq_node(working_fluid = :pbli, T = 600, P = 45,     mass_flow = 1) 
liq_node_P(pin::Float64; kw...)                 = liq_node(working_fluid = :pbli,  P = pin; kw...)  
iwnit_liq_nodes(numN::Int64)                    = [liq_node(node_num = i, working_fluid = :pbli) for i=1:numN]
init_liq_nodes(numN::Int64,Pmin::Float64)       = [liq_node_P(Pmin; node_num = i) for i=1:numN]


init_water_node(;kw...)                       = liq_node(working_fluid = :pbli, kw...)
init_water_nodes(numN::Int64)                 = [liq_node(working_fluid = :water ,node_num = i) for i=1:numN]


liq_properties(n::component_nodes) = liq_properties(n.working_fluid;P = n.P, T = n.T, x = n.x)
gas_properties(g::component_nodes) = gas_properties(g.working_fluid;P=g.P,T=g.T)


function gastype(working_fluid::Symbol)

    if working_fluid in [:air, :helium] 
        return :gas
    else
        return :liq
    end


end

function calculate_fluid_properties(n::component_nodes)
    fluid =  n.working_fluid
    fluid = enforce_lowercase(fluid)
    ftype = gastype(fluid)
    if ftype == :liq
        return liq_properties(n)
    else
        return gas_properties(n)
    end
end

function calculate_fluid_properties(nds::siso_nodes)
    fluid_props_a = calculate_fluid_properties(nds.A)
    fluid_props_b = calculate_fluid_properties(nds.B)
    return fluid_props_a, fluid_props_b
end

function get_inlet_duty(nds::siso_nodes; verbose = false)
    fluid_props = calculate_fluid_properties(nds.A)
    mflow = get_mass_flow(nds.A)
    if verbose
        @printf "\t massflow = %.2f, cp = %.2f\n" mflow fluid_props.cp
    end
    return mflow * fluid_props.cp
end
"""
    nodal_changes(nds::siso_nodes)

Calculates specific energy transfer, temperature, and pressure across 2 nodes
"""
function nodal_changes(nds::siso_nodes)
    # calculates the dh, du, ds of specific 
    fluid_props_a, fluid_props_b = calculate_fluid_properties(nds)
    cp_ave = (fluid_props_a.cp + fluid_props_b.cp) /2

    T1,P1,T2,P2 = get_TP(nds)
    ΔT = T2 - T1
    ΔP = P2 - P1
    ΔP = ΔP * 1e5
    ftype = get_fluid_type(nds.A)
    fluid = get_working_fluid(nds.A)

    if fluid == :pbli
        # EVALUATION INCOMPRESSIBLE FLUID
        inpr    = PbLi_properties(nds.A.T, nds.A.P)
        opr     = PbLi_properties(nds.B.T, nds.B.P)
        v_ave = (fluid_props_a.v + fluid_props_b.v)/2
        Δu = cp_ave * ΔT
        vΔP = v_ave * ΔP
        Δh = Δu + vΔP
        Δs = cp_ave * log(T2/T1)
        return ΔT,ΔP,Δu,Δh,Δs
    end
    
    if ftype == :gas
        # IDEAL GAS EVALUATION
        R = fluid_props_a.R
        cv_ave = (fluid_props_a.cv + fluid_props_b.cv)/2
        Δu = cv_ave * ΔT
        Δh = cp_ave * ΔT
        Δs = cp_ave * log(T2/T1) - R * log(P2/P1)
        return ΔT,ΔP,Δu,Δh,Δs
    end

    if fluid == :water
        Δu = fluid_props_b.u - fluid_props_a.u
        Δh = nds.B.h - nds.A.h
        Δs = fluid_props_b.s - fluid_props_a.s
        return ΔT,ΔP,Δu,Δh,Δs
    end
end

function get_vP(nds::siso_nodes)
    inprops, outprops = calculate_fluid_properties(nds::siso_nodes)
    v1 = inprops.v
    P1 = get_pressure(nds.A)
    v2 = outprops.v
    P2 = get_pressure(nds.B)
    return v1,P1,v2,P2
end

function dh_nodes(nds::siso_nodes)
    in,out = calculate_fluid_properties(nds)
    return out.h - in.h
end