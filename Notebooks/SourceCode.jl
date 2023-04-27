
using Roots, DataFrames, XSteam, Printf


abstract type model_variables end
abstract type model_component end
abstract type fluidObj end
abstract type gas_component   <: model_component end
abstract type thru_variable   <: model_variables end
abstract type across_variable <: model_variables end
abstract type component_nodes end
abstract type component_parameters end
abstract type energy_reservoir end

#utilities
"""
    convert_pressure(val::Real,from_unit::Symbol,to_unit::Symbol)

    convert pressure units, optional to specify to_unit, default = bar because that is what XSTeam uses
"""
function convert_pressure(val::Real,from_unit::Symbol,to_unit::Symbol)
    # default = SI
    # converting to lowercase to avoid errors
    from_unit = Symbol(lowercase(string(from_unit)))
    to_unit = Symbol(lowercase(string(to_unit)))

    # unit vector (pun)
    units = (pa = 1.0, kpa = 1.0e3, mpa = 1.0e6, gpa = 1.0e9, psi = 6.8947572932e3, ksi = 6.8947572932e3 * 1000 ,bar = 100000.0, atm = 101325.0, inh2o = 249.082);

    inUnit = units[from_unit]
    outUnit = units[to_unit]
    return val/(outUnit/inUnit)
end

function convert_pressure(val::Real,from_unit::Symbol)
    # default = SI
    # converting to lowercase to avoid errors
    from_unit = Symbol(lowercase(string(from_unit)))
    to_unit = :bar

    # unit vector (pun)
    units = (pa = 1.0, kpa = 1.0e3, mpa = 1.0e6, gpa = 1.0e9, psi = 6.8947572932e3, ksi = 6.8947572932e3 * 1000 ,bar = 100000.0, atm = 101325.0, inh2o = 249.082);

    inUnit = units[from_unit]
    outUnit = units[to_unit]
    return val/(outUnit/inUnit)
end

function enforce_lowercase(x = insym) 
    return x = Symbol(lowercase(String(x)))
end

function hxeff(Thi::Real, Ch::Real, Tci::Real, Cc::Real, eff::Real)
    # if Ch and Cc are just cp, Qact = Qspecific
    Cmin = min(Ch, Cc)
    dt_mx = Thi - Tci
    Q_act = Cmin .* dt_mx .* eff
    dtH = Q_act ./ Ch
    dtC = Q_act ./ Cc
    Tho = Thi - dtH
    Tco = Tci + dtC
    return Q_act, Tho, Tco
end

function hx_pinch(Thi::Real, Ch::Real, Tci::Real, Cc::Real, approach::Real)
    if Ch >= Cc
        #dT2 = minimum approach
        Tco = Thi - approach
        Q   = Cc*(Tco - Tci)
        dtH =   Q/Ch
        Tho =   Thi - dtH
        return Q, Tho, Tco
    elseif Cc > Ch
        Tho = Tci + approach
        Q   = Ch*(Thi - Tho)
        dtC = Q/Cc
        Tco = Tci + dtC
        return Q, Tho, Tco
    end
end

function lmtd(Thi,Tho,Tci,Tco)
    # log mean temperature difference, cross flow only
    ΔT_1 = Thi - Tco
    ΔT_2 = Tho - Tci
    ΔT_logmean = (ΔT_1-ΔT_2)/log(ΔT_1/ΔT_2);
    return ΔT_1,ΔT_2,ΔT_logmean
end

repmat(v,n)  =  repeat(v;inner=n)
repelem(v,n) =  repeat(v;outer=n)


#===============================================================#
#                       NODE STRUCTURE                          #
#   node type - calculated    - will change during analysis
#               controlled    - will not change during analysis
#
#===============================================================#
Base.@kwdef mutable struct gas_node <: component_nodes
    working_fluid::Symbol               = :helium             # fluid [:helium]
    mass_flow::Union{Int64,Float64}     = 1     #   mass flow kg/s
    T::Union{Float64,Int64}             = 300   #   temperature (K)
    P::Union{Float64,Int64}             = 1     #   pressure (Bar)

    Ptype::Symbol                       = :calculated
    Ttype::Symbol                       = :calculated
    fluid_type::Symbol                  = :gas
    node_num::Int                       =   1
end

Base.@kwdef mutable struct liq_node <: component_nodes
    working_fluid::Symbol               = :water            # fluid [:helium]
    mass_flow::Union{Int64,Float64}     = 1     #   mass flow kg/s
    T::Union{Float64,Int64}             = 300   #   temperature (K)
    P::Union{Float64,Int64}             = 1     #   pressure (Bar)
    x::Union{Float64,Int64}             = 1     #   vapour fraction

    Ptype::Symbol                       = :calculated
    Ttype::Symbol                       = :calculated
    fluid_type::Symbol                  = :liq
    node_num::Int                       =   1
end

mutable struct siso_nodes <: component_nodes
    A::component_nodes #outlet
    B::component_nodes #inlet
end
siso_nodes(AB::Vector{<:component_nodes}) = siso_nodes(AB[1],AB[2])

# getters and setters
function get_pressure(n::component_nodes)
    return n.P
end
function get_temperature(n::component_nodes)
    return n.T
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

function set_pressure(n::component_nodes,Pin::Union{Float64,Int64})
    if n.Ptype == :calculated
        n.P = Pin
    else
        println("Attempt to set pressure of fixed node")
    end
    return n.P
end
function set_temperature(n::component_nodes,Tin::Union{Float64,Int64})
    if n.Ttype == :calculated
        n.T = Tin
    else
        println("Attempt to set temperature of fixed node")
    end
end
function set_mass_flow(n::component_nodes, mflow::Union{Float64,Int64}; use_existing_mflow_as_reference = false)
    #   if use_existing_mflow_as_reference, then mass flow - 0->1 @ all compoenents
    n.mass_flow = use_existing_mflow_as_reference ? c.mass_flow * mflow : mflow
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
gas_node_P(pin::Float64)                    = gas_node(working_fluid = :helium, T = 300, P = pin, mass_flow = 1) 
passive_gas_node(;kw...)                    = gas_node(Ptype = :calculated, Ttype = :calculated, kw...)
init_gas_nodes(numN::Int64)                 = [gas_node(node_num = i) for i=1:numN]
init_gas_nodes(numN::Int64,Pmin::Float64)   = [gas_node_P(Pmin) for i=1:numN]



# PROPERTY UTILITIES
#===============================================================#
#                       PROPERTY UTILITIES                      #
#===============================================================#

mutable struct gas_prop_obj 
    cp::Union{Float64,Int64}    # specific heat
    cv::Union{Float64,Int64}    # specific heat constant volume
    gam::Union{Float64,Int64}   # ratio of specific heats
    v::Union{Float64,Int64}     # specific volume
    R::Union{Float64,Int64}     # gas constant
end

mutable struct liq_prop_obj
    P::Real
    T::Real
    cp::Real    # specific heat
    v::Real     #   specific volume
    h::Real     # enthalpy (kJ/kg)
    u::Real     # specific energy
    s::Real     # entropy (kJ/kg)
    x::Real     # vapor fraction
end
pbli_fluid_obj(Tin::Union{Float64,Int64},Pin::Union{Float64,Int64},cpin::Union{Float64,Int64},vin::Union{Float64,Int64}) = liq_prop_obj(P=Pin,T=Tin,cp=cpin,v=vin,h = -1,s = -1, x = 0.0)

function gas_properties(gas::Symbol;P=1,T=300)
    #default pressure and temperature (1 bar, 300K)
    gas     = Symbol(lowercase(String(gas)))
    R_dict  = Dict(:helium => 2076.9,   :air => 287.0)  #"J/kg/K"
    cp_dict = Dict(:helium => 5192.6,   :air => 1005.0) #"J/kg/K"
    
    try
         cp = cp_dict[gas]    # constant pressure specific heat
         R  = R_dict[gas]     # Gas Constant
         cv = cp-R                          # constant volume specific heat
         v = R*T/(P*1e5)
         k = cp/cv
         return gas_prop_obj(cp, cv,k, v,R)
    catch e
        println("uknown fluid type in dependent_gas_props()")
        display(e)
    end
end
function PbLi_properties(T::Union{Float64,Int64},P::Union{Float64,Int64})
    if T>1615.15
        println("Exceeded lead lithium boiling point $(T)")
    elseif T<453
        println("Warning - lead lithium temperature below melting point $(T)")
    end
    # data from [2]
    cp   = (0.195 - 9.116e-6 .* T) .* 1000  # J/kgK
    v    = 1/(10520.35 - 1.19051 .* T)      # m3/kg
    # λ    = (.1451 + 3.46* T)*100            # thermal conductivity
    return pbli_fluid_obj(T,P,cp,v)
end
# superheated water (above dome)
function initSuperHeatedWater(P::Real,T::Real)
    #  assumes vapor fraction = 1, i.e. fluid is 100% superheated 
    cp = Cp_pT(P,T-273.15);
    h = h_pT(P,T-273.15);
    v = v_ph(P,h)
    s = s_ph(P,h)
    u = u_ph(P,h)
    return water_prop_obj(P,T,v, cp*1e3,h*1e3,u*1e3,s*1e3,1);
end
# any water
function initWater_PT(P::Real,T::Real)
    #Tin in Kelvin, Pin in Bar
    Tref = Tin - 273.15;
    cp = Cp_pT(p,Tref)
    h = h_pT(p,Tref)
    v = v_pT(p,Tref)
    s = s_ph(p,h)
    u = u_pT(P,Tref)
    x = x_ph(p,h)
    return liq_prop_obj(P,T,v,cp*1e3,h*1e3,u*1e3,s*1e3,x)
end

"""
    initLiquidWater(P::Real)

    Initializes saturated water properties
"""
function initLiquidWater(P::Real)
    cp = CpV_p(P);
    h = h_px(P,0);
    v = vL_p(P);
    s = s_ph(P,h);
    x = x_ph(P,h);
    T = Tsat_p(P);
    u = u_pT(P,T)
    T = T+273.15
    return liq_prop_obj(P,T,v,cp*1e3,h*1e3,u*1e3,s*1e3,x);
end

function liq_properties(liq::Symbol;P = 1, T = 500, x=-1.0)
    liq = enforce_lowercase(liq)
    if liq == :pbli
        return PbLi_properties(T,P)
    elseif liq == :water
        if x != 1.0 && x != 0.0
            return initWater_PT(P,T)
        elseif x == 1
            return initSuperHeatedWater(P,T)
        elseif x == 0
            return initLiquidWater(P)
        end
        
    end

end

liq_properties(n::component_nodes) = liq_properties(n.working_fluid;P = n.P, T = n.T)
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

function get_inlet_duty(nds::siso_nodes)
    fluid_props = calculate_fluid_properties(nds.A)
    mflow = get_mass_flow(nds.A)
    return mflow * fluid_props.cp
end

"""
    nodal_changes(nds::siso_nodes)

Calculates specific energy transfer, temperature, and pressure across 2 nodes
"""
function nodal_changes(nds::siso_nodes)
    # calculates the dh, du, ds of specific 
    fluid_props_a, fluid_props_b = calculate_fluid_properties(nds)
    cp_ave = (fluid_props_a.cp + fluid_props_b.cp)/2

    T1,P1,T2,P2 = get_TP(nds)
    ΔT = T2 - T1
    ΔP = P2 - P1

    ftype = get_fluid_type(nds.A)
    fluid = get_working_fluid(nds.A)

    if fluid == :pbli
        # EVALUATION INCOMPRESSIBLE FLUID
        v_ave = (fluid_props_a.v + fluid_props_b.v)/2
        Δu = cp_ave * ΔT
        Pv = v_ave * ΔP
        Δh = Δu + Pv
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
        Δh = fluid_props_b.h - fluid_props_a.h
        Δs = fluid_props_b.s - fluid_props_a.s
        return ΔT,ΔP,Δu,Δh,Δs
    end
end

#===============================================================#
#                       COMPONENTS                              #
# parameter struct
# component struct
#   Nodes - sometimes multiple sets
#   All components have -   params::ParametersActor
#                           these are constants/qualities of the component, effeciency, compression ratio etc
#   Across variables ΔP, ΔT, ΔH, ΔU, ΔS
#   Through variables Q, W, mass_flow
#   
#   advance(c::component_type)  Computes the next state of the component (outlet) - 
#           acts ONLY on the node struct
#           sets mass flow of the nodes
#
#   process(c::component_type)  Computes the component variables        - acts ONLY on the comoponent struct
#   check(c::compoennt_type)    Ensure energy balance
#   
#   print funciton
#
#===============================================================#

##  HIGH LEVEL UTILITIES


### MASS FLOW SOURCE

Base.@kwdef mutable struct mass_flow_source_paremeter <: component_parameters
    working_fluid::Symbol       = :helium
    mass_flow::Float64  = 1.0
end
mass_flow_source_paremeter(massFlow::Union{Float64,Int64}) = mass_flow_source_paremeter(mass_flow = massFlow)
mass_flow_source_paremeter(massFlow::Union{Float64,Int64}, fluid::Symbol) = mass_flow_source_paremeter(mass_flow = massFlow, working_fluid = fluid)
function print_params(para::mass_flow_source_paremeter)
    @printf "\t Fluid = %s, Mass flow [kg/s] = %.2f \n" String(para.fluid)  para.mass_flow
end

Base.@kwdef mutable struct mass_flow_source <: model_component
    nodes::siso_nodes                   #  inlet and outlet
    params::mass_flow_source_paremeter = mass_flow_source_paremeter()

    #across variables
    ΔT::Union{Float64,Int64} = 0       # temperature
    ΔP::Union{Float64,Int64} = 0       # pressure change

    ΔH::Union{Float64,Int64} = 0       # enthalpy
    ΔU::Union{Float64,Int64} = 0       # internal energy
    ΔS::Union{Float64,Int64} = 0       # entropy

    Q::Union{Float64,Int64} = 0 
    W::Union{Float64,Int64} = 0
    #through variables
    mass_flow::Union{Float64,Int64}  = 1
    name::String = "mass_flow_source"
end

function advance(c::mass_flow_source)
    return
end

function process(c::mass_flow_source)
    return
end
mass_flow_source(AB::Vector{gas_node}; kw...)                           =   mass_flow_source(nodes = siso_nodes(AB); kw...)
mass_flow_source(AB::Vector{gas_node},massFlow::Union{Float64,Int64})   =   mass_flow_source(AB; params = mass_flow_source_paremeter(massFlow), mass_flow = massFlow)


#   set_mass_flow = sets the mass flow of a fluid 
function set_mass_flow(c:: model_component, mflow::Union{Float64,Int64}; use_existing_mflow_as_reference = false) 
    #   if use_existing_mflow_as_reference, then mass flow - 0->1 @ all compoenents
    c.mass_flow = use_existing_mflow_as_reference ? c.mass_flow * mflow : mflow
end

function set_mass_flow(cvec::Vector{model_component},mflow::Union{Float64,Int64};use_ref = false)
    for c in cvec
        if typeof(c) != mass_flow_source
            set_mass_flow(c,    mflow;   use_existing_mflow_as_reference = use_ref)
        end
    end
end
# ===================  TURBINE AND COMPRESSORS==============================
#     
#       parameters  rp, η
#
#       Constructors
#
# =============================================================================

Base.@kwdef mutable struct gas_working_element_parameters <: component_parameters
    rp::Float64 = 3.5   # compression ratio
    η::Float64  = 0.9   # isentropic effeciency
end

function print_params(para::gas_working_element_parameters)
    @printf " \t rp = %.2f, \t η =  %.2f\n " para.rp para.η
end

turb_params(;rpin::Float64) = gas_working_element_parameters(rp=rpin)
comp_params(;rpin::Float64) = gas_working_element_parameters(rp=rpin)

Base.@kwdef mutable struct gas_compressor <: model_component
    nodes::siso_nodes   # inlet and outlet
    params::gas_working_element_parameters = gas_working_element_parameters()

    #across variables
    ΔH::Union{Float64,Int64} = 0       #enthalpy (specific)
    ΔT::Union{Float64,Int64} = 0       #temperature
    ΔP::Union{Float64,Int64} = 0       #pressure change

    Q::Union{Float64,Int64} = 0        # Q
    W::Union{Float64,Int64} = 0        # W

    #through variables
    mass_flow::Union{Float64,Int64}  = 1

    name::String = "compressor"
end

function advance(c::gas_compressor)
    set_mass_flow(c.nodes,c.mass_flow)
    gprops = gas_properties(c.nodes.A)
    k = gprops.gam
    kcoeff = (k-1)/k
    ηc = c.params.η
    rp = c.params.rp
    a1c = (-1 + ηc + rp^kcoeff)/ηc;
    set_temperature(c.nodes.B,  a1c * c.nodes.A.T )  #outlet temperature
    set_pressure(c.nodes.B, c.nodes.A.P * rp)
end

Base.@kwdef mutable struct gas_turbine <: model_component
    nodes::siso_nodes   # inlet and outlet
    params::gas_working_element_parameters = gas_working_element_parameters()

    #across variables
    ΔH::Union{Float64,Int64} = 0       #enthalpy (specific)
    ΔT::Union{Float64,Int64} = 0       #temperature
    ΔP::Union{Float64,Int64} = 0       #pressure change

    Q::Union{Float64,Int64} = 0        # Q
    W::Union{Float64,Int64} = 0        # W
    
    #through variables
    mass_flow::Union{Float64,Int64}  = 1
    name::String = "Turbine"
end

function advance(t::gas_turbine)
    set_mass_flow(t.nodes,t.mass_flow)
    gprops = gas_properties(t.nodes.A)
    k = gprops.gam
    kcoeff = (k-1)/k
    ηt = t.params.η
    rp = t.params.rp
    a1t = (ηt + rp^kcoeff - ηt*(rp^kcoeff)) / (rp^kcoeff)

    set_temperature(t.nodes.B , a1t * t.nodes.A.T)       # outlet temperature
    set_pressure(t.nodes.B, t.nodes.A.P / rp)
end

function process(c::Union{gas_compressor,gas_turbine})
    gprops = gas_properties(c.nodes.A)
    cp  = gprops.cp
    c.ΔT =      c.nodes.B.T - c.nodes.A.T
    c.ΔH =      cp*(c.nodes.B.T-c.nodes.A.T)
    c.ΔP =      (c.nodes.B.P - c.nodes.A.P)
    c.W =       c.mass_flow * cp * (c.nodes.B.T-c.nodes.A.T)
end
gas_turbine(AB::Vector{gas_node}; kw...)                    =   gas_turbine(nodes = siso_nodes(AB); kw...)
gas_compressor(AB::Vector{gas_node};kw...)                  =   gas_compressor(nodes=siso_nodes(AB);kw...)

# =================== HEATERS==============================
#  All heater and cooler elements
#   Parameters
#       Heat_Modes -    Tset , (Default) heat load calculated from dT  
#                       Qset , Output temperature calculated from prescribed heat load
#                               Use Qset for regenerators
#                     
#       pressure_mode   :passive = (Default) pressure change calculated from nodes
#       pressure_mode   :active  = pressure change enforced by default_ΔP 
#       
#
#   Constructors       
#           intercooler(AB,Tout; ΔP) -> intercooler, forced pressure drop, outlet node = fixed temp
#           heat_sink(AB,Tout,Pout)  -> Heat rejection, passive component, outlet node = fixed temp, fixed pressure
#           ideal_heat_source(AB,Qin;ΔP) -> Prescribed heat load source
#           regen_hot(AB)
#           regen_cold(AB)
#           heat_exchanger_heater(AB,Qin;kw)
#                           kw...   ΔP = 0.2, 
#                                   fixed_T::Bool = false,              -> prescribed temperature difference
#                                   Ti::Float64 = -1, 
#                                   To::Float64 = -1 ->                 if fixed_T is true, either change the node temps prior, or set To, Ti in the constructor to define desired temperatures
# =============================================================================


Base.@kwdef mutable struct heating_element_parameters <: component_parameters
    heat_mode::Symbol = :Tset
    pressure_mode::Symbol = :passive
    default_ΔP::Float64 = 0.2
end
regen_heat_params() = heating_element_parameters(heat_mode = :Qset, pressure_mode = :active)

function print_params(para::heating_element_parameters)
    @printf "\t Heat Mode %s, \t pressure mode: %s, default ΔP: %.2f \n" String(para.heat_mode) String(para.pressure_mode) para.default_ΔP
end

Base.@kwdef mutable struct heat_source <: model_component
    nodes::siso_nodes                   #  inlet and outlet
    params::heating_element_parameters = heating_element_parameters()

    #across variables
    ΔT::Union{Float64,Int64} = 0       # temperature
    ΔP::Union{Float64,Int64} = 0       # pressure change

    ΔH::Union{Float64,Int64} = 0       # enthalpy
    ΔU::Union{Float64,Int64} = 0       # internal energy
    ΔS::Union{Float64,Int64} = 0       # entropy

    Q::Union{Float64,Int64} = 0 
    W::Union{Float64,Int64} = 0
    
    #through variables
    mass_flow::Union{Float64,Int64}  = 1
    name::String = "heat_source"
end
heat_source(nds::siso_nodes;kw...)                      = heat_source(nodes= nds; kw...)
heat_source(AB::Vector{<:component_nodes}; kw...)       = heat_source(nodes = siso_nodes(AB); kw...)
heat_source(AB::Vector{gas_node};kw...)                 = heat_source(nodes = siso_nodes(AB);kw...)

init_regen_hot(AB::Vector{<:component_nodes})           = heat_source(AB;params = regen_heat_params(), name = "regen_hot")
init_regen_cold(AB::Vector{<:component_nodes})          = heat_source(AB;params = regen_heat_params(), name = "regen_cold")

# default constructors
function intercooler(AB::Vector{<:component_nodes},Tout::Float64; ΔP = 0.2, name_string = "intercool")
    para = heating_element_parameters(heat_mode = :Tset, pressure_mode = :active, default_ΔP = ΔP)
    nds = siso_nodes(AB)
    fix_temperature(nds.B,Tout)
    return heat_source(nodes = nds; params = para, name = name_string)
end

function reheater(AB::Vector{<:component_nodes},Tout::Float64; ΔP = 0.2, name_string = "reheat")
    para = heating_element_parameters(heat_mode = :Qset, pressure_mode = :active, default_ΔP = ΔP)
    nds = siso_nodes(AB)
    return heat_source(nodes = nds; params = para, name = name_string)
end

function simple_reheater(AB::Vector{<:component_nodes},Tout::Float64; ΔP = 0.2, name_string = "reheat")
    para = heating_element_parameters(heat_mode = :Tset, pressure_mode = :active, default_ΔP = ΔP)
    nds = siso_nodes(AB)
    return heat_source(nodes = nds; params = para, name = name_string)
end

#   heat sink - very simple component, passive pressrue drop, set outlet conditions
function heat_sink(AB::Vector{<:component_nodes},Tout::Float64,Pout::Float64; name_string = "heat_sink")
    para = heating_element_parameters(heat_mode = :Tset, pressure_mode = :passive, default_ΔP = 0.0)
    nds = siso_nodes(AB)
    fix_temperature(nds.B,Tout)
    fix_pressure(nds.B,Pout)
    return heat_source(nodes = nds; params = para, name = name_string)
end

function ideal_heat_source(AB::Vector{<:component_nodes},Qin::Float64; ΔP = 0.2)
    para = heating_element_parameters(heat_mode = :Qset, pressure_mode = :active, default_ΔP = ΔP)
    return heat_source(AB; params = para, Q = Qin)
end

ideal_heat_source(AB::Vector{<:component_nodes},Qin::Float64) = heat_source(nodes = siso_nodes(AB); params = heating_element_parameters(heat_mode = :Qset, pressure_mode = :active, default_ΔP = 0.2), Q = Qin)

function heat_exchanger_heater(AB::Vector{<:component_nodes}; Qin = 0.0, ΔP = 0.2, fixed_T::Bool = false, Ti::Float64 = -1.0, To::Float64 = -1.0, name_string = "hx_noname")
    
    # fixed T - means that this side of the heat exchanger has fixed temperature differential
    nds = siso_nodes(AB)
    if fixed_T == false
        # not fixed T, meaning Q is calculated during hx evaluation - which controls temperature distribution , this is the more robust version
        para = heating_element_parameters(heat_mode = :Qset, pressure_mode = :active, default_ΔP = ΔP)
        return heat_source(nodes = nds; params = para, Q = Qin,name = name_string)

    else
        if Ti != -1
            fix_temperature(nds.A,Ti)
        end
        if To != -1
            fix_temperature(nds.B,To)
        end

        para = heating_element_parameters(heat_mode = :Tset, pressure_mode = :active, default_ΔP = ΔP)
        return heat_source(nodes = nds; params = para,name = name_string)
    end
end


# utitlities/checks
check_q_type(hs::heat_source)   = hs.params.heat_mode
check_p_type(hs::heat_source)   = hs.params.pressure_mode
get_Q(hs::heat_source)          = hs.Q
get_inlet_temp(hs::heat_source) = get_inlet_temp(hs.nodes)
get_inlet_duty(hs::heat_source) = get_inlet_duty(hs.nodes)
get_TP(hs::heat_source)         = get_TP(hs.nodes)

function set_Q(hs::heat_source,Q; mode::Symbol = :in, overide::Bool = false)
    if mode == :in
        Q = -abs(Q) # heat additino
    else
        Q = abs(Q)  # heat loss
    end
    q_type = check_q_type(hs)
    if q_type != :Qset && overide == false
        println("Attempt to change heat load manually, change kw args to overide = true")
        return
    else
        hs.Q = Q
    end
end

function advance(hs::heat_source)
    #calculate the next state of the heat source

    #set node mass flow
    set_mass_flow(hs.nodes,hs.mass_flow)

    if check_q_type(hs) == :Qset
        # set heat load, calcualating temperature rise
        inprops,outprops = calculate_fluid_properties(hs.nodes)
        Ti = get_inlet_temp(hs.nodes)
        Q = hs.Q                                # (-) means heat ADDITION mdot cp (T2-T1) = Q,  T2 = Q/mcp+T1
        mdot_cp = get_mass_flow(hs.nodes.A)
        Tout = Q/mdot_cp + Ti
        set_outlet_temp(hs.nodes,Tout)
    end

    if check_p_type(hs) != :passive
        p1 = get_inlet_press(hs.nodes)
        set_outlet_press(hs.nodes,p1 - hs.params.default_ΔP )
    end

end

function process(hs::heat_source)
    ptype = check_p_type(hs)
    qtype = check_q_type(hs)
    inprops,outprops    = calculate_fluid_properties(hs.nodes)
    T1,P1,T2,P2         = get_TP(hs.nodes)
    ΔT,ΔP,Δu,Δh,Δs      = nodal_changes(hs.nodes)
    hs.ΔT =      ΔT
    hs.ΔP =      ΔP
    hs.ΔU =      Δu * hs.mass_flow
    hs.ΔH =      Δh * hs.mass_flow
    hs.ΔS =      Δs * hs.mass_flow
    if qtype == :Tset
        hs.Q =   hs.ΔH
    end
end

#===============================================================#
#                       heat exchangers   - meta component      #
#       Modeled as 1x cooler, 1x heater
#       gas heater and cooler elements
#           hx_mode = :approach, hot fixed, cold fixed, effectiveness
#           Heat_Modes - Tset , heat load calculated from dT  
#           Qset , Output temperature calculated from prescribed heat load
#===============================================================#

Base.@kwdef mutable struct heat_exchanger_parameters <: component_parameters
    hx_mode = :eff        #Other options, :eff = effectiveness NTU method, :pinch
    ϵ::Float64 = .95
    min_approach::Float64 = 10.0
    UA = 0.0
end
get_hx_mode(hpar::heat_exchanger_parameters) = hpar.hx_mode

Base.@kwdef mutable struct heat_exchanger <: model_component
    #   Cases
    #   Thi,Tci, effectiveness
    hot_stream::heat_source
    cold_stream::heat_source

    params::heat_exchanger_parameters = heat_exchanger_parameters()

    #   across variables
    Q_hx::Union{Float64} = 0.0

    Cc::Union{Float64,Int64} = 0.0
    Tci::Real = 0.0
    Tco::Real = 0.0
    ΔT_c::Real = 0

    Ch::Union{Float64,Int64} = 0.0
    Thi::Real = 0.0
    Tho::Real = 0.0
    ΔT_h::Real = 0.0
    ΔT_lm::Real = 0.0

    UA::Real = 0.0
    name::String = "hx"
end
get_hx_mode(hpar::heat_exchanger) = get_hx_mode(hpar.params)

#       joins 2 heater elements
function heat_exchanger(cold_heater::heat_source,hot_heater::heat_source)
    return heat_exchanger(cold_stream = cold_heater,hot_stream = hot_heater)
end

function advance(hx::heat_exchanger)
    cold_type   = check_q_type(hx.cold_stream)
    hot_type    = check_q_type(hx.hot_stream)

    if cold_type == hot_type == :Tset
        # this cant be the case unless
        println("Heat exchanger joined to incompattible heating elements, both with prescribed temperature differentials")
        error("Advance heat exchagner error")
    end

    if cold_type == :Tset
        advance(hx.cold_stream)
        process(hx.cold_stream)
        @assert(hot_type == :Qset)
        set_Q(hx.hot_stream,hx.cold_stream.Q; mode = :out)

    elseif hot_type == :Tset
        advance(hx.hot_stream)
        process(hx.hot_stream)
        @assert(cold_type == :Qset)
        set_Q(hx.cold_stream,hx.hot_stream.Q; mode = :in)
        advance(hx.cold_stream)
        process(hx.cold_stream)
    else
        # both are Qset
        Cc = get_inlet_duty(hx.cold_stream)
        Ch = get_inlet_duty(hx.hot_stream)

        Tci = get_inlet_temp(hx.cold_stream)
        Thi = get_inlet_temp(hx.hot_stream)
        mode = get_hx_mode(hx)
        if mode == :approach
            Qx,Tho,Tco = hx_pinch(Thi,Ch,Tci,Cc,hx.params.min_approach)
            set_Q(hx.cold_stream,Qx; mode = :in)
            set_Q(hx.hot_stream,Qx; mode = :out)
        else
            Qx,Tho,Tco = hxeff(Thi,Ch,Tci,Cc,hx.params.ϵ)
            set_Q(hx.cold_stream,Qx; mode = :in)
            set_Q(hx.hot_stream,Qx; mode = :out)
        end
        advance(hx.cold_stream)
        advance(hx.hot_stream)
    end
end

function process(hx::heat_exchanger)
    Qh = get_Q(hx.hot_stream)
    Qc = get_Q(hx.cold_stream)

    if abs(abs(Qh) -  abs(Qc))>1.0
        println("difference in stream loads")
        error("process heat exchanger")
    end
    process(hx.hot_stream)
    process(hx.cold_stream)
    
    hx.Q_hx = Qh

    Thi,Phi,Tho,Pho         = get_TP(hx.hot_stream)
    Tci,Pci,Tco,Pco         = get_TP(hx.cold_stream)

    hx.Ch = get_inlet_duty(hx.hot_stream)
    hx.Cc = get_inlet_duty(hx.cold_stream)
    hx.Tci = Tci
    hx.Tco = Tco
    hx.ΔT_c = Tco - Tci

    hx.Thi = Thi
    hx.Tho = Tho
    hx.ΔT_h = Tho - Thi

    ΔT_1,ΔT_2,ΔT_logmean = lmtd(Thi,Tho,Tci,Tco)
    hx.ΔT_lm = ΔT_logmean

    hx.UA = hx.Q_hx / hx.ΔT_lm

end
#===============================================================#
#     LIQUID COMPONENTS
#       Modeled as 1x cooler, 1x heater
#       gas heater and cooler elements
#           hx_mode = :approach, hot fixed, cold fixed, effectiveness
#           Heat_Modes - Tset , heat load calculated from dT  
#           Qset , Output temperature calculated from prescribed heat load
#===============================================================#
Base.@kwdef mutable struct simple_pump_parameters <: component_parameters
    η::Float64  = 0.9   # isentropic effeciency
end
pump_def_parameters() = simple_pump_parameters(η = 0.9)

Base.@kwdef mutable struct simple_pump <: model_component
    #note - this parameter will compute based on the decribed tempreatutres of the nodes - they should be fixed
    nodes::siso_nodes                                # inlet and outlet
    params::simple_pump_parameters = simple_pump_parameters()

    #across variables
    ΔH::Union{Float64,Int64} = 0       #enthalpy (specific)
    ΔT::Union{Float64,Int64} = 0       #temperature
    ΔP::Union{Float64,Int64} = 0       #pressure change

    Q::Union{Float64,Int64} = 0        #    Q
    W::Union{Float64,Int64} = 0        #     W
    
    #through variables
    mass_flow::Union{Float64,Int64}     =   1
    name::String                        =   "LiquidPump"
end

function init_simple_pump(AB::Vector{<:component_nodes}; name_string = "pump")
    param = pump_def_parameters();
    nodes = siso_nodes(AB)
    return simple_pump(nodes; params = param, ΔH = 0.0, ΔT = 0.0, ΔP = 0.0, Q = 0.0, W = 0.0, mass_flow = 1.0, name = name_string)
end



function show_component(c::comp) where comp<:model_component
    str = String(nameof(typeof(c)))
    @printf "%s \t Q = %.2f, \t W = %.2f, \t ΔT = %.2f , ΔH = %.2f ,ΔP = %.2f \n" str c.W c.W c.ΔT c.ΔH c.ΔP
end
function show_node(n::node) where node <: component_nodes
    ptyp,ttyp = get_node_fix_type(n)
    
    @printf "Node %i \t T [C] = %.2f %s , P [Bar] = %.2f %s \n" n.node_num (n.T-273.15) ttyp (n.P) (ptyp)
end

function show_node(cvec::Vector{<:component_nodes})
    for nod in cvec
        show_node(nod)
    end
end

function show_comp_nodes(c::comp) where comp<:model_component
    inP, inT, oP,oT = get_node_fix_type(c.nodes)
    f1str = (inP,inT)
    f2str = (oP,oT)
    str = c.name
    @printf "%s \t \tNode = %s %i => %i %s  \n" str f1str c.nodes.A.node_num c.nodes.B.node_num f2str
end

function show_comp_nodes(cvec::Vector{comp}) where comp<:model_component
    for c in cvec
        show_comp_nodes(c)
    end
end


function component_dp(c::comp) where comp<:model_component
    if typeof(c) != gas_regenerator                  
        return (c.nodes.B.P - c.nodes.A.P)
    else
        #evaluating regenerator
        dp1 = (c.hot_nodes.B.P - c.hot_nodes.A.P)
        dp2 =(c.cold_nodes.B.P - c.cold_nodes.A.P)
        return dp1 + dp2
    end
end

Base.@kwdef mutable struct component_network_params
    Tmin::Float64 = 650.0     # temperature
    Pmin::Float64 = 80.0      # pressure
    Tmax::Float64 = 900.0     # max allowable temperatures
    working_fluid::Symbol = :helium
end

breeder_network_params() = component_network_params(700 + 273.15, 32, 1273.15, :pbli)
blanket_network_params() = component_network_params(450 + 273.15, 32, 750+273.15, :helium)
divertor_network_params() = component_network_params(350 + 273.15, 32, 650+273.15, :helium)

mutable struct component_network
    params::component_network_params
    nodes::Vector{<:component_nodes}
    elements::Vector{model_component}
end

function match_compression_ratio(Pmin,Nc,Nt;    ic_dp = 0.2, rh_dp = 0.2,   regen_dp = 0.2,    heating_dp = 0.2,   cooling_dp = 0.2)
    Pout(pin,ncomp,rp,dp) = pin*rp^ncomp - dp * sum(rp .^ (1:(ncomp-1)))
    Pout_comp(rpc_x) = Pout(Pmin,Nc,rpc_x,ic_dp)
    ΔP_1 = regen_dp + heating_dp;
    ΔP_2 = regen_dp + cooling_dp;

    fzer(rp1,rp2) = Pout(Pout_comp(rp1)-ΔP_1, Nt, 1/rp2, rh_dp) - ΔP_2 - Pmin
    if Nc > Nt
        rpt = 3.75;
        f_c(rp) = fzer(rp,rpt)
        rpc = find_zero(f_c,3)
    else 
        rpc = 3.75;
        f_t(rp) = fzer(rpc,rp)
        rpt = find_zero(f_t,3)
    end
    return rpt, rpc
end

function brayton_network(Tmax_guess; Tmin = 300.0, Nc=1, Nt=1, Nhx = 3, regen=true, Pmin = 50.0,cycle_fluid = :helium)
    params = component_network_params(Tmin,Pmin,Tmax_guess,cycle_fluid)
    # flow source and initial conditinos
    addl_states = 2;
    states_in_between_comp_and_turb = 2 + Nhx   # compout, regen_out, 
     
    totStates   = 1+ 2*Nc + 2*Nt+ 2 + Nhx
    # totComp     = Nt+Nc+numIC+numRH+addl_comp+2 # +2 for primaryheat and primary cooling

    heating_dp = 0.2*Nhx;
    cooling_dp = 0.2;

    rpt,rpc = match_compression_ratio(Pmin,Nc,Nt)



    # Initializing
    nodes = init_gas_nodes(totStates)

    component_vector = Vector{model_component}(undef,1)
    component_vector[1] = mass_flow_source(nodes[1:2];params = mass_flow_source_paremeter(1.0,cycle_fluid))
    
    compIdx = 1:(2*Nc-1);           # component indexes 1 - last compressor
    idcount = repmat(1:Nc,Nc);      # for labels
    for i in compIdx
        idx = i+1
        node_idx = [idx,idx+1]
        if iseven(idx)
            cparams = comp_params(;rpin = rpc)
            push!(component_vector,gas_compressor(nodes[node_idx];params = cparams, name = "compressor_$(idcount[idx])"))
        else
            push!(component_vector, intercooler(nodes[node_idx],Tmin; name_string =  "intercool_$(idcount[idx])"))
        end
    end

    idx += 1
    node_idx = [idx,idx+1]
    push!(component_vector, init_regen_cold(nodes[node_idx]))

    for i =1:Nhx
        idx += 1
        node_idx = [idx,idx+1]
        push!(component_vector, heat_exchanger_heater(nodes[node_idx]; name_string = "HX_$(i)"))
    end

    turbIdx = 1:(2*Nt-1);           # component indexes 1 - last compressor
    idcount = repmat(1:Nt,Nt);      # for labels
    for i =1:(2*Nt-1)
        idx += 1
        node_idx = [idx,idx+1]
        if isodd(i)
            tparams = comp_params(;rpin = rpc)
            push!(component_vector,gas_turbine(nodes[node_idx];params = tparams, name = "turbine_$(idcount[i])"))
        else
            push!(component_vector, simple_reheater(nodes[node_idx],Tmin; name_string =  "reheat_$(idcount[i])"))
        end
    end

    
    idx += 1
    node_idx = [idx,idx+1]
    push!(component_vector, init_regen_hot(nodes[node_idx]))

    idx += 1
    node_idx = [idx,1]
    push!(component_vector, heat_sink(nodes[node_idx],Tmin,Pmin))

    return component_network(params,nodes,component_vector)

end

function cooling_network(para::component_network_params)
    numNodes = 5
    fluid_type  = gastype(para.working_fluid)
    component_vector = Vector{model_component}(undef,1)

    if fluid_type == :liq
        nodes = init_liq_nodes(numNodes,para.Pmin)
        component_vector[1] = mass_flow_source(nodes[1:2];  params = mass_flow_source_paremeter(1.0,  para.working_fluid))
        idx = [2,3]
        nxt_nodes = nodes[2:3]
        # fix_pressure(nodes[1],para.Pmin)
        # fix_pressure(nodes[2],para.Pmin)
        @show nxt_nodes
        @show typeof(siso_nodes(nxt_nodes))

        push!(component_vector, init_simple_pump(nxt_nodes))
        nxt_nodes = nodes[3:4]
        push!(component_vector, ideal_heat_source(nxt_nodes; Qin = 0.0, ΔP = 12.0, name_string = "Fusion Heat"))
        nxt_nodes =  nodes[4:5]
        push!(component_vector, heat_exchanger_heater(nxt_nodes; Qin = 0.0, ΔP = 12.0, name_string = "To_cycle"))
        idx = [5,1]
        nxt_nodes = nodes[idx]
        push!(component_vector, heat_sink(nxt_nodes, para.Tmin, para.Pmin; name_string = "Sink"))
    else
        nodes = init_gas_nodes(numNodes,para.Pmin)
        fix_pressure(nodes[1],para.Pmin)
        fix_pressure(nodes[2],para.Pmin)
        rpc = (para.Pmin+3.0) / para.Pmin
        component_vector[1] = mass_flow_source(nodes[1:2];params = mass_flow_source_paremeter(1.0,  para.working_fluid))
        nxt_nodes = nodes[2:3]
        push!(component_vector, gas_compressor(nxt_nodes; params = comp_params(;rpin = rpc)))
        nxt_nodes = nodes[3:4]
        push!(component_vector, ideal_heat_source(nxt_nodes; Qin = 0.0, name_string = "Fusion Heat"))
        nxt_nodes =  nodes[4:5]
        push!(component_vector, heat_exchanger_heater(nxt_nodes; Qin = 0.0, name_string = "To_Cycle"))
        idx = [5,1]
        nxt_nodes = nodes[idx]
        push!(component_vector, heat_sink(nxt_nodes, para.Tmin, para.Pmin; name_string = "Sink"))
    end
    return component_network(para,nodes,component_vector)
end

function default_coolant_networks()
    breeder_network     = cooling_network(breeder_network_params())
    divertor_network    = cooling_network(divertor_network_params())
    blanket_network     = cooling_network(blanket_network_params())
    return breeder_network, divertor_network,blanket_network
end


function process(cn::component_network)

    for comp in cn.elements
        process(comp)
    end
end


function advance(cn::component_network)
    for comp in cn.elements
        advance(comp)
    end
end

function show_comp_details(c::model_component) show_component(c::comp) where comp<:model_component
    str = String(nameof(typeof(c)))
    @printf "%s => ΔT = %.2f , ΔH = %.2f ,ΔP = %.2f \n" str,c.ΔT c.ΔH c.ΔP
end

function show_comp_details(cn::component_network)
    for c in cn.elements
        if typeof(c) != [mass_flow_source]
            show_component(c)
        end
    end

end