
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
#           sets mass flow of the 
#   process(c::component_type)  Computes the component variables        - acts ONLY on the comoponent struct
#   check(c::compoennt_type)    Ensure energy balance
#   
#   print funciton
#===============================================================#
#                    EQUATIONS
#   V̄ = velocity
#   v = specific volume
#   V = volume
#
#   e = P/ρ + V̄^2 / 2 + gz  =>  u + ke + pe
#   E = m * e
#   Ė = ṁ * e
#
#   Wflow = pressure * volume = Pv
#   Flow total energy
#        θ = Pv + e
#           =>  (Pv + u) + (V^2/2) + (gz)    
#           =>  h + ke + pe
#           =>  h + V^2 / 2 + gz 
#
#   Conservation of Mass => always applies
#
#       ṁ_in = ṁ_out
#
#   Conservation of energy
#   
#       Ėin     = Q̇in  + Ẇin  + ∑ (ṁ_in  * θ_in)  
#       Ėout    = Q̇out + Ẇout + ∑ (ṁ_out * θ_out)
#
#   Steady State ΔĖsys = 0
#                       Ėin = Ėin
#           (Q̇in - Q̇out) + (Ẇin - Ẇout) = ∑ (ṁ_out * θ_out) - ∑ (ṁ_in  * θ_in) 
#   SIGN CONVENTIONS
#       Q, W 
#       
#    Always true
#       Ėin - Ėout = ΔĖsys
#===============================================================#

using Graphs

abstract type model_component end
abstract type fluid_component     <: model_component end      # functinoal elemen
abstract type source_element      <: model_component end      # mass fow, 
abstract type reservoir           <: model_component end      # heat and work reservoirs
abstract type component_parameters end

##      HIGH LEVEL UTILITIES
###     MASS FLOW SOURCE

Base.@kwdef mutable struct mass_flow_source_paremeter <: component_parameters
    working_fluid::Symbol       = :helium
    mass_flow::Float64          = 1.0
end

mass_flow_source_paremeter(massFlow::Union{Float64,Int64}) = mass_flow_source_paremeter(mass_flow = massFlow)
mass_flow_source_paremeter(massFlow::Union{Float64,Int64}, fluid::Symbol) = mass_flow_source_paremeter(mass_flow = massFlow, working_fluid = fluid)

function print_params(para::mass_flow_source_paremeter)
    @printf "\t Fluid = %s, Mass flow [kg/s] = %.2f \n" String(para.fluid)  para.mass_flow
end

Base.@kwdef mutable struct mass_flow_source <: source_element
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
    name::String = "M-FLOW SRC"
end

function advance(c::mass_flow_source; verbose = false)
    set_outlet_temp(c.nodes,get_inlet_temp(c.nodes))
    set_outlet_press(c.nodes,get_inlet_press(c.nodes))

    if get_working_fluid(c.nodes.A) == :water
        if c.nodes.A.h == 0
            set_h(c.nodes.A,c.nodes.B.h)
        end
    end
    return
end

function process(c::mass_flow_source; verbose = false)
    # advance(c)
    return
end

mass_flow_source(AB::Vector{<:component_nodes}; kw...)                           =   mass_flow_source(nodes = siso_nodes(AB); kw...)
mass_flow_source(AB::Vector{<:component_nodes}, massFlow::Union{Float64,Int64})   =   mass_flow_source(AB; params = mass_flow_source_paremeter(massFlow), mass_flow = massFlow)

#   set_mass_flow = sets the mass flow of a fluid 
function set_mass_flow(c::fluid_component, mflow::Union{Float64,Int64}; use_existing_mflow_as_reference = false) 
    #   if use_existing_mflow_as_reference, then mass flow - 0->1 @ all compoenents
    c.mass_flow = use_existing_mflow_as_reference ? c.mass_flow * mflow : mflow
end

function component_info(c::mass_flow_source; kw...)
    @printf "\t%s\n\t\t mflow = %.2f [kg/s] \t fliud = %s \n" c.name c.mass_flow c.params.working_fluid
end

# ===================  TURBINE AND COMPRESSORS ==============================
#     
#       parameters  rp, η
#       Constructors
#
# =============================================================================

Base.@kwdef mutable struct gas_working_element_parameters <: component_parameters
    rp::Float64 = 3.5   # compression ratio
    η::Float64  = 0.9   # isentropic effeciency
end

gas_circulator_params() = gas_working_element_parameters(1.0,0.9)

function print_params(para::gas_working_element_parameters)
    @printf " \t rp = %.2f, \t η =  %.2f\n " para.rp para.η
end

turb_params(;rpin::Float64) = gas_working_element_parameters(rp=rpin)
comp_params(;rpin::Float64) = gas_working_element_parameters(rp=rpin)

Base.@kwdef mutable struct gas_compressor <: fluid_component
    nodes::siso_nodes   # inlet and outlet
    params::gas_working_element_parameters = gas_working_element_parameters()

    #across variables
    ΔT::Union{Float64,Int64} = 0       #temperature
    ΔP::Union{Float64,Int64} = 0       #pressure change

    
    ΔH::Union{Float64,Int64} = 0       # enthalpy
    ΔU::Union{Float64,Int64} = 0       # internal energy
    ΔS::Union{Float64,Int64} = 0       # entropy

    Q::Union{Float64,Int64} = 0        # Q
    W::Union{Float64,Int64} = 0        # W

    #through variables
    mass_flow::Union{Float64,Int64}  = 1

    name::String = "compressor"
end

function advance(c::gas_compressor; verbose = true)
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

Base.@kwdef mutable struct gas_turbine <: fluid_component
    nodes::siso_nodes   # inlet and outlet
    params::gas_working_element_parameters = gas_working_element_parameters()

    #across variables
    ΔT::Union{Float64,Int64} = 0       #temperature
    ΔP::Union{Float64,Int64} = 0       #pressure change

    
    ΔH::Union{Float64,Int64} = 0       # enthalpy
    ΔU::Union{Float64,Int64} = 0       # internal energy
    ΔS::Union{Float64,Int64} = 0       # entropy
    Q::Union{Float64,Int64} = 0        # Q
    W::Union{Float64,Int64} = 0        # W
    #through variables
    mass_flow::Union{Float64,Int64}  = 1
    name::String = "Turbine"
end

function advance(t::gas_turbine; verbose = false)
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
    inprops,outprops    = calculate_fluid_properties(c.nodes)
    T1,P1,T2,P2         = get_TP(c.nodes)
    ΔT,ΔP,Δu,Δh,Δs      = nodal_changes(c.nodes)
    c.ΔT =      ΔT
    c.ΔP =      ΔP
    c.ΔU =      Δu * c.mass_flow
    c.ΔH =      Δh * c.mass_flow
    c.ΔS =      Δs * c.mass_flow
    c.W =       c.ΔH
end
gas_turbine(AB::Vector{gas_node}; kw...)                    =   gas_turbine(nodes = siso_nodes(AB); kw...)
gas_compressor(AB::Vector{gas_node};kw...)                  =   gas_compressor(nodes=siso_nodes(AB);kw...)

function component_info(c::Union{gas_compressor,gas_turbine}; kw...)
    @printf "\t%s\n\t\t mflow = %.2f [kg/s] \t fliud = %s \t W = %.2f \n" c.name c.mass_flow c.nodes.A.working_fluid c.W/10^6
    @printf "\t\t Tin = %.2f  [K] \t Tout = %.2f \n" get_inlet_temp(c.nodes) get_outlet_temp(c.nodes)
end

# =================== HEATERS==============================
#  All heater and cooler elements
#   Parameters
#       Heat_Modes -    Tset , (Default) heat load calculated from dT  
#                       Qset , Output temperature calculated from prescribed heat load
#                               Use Qset for regenerators
#                       continous => sets outlet to inlet temp, Q = 0
#                       Tdesired  => 
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
    source_mode::Symbol = :env  #env or transfer, env = heat exchange or internal exchange
end

regen_heat_params() = heating_element_parameters(heat_mode = :Qset, pressure_mode = :active, source_mode = :transfer)

function print_params(para::heating_element_parameters)
    @printf "\t\t Heat Mode %s, \t pressure mode: %s, \tdefault ΔP: %.2f \t Source Mode = %s \n" String(para.heat_mode) String(para.pressure_mode) para.default_ΔP String(para.source_mode)
end

Base.@kwdef mutable struct heat_source <: fluid_component
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
    para = heating_element_parameters(heat_mode = :Tset, pressure_mode = :active, default_ΔP = ΔP, source_mode = :env)
    nds = siso_nodes(AB)
    fix_temperature(nds.B,Tout)
    return heat_source(nodes = nds; params = para, name = name_string)
end

function reheater(AB::Vector{<:component_nodes},Tout::Float64; ΔP = 0.2, name_string = "reheat")
    para = heating_element_parameters(heat_mode = :Qset, pressure_mode = :active, default_ΔP = ΔP, source_mode = :transfer)
    nds = siso_nodes(AB)
    return heat_source(nodes = nds; params = para, name = name_string)
end

function simple_reheater(AB::Vector{<:component_nodes},Tout::Float64; ΔP = 0.2, name_string = "reheat")
    para = heating_element_parameters(heat_mode = :Tset, pressure_mode = :active, default_ΔP = ΔP, source_mode = :env)
    nds = siso_nodes(AB)
    return heat_source(nodes = nds; params = para, name = name_string)
end

#   heat sink - very simple component, passive pressrue drop, set outlet conditions
function heat_sink(AB::Vector{<:component_nodes},Tout::Float64,Pout::Float64; name_string = "heat_sink", hmode = :Tset)
    para = heating_element_parameters(heat_mode = hmode, pressure_mode = :passive, default_ΔP = 0.0, source_mode = :env)
    nds = siso_nodes(AB)
    if hmode == :Tset
        fix_temperature(nds.B,Tout)
    end
    fix_pressure(nds.B,Pout)
    return heat_source(nodes = nds; params = para, name = name_string)
end

function simple_condensor(AB::Vector{<:component_nodes}; name_string = "condensor")
    para = heating_element_parameters(heat_mode = :Tset, pressure_mode = :passive, default_ΔP = 0.0, source_mode = :env)
    nds = siso_nodes(AB)
    heat_source(nodes = nds; params = para, name = name_string)
end

function simple_boiler(AB::Vector{<:component_nodes}; name_string = "boiler")
    para = heating_element_parameters(heat_mode = :Tset, pressure_mode = :passive, default_ΔP = 0.0, source_mode = :env)
    nds = siso_nodes(AB)
    heat_source(nodes = nds; params = para, name = name_string)
end

function ideal_heat_source(AB::Vector{<:component_nodes};Qin::Float64=0.0, ΔP = 0.2, name_string = "ideal_source")
    para = heating_element_parameters(heat_mode = :Qset, pressure_mode = :active, default_ΔP = ΔP, source_mode = :env)
    return heat_source(AB; params = para, Q = Qin, name = name_string)
end

ideal_heat_source(AB::Vector{<:component_nodes},Qin::Float64) = heat_source(nodes = siso_nodes(AB); params = heating_element_parameters(heat_mode = :Qset, pressure_mode = :active, default_ΔP = 0.2), Q = Qin)

function heat_exchanger_heater(AB::Vector{<:component_nodes}; Qin = 0.0, ΔP = 0.2, fixed_T::Bool = false, Ti::Float64 = -1.0, To::Float64 = -1.0, name_string = "hx_noname")
    
    # fixed T - means that this side of the heat exchanger has fixed temperature differential
    nds = siso_nodes(AB)
    if fixed_T == false
        # not fixed T, meaning Q is calculated during hx evaluation - which controls temperature distribution , this is the more robust version
        para = heating_element_parameters(heat_mode = :Qset, pressure_mode = :active, default_ΔP = ΔP, source_mode = :transfer)
        return heat_source(nodes = nds; params = para, Q = Qin,name = name_string)

    else
        if Ti != -1
            fix_temperature(nds.A,Ti)
        end
        if To != -1
            fix_temperature(nds.B,To)
        end

        para = heating_element_parameters(heat_mode = :Tset, pressure_mode = :active, default_ΔP = ΔP, source_mode = :transfer)
        return heat_source(nodes = nds; params = para,name = name_string)
    end
end

function fix_outlet_temp(hs::heat_source,Tmin::Real)
    hs.params.heat_mode = :Tset
    fix_temperature(hs.nodes.B,Tmin)
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
        Q = abs(Q)  # heat addition
    else
        Q = -abs(Q)  # heat loss
    end
    q_type = check_q_type(hs)
    if q_type != :Qset && overide == false
        println("Attempt to change heat load manually, change kw args to overide = true")
        return
    else
        hs.Q = Q
    end
end

function advance(hs::heat_source; verbose = false)
    #calculate the next state of the heat source

    #set node mass flow
    # set_mass_flow(hs.nodes,hs.mass_flow)
    ΔT,ΔP,Δu,Δh,Δs      = nodal_changes(hs.nodes)

    if check_q_type(hs) == :Qset
        # set heat load, calcualating temperature rise
        inprops,outprops = calculate_fluid_properties(hs.nodes)
        inprops.cp
        cp_ave = (inprops.cp + outprops.cp)/2
        if isnan(cp_ave)
            @show cp_ave = inprops.cp
        end
        
        Ti = get_inlet_temp(hs.nodes)
        Q = hs.Q                                # (+) means heat ADDITION mdot cp (T2-T1) = Q1 - Q2,  T2 = Q/mcp+T1
        mdot_cp = get_mass_flow(hs.nodes.A) * cp_ave
        if get_working_fluid(hs.nodes.A) == :water
            mdot_cp = get_inlet_duty(hs.nodes)
        end
        Tout = Q/mdot_cp + Ti
        set_outlet_temp(hs.nodes,Tout)

    elseif check_q_type(hs) == :cont
        set_outlet_temp(hs.nodes, get_inlet_temp(hs.nodes))
    end

    if check_p_type(hs) != :passive
        p1 = get_inlet_press(hs.nodes)
        set_outlet_press(hs.nodes, p1 - hs.params.default_ΔP)
    end

    if verbose == true
        println(hs.name)
    end
end

function process(hs::heat_source)
    ptype = check_p_type(hs)
    qtype = check_q_type(hs)
    inprops,outprops    = calculate_fluid_properties(hs.nodes)
    T1,P1,T2,P2         = get_TP(hs.nodes)
    mass_flow           = get_mass_flow(hs.nodes.A)
    ΔT,ΔP,Δu,Δh,Δs      = nodal_changes(hs.nodes)
    hs.ΔT =      ΔT
    hs.ΔP =      ΔP
    hs.ΔU =      Δu * mass_flow
    hs.ΔH =      Δh * mass_flow
    hs.ΔS =      Δs * mass_flow
    if qtype == :Tset
        # Heat in
        mode = :in
        if ΔT <= 0 || Δh < 0
            mode =  :out
        end
        set_Q(hs, hs.ΔH; mode = mode, overide = true)
    end
end

function component_info(c::heat_source,; kw...)
    @printf "\t%s\n\t\t mflow = %.2f [kg/s] \t fliud = %s \t Q = %.2f\n" c.name c.nodes.A.mass_flow c.nodes.A.working_fluid c.Q/10^6
    @printf "\t\t Duty = %.2f [MW]\t Tin = %.2f  [K] \t Tout = %.2f ΔH = %2s\n" get_inlet_duty(c.nodes)/10^6 get_inlet_temp(c.nodes) get_outlet_temp(c.nodes) c.ΔH
    print_params(c.params)
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
    hx_mode = :eff          #Other options, :eff = effectiveness NTU method, :pinch
    ϵ::Float64 = .95
    min_approach::Float64 = 10.0
    UA = 0.0
end

get_hx_mode(hpar::heat_exchanger_parameters) = hpar.hx_mode

Base.@kwdef mutable struct heat_exchanger <: fluid_component
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

function heat_exchanger(cold_heater::heat_source,hot_heater::heat_source; kw...)
    return heat_exchanger(cold_stream = cold_heater,hot_stream = hot_heater;kw...)
end

function heat_exchanger(hx_vec::Vector{model_component};kw...)
    heat_exchanger(hx_vec[1],hx_vec[2])
end

function advance(hx::heat_exchanger; kw...)
    cold_type   = check_q_type(hx.cold_stream)
    hot_type    = check_q_type(hx.hot_stream)
    # advance(hx.cold_stream)
    # advance(hx.hot_stream)
    if cold_type == hot_type == :Tset
        # this cant be the case unless
        println("Heat exchanger joined to incompattible heating elements, both with prescribed temperature differentials")
        error("Advance heat exchagner error")
    end

    if cold_type == :Tset
        advance(hx.cold_stream)
        # process(hx.cold_stream)
        @assert(hot_type == :Qset)
        set_Q(hx.hot_stream,hx.cold_stream.Q; mode = :out)
        # advance(hx.hot_stream)
        # process(hx.hot_stream)
    elseif hot_type == :Tset
        advance(hx.hot_stream)
        # process(hx.hot_stream)
        @assert(cold_type == :Qset)
        set_Q(hx.cold_stream,hx.hot_stream.Q; mode = :in)
        # advance(hx.cold_stream)
        # process(hx.cold_stream)
    else
        # both are Qset
        Cc = get_inlet_duty(hx.cold_stream)
        Ch = get_inlet_duty(hx.hot_stream)

        Tci = get_inlet_temp(hx.cold_stream)
        Thi = get_inlet_temp(hx.hot_stream)
        mode = get_hx_mode(hx)
        # display(mode)
        if mode == :approach
            Qx,Tho,Tco = hx_pinch(Thi,Ch,Tci,Cc,hx.params.min_approach)
            set_Q(hx.cold_stream,Qx; mode = :in)
            set_Q(hx.hot_stream,Qx; mode = :out)
            println("Approach_mode")
        else
            if Thi >= Tci
                Qx,Tho,Tco = hxeff(Thi,Ch,Tci,Cc,hx.params.ϵ; verbose = false)
                set_Q(hx.cold_stream, Qx; mode = :in)
                set_Q(hx.hot_stream, Qx; mode = :out)
                if Tco > Thi
                    Qx,Tho,Tco = hxeff(Tci,Ch,Tci,Cc,hx.params.ϵ; verbose = true)
                end
                # @printf "Que Bueno \t\t  Tci = %.2f Thi = %.2f\n" Tci Thi
                # @printf "Que Bien: \t\t  Tco = %.2f Tho = %.2f\n" Tco Tho
            else
                Qx,Tco,Tho = hxeff(Tci,Cc,Thi,Ch,hx.params.ϵ; verbose = false)
                set_Q(hx.cold_stream, abs(Qx); mode = :out)
                set_Q(hx.hot_stream, abs(Qx); mode = :in)
                # println("$(hx.name) Tci >= Thi")
                # @printf "TROUBLE Tco = %.2f Tho = %.2f\n" Tco Tho
            end
        end
        # advance(hx.cold_stream)
        # advance(hx.hot_stream)
    end
end

function process(hx::heat_exchanger; verbose = false)
    Qh = get_Q(hx.hot_stream)
    Qc = get_Q(hx.cold_stream)

    if abs(abs(Qh) -  abs(Qc))>1.0
        println("difference in stream loads $(abs(abs(Qh) -  abs(Qc)))")
        # error("process heat exchanger")
    end
    # advance(hx.cold_stream)
    # advance(hx.hot_stream)
    # process(hx.hot_stream)
    # process(hx.cold_stream)
    
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
    # if Tco > Thi
    #     println("$(hx.name) Tco > Thi")
    # end

end

function show_details(hx::heat_exchanger)
    process(hx)
    @printf "%s HEAT EXCHANGER\n" hx.name
    @printf "\tHOT SIDE\t %s \t Thi = %.2f, \t Tho = %.2f, Ch = %.2f,, Qh = %.2f \n" hx.hot_stream.nodes.A.working_fluid hx.Thi hx.Tho hx.Ch/10^6 get_Q(hx.hot_stream)/10^6
    @printf "\tCLD SIDE\t %s \t Tci = %.2f, \t Tco = %.2f, Cc = %.2f,  Qc = %.2f\n"     hx.cold_stream.nodes.A.working_fluid hx.Tci  hx.Tco hx.Cc/10^6  get_Q(hx.cold_stream)/10^6
    @printf "\tTOTALS: Qhx [MW] = %.2f \n" hx.Q_hx/10^6
end

function component_info(hx::heat_exchanger; kw...)
    @printf "\t %s HEAT EXCHANGER \t MODE = %s \n" hx.name get_hx_mode(hx)
    print_params(hx.hot_stream.params)
    print_params(hx.cold_stream.params)
    @printf "\t\tHOT SIDE\t %s \t Thi = %.2f, \t Tho = %.2f, Ch = %.2f,, Qh = %.2f \n" hx.hot_stream.nodes.A.working_fluid hx.Thi hx.Tho hx.Ch/10^6 get_Q(hx.hot_stream)/10^6
    @printf "\t\tCLD SIDE\t %s \t Tci = %.2f, \t Tco = %.2f, Cc = %.2f,  Qc = %.2f\n"  hx.cold_stream.nodes.A.working_fluid hx.Tci  hx.Tco hx.Cc/10^6  get_Q(hx.cold_stream)/10^6
    @printf "\t\tTOTALS: Qhx [MW] = %.2f \n" hx.Q_hx/10^6
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

Base.@kwdef mutable struct simple_pump <: fluid_component
    #note - this parameter will compute based on the decribed tempreatutres of the nodes - they should be fixed
    nodes::siso_nodes                                # inlet and outlet
    params::simple_pump_parameters = pump_def_parameters()

    #across variables
    ΔT::Union{Float64,Int64} = 0       # temperature
    ΔP::Union{Float64,Int64} = 0       # pressure change

    ΔH::Union{Float64,Int64} = 0       # enthalpy
    ΔU::Union{Float64,Int64} = 0       # internal energy
    ΔS::Union{Float64,Int64} = 0       # entropy

    Q::Union{Float64,Int64} = 0 
    W::Union{Float64,Int64} = 0
    
    #through variables
    mass_flow::Union{Float64,Int64}     =   1
    name::String                        =   "LiquidPump"
end
show_params(p::simple_pump) = println(p.params.η)

function init_simple_pump(AB::Vector{<:component_nodes}; name_string = "Circ. Pump")
    param = pump_def_parameters();
    nodes = siso_nodes(AB)
    return simple_pump(nodes, param, 0.0,  0.0, 0.0, 0.0, 0.0, 0.0,0.0,1.0,  name_string)
end

function advance(p::simple_pump;verbose = false)
    # calculating us Δh = cΔT + vΔP
    # initializing
    v1,P1,v2,P2 = get_vP(p.nodes)
    P1 = bar2pascal(P1)
    P2 = bar2pascal(P2)
    v_ave = (v1+v2)/2
    fluid_props_a, fluid_props_b = calculate_fluid_properties(p.nodes)
    cp_ave = (fluid_props_a.cp + fluid_props_b.cp)/2

    # pressure change
    ΔP  = (P2-P1)
    vΔP = (v_ave * ΔP)              # reversible work
    Δh  = (vΔP)/p.params.η  # actual specific work
    
    dt = (Δh - vΔP) / cp_ave;
    T1 = get_inlet_temp(p.nodes)
    T2 = T1 + dt;
    set_outlet_temp(p.nodes,T2)
end

function process(hs::simple_pump)
    inprops,outprops    = calculate_fluid_properties(hs.nodes)
    T1,P1,T2,P2         = get_TP(hs.nodes)
    ΔT,ΔP,Δu,Δh,Δs      = nodal_changes(hs.nodes)
    hs.ΔT =      ΔT
    hs.ΔP =      ΔP
    hs.ΔU =      Δu * hs.mass_flow
    hs.ΔH =      Δh * hs.mass_flow
    hs.ΔS =      Δs * hs.mass_flow
    hs.W =       hs.ΔH
end

function component_info(c::simple_pump; kw...)
    @printf "\t%s\n\t\t => Tin = %.2f Tout = %.2f W = %.2f \n" c.name c.nodes.A.T c.nodes.B.T c.W/10^6
end


#===============================================================#
#    Water components
#===============================================================#

#       PUMP
Base.@kwdef mutable struct h20_pump_parameters <: component_parameters
    η::Float64  = 0.9   # isentropic effeciency
end

Base.@kwdef mutable struct h20_pump <: fluid_component
    #note - this parameter will compute based on the decribed tempreatutres of the nodes - they should be fixed
    nodes::siso_nodes                                # inlet and outlet
    params::h20_pump_parameters = h20_pump_parameters()

    #across variables
    ΔT::Union{Float64,Int64} = 0       # temperature
    ΔP::Union{Float64,Int64} = 0       # pressure change

    ΔH::Union{Float64,Int64} = 0       # enthalpy
    ΔU::Union{Float64,Int64} = 0       # internal energy
    ΔS::Union{Float64,Int64} = 0       # entropy

    Q::Union{Float64,Int64} = 0 
    W::Union{Float64,Int64} = 0
    
    #through variables
    mass_flow::Union{Float64,Int64}     =   1
    name::String                        =   "LiquidPump"
end

function initialize(p::h20_pump, Pinlet::Real, Poutlet::Real)
    p1 = Pinlet;
    set_x(p.nodes.A,0.0)
    set_x(p.nodes.B,0.0)
    lprops = initLiquidWater(p1)

    set_inlet_temp(p.nodes,lprops.T)
    set_inlet_press(p.nodes,Pinlet)
    set_outlet_press(p.nodes,Poutlet)
    fix_temperature(p.nodes.A)

    fix_pressure(p.nodes.A)
    fix_pressure(p.nodes.B)
    set_h(p.nodes.A,lprops.h)
end

function advance(p::h20_pump; verbose = false)
    # calculating us Δh = cΔT + vΔP
    # initializing
    v1,P1,v2,P2 = get_vP(p.nodes)
    steam_in = calculate_fluid_properties(p.nodes.A)   # isentropic
    v_ave = (v1 + v2) / 2
    
    vΔP = 1e5 * (P2 - P1) * v_ave
    Δh  = vΔP / p.params.η
    h2  = (steam_in.h + Δh )/1e3
    T2  = T_ph(P2, h2)
    x_out = x_ph(P2,h2)
    set_outlet_temp(p.nodes,T2 + 273.15)
    set_x(p.nodes.B,x_out)
    set_h(p.nodes.B,h2 * 1e3)

    if verbose
        @printf "pump evaluataion %s \n" p.name
        @printf "\t h1 = %.2f, \t h2a = %.2f \n" steam_in.h/1e3 h2
    end
end

function process(hs::h20_pump)
    inprops,outprops    = calculate_fluid_properties(hs.nodes)
    T1,P1,T2,P2         = get_TP(hs.nodes)
    ΔT,ΔP,Δu,Δh,Δs      = nodal_changes(hs.nodes)
    hs.ΔT =      ΔT
    hs.ΔP =      ΔP
    hs.ΔU =      Δu * hs.mass_flow
    hs.ΔH =      Δh * hs.mass_flow
    hs.ΔS =      Δs * hs.mass_flow
    hs.W =       hs.ΔH
end


function component_info(c::h20_pump; kw...)
    @printf "\t%s\n\t\t Tin = %.2f Pin = %.2f Tout = %.2f Pout = %.2f W = %.2f \t ṁ = %2s \n" c.name c.nodes.A.T c.nodes.A.P c.nodes.B.T c.nodes.B.P c.W/10^6 get_mass_flow(c.nodes.A)
end

#       TURBINE
Base.@kwdef mutable struct h20_turbine_parameters <: component_parameters
    η::Float64  = 0.9   # isentropic effeciency
end

# siso turbine
Base.@kwdef mutable struct h20_turbine <: fluid_component
    #note - this parameter will compute based on the decribed tempreatutres of the nodes - they should be fixed
    nodes::siso_nodes                                # inlet and outlet
    params::h20_turbine_parameters = h20_turbine_parameters()

    #across variables
    ΔT::Union{Float64,Int64} = 0       # temperature
    ΔP::Union{Float64,Int64} = 0       # pressure change

    ΔH::Union{Float64,Int64} = 0       # enthalpy
    ΔU::Union{Float64,Int64} = 0       # internal energy
    ΔS::Union{Float64,Int64} = 0       # entropy

    Q::Union{Float64,Int64} = 0 
    W::Union{Float64,Int64} = 0
    
    #through variables
    mass_flow::Union{Float64,Int64}     =   1
    name::String                        =   "Steam Turbine"
end

function initialize(p::h20_turbine, Pinlet::Real, Poutlet::Real)
    set_inlet_press(p.nodes,Pinlet)
    set_outlet_press(p.nodes,Poutlet)
end

function advance(p::h20_turbine; verbose = false)
    # calculating us Δh = cΔT + vΔP
    # initializing
    v1,P1,v2,P2 = get_vP(p.nodes)
    steam_in = calculate_fluid_properties(p.nodes.A)   # isentropic
    set_h(p.nodes.A, steam_in.h)
    Pout = P2
    h2s = h_ps(P2, steam_in.s / 1e3)
    w   = p.params.η * (steam_in.h/1e3 - h2s)
    h   = steam_in.h/1e3 - w
    s2  = s_ph(Pout, h)
    x2  = x_ps(Pout, s2)
    T   = T_ps(Pout, s2)
    set_temperature(p.nodes.B,T+273.15)
    set_x(p.nodes.B,x2)
    set_h(p.nodes.B.h * 1e3)
    if verbose
        @printf "Turbine evaluataion %s \n" p.name
        @printf "\t h1 = %.2f, \t h2s = %.2f \t h2a = %.2f \n" steam_in.h/1e3 h2s h
    end

end

function process(hs::h20_turbine)
    inprops,outprops    = calculate_fluid_properties(hs.nodes)
    T1,P1,T2,P2         = get_TP(hs.nodes)
    ΔT,ΔP,Δu,Δh,Δs      = nodal_changes(hs.nodes)
    hs.ΔT =      ΔT
    hs.ΔP =      ΔP
    hs.ΔU =      Δu * hs.mass_flow
    hs.ΔH =      Δh * hs.mass_flow
    hs.ΔS =      Δs * hs.mass_flow
    hs.W =       hs.ΔH
end


function component_info(c::h20_turbine; kw...)
    @printf "%s\n\t\t\t Tin = %.2f Tout = %.2f \n" c.name c.nodes.A.T c.nodes.B.T
end

#simo turbine
Base.@kwdef mutable struct h20_simo_turbine <: fluid_component
    #note - this parameter will compute based on the decribed tempreatutres of the nodes - they should be fixed
    nodes::simo_nodes                                # inlet and outlet
    params::h20_turbine_parameters = h20_turbine_parameters()

    #across variables
    ΔT::Vector{Union{Float64,Int64}} = []       # temperature
    ΔP::Vector{Union{Float64,Int64}} = []       # pressure change

    ΔH::Vector{Union{Float64,Int64}} = []       # enthalpy
    ΔU::Vector{Union{Float64,Int64}} = []       # internal energy
    ΔS::Vector{Union{Float64,Int64}} = []       # entropy

    Q::Union{Float64,Int64} = 0.0 
    W::Union{Float64,Int64} = 0.0
    
    #through variables
    mass_flow::Union{Float64,Int64}     =   0.0
    name::String                        =   "Steam Turbine"
end

function initialize(p::h20_simo_turbine, Pinlet::Float64, Poutlet::Vector{Float64})
    set_pressure(p.nodes.A,Pinlet)
    fix_pressure(p.nodes.A,Pinlet)
    tot_p = length(Poutlet)

    for (idx,P) in enumerate(Poutlet)
        set_pressure(p.nodes.B[idx],P)
        push!(p.ΔT,0)
        push!(p.ΔP,P-Pinlet)
        push!(p.ΔU,0)
        push!(p.ΔH,0)
        push!(p.ΔS,0)
        set_mass_flow(p.nodes.B[idx],1/tot_p)
    end
end

function advance(p::h20_simo_turbine; verbose = false)
    # calculating us Δh = cΔT + vΔP
    # initializing
    P1 = get_pressure(p.nodes.A)

    steam_in = calculate_fluid_properties(p.nodes.A)   # isentropic
    
    s1 = steam_in.s
    set_h(p.nodes.A, steam_in.h)

    Pout = get_outlet_press(p.nodes)
    if verbose
        @printf "Turbine evaluataion %s \n" p.name
    end
    for (idx,P2) in enumerate(Pout)
        # @show steam_in.s
        # @show steam_in.x
        steam_in.T
        h2s = h_ps(P2, steam_in.s / 1e3)                # kJ / kg
        if isnan(h2s)
            @printf "\t Nan XSTEAM query h2s for P = %.2f s = %.2f\n" P2 steam_in.s/1e3
        end
        w   = p.params.η * (steam_in.h/1e3 - h2s)
        h   = steam_in.h/1e3 - w

        s2  = s_ph(P2, h)
        x2  = x_ps(P2, s2)
        T   = T_ps(P2, s2)

        set_temperature(p.nodes.B[idx],T+273.15)
        set_x(p.nodes.B[idx],x2)
        set_h(p.nodes.B[idx],h*1e3)
        p.nodes.B[idx].T
        if verbose
            @printf "\t h1 = %.2f, \t h2s = %.2f \t h2a = %.2f \n" steam_in.h/1e3 h2s h
        end
    end

end

function process(turb::h20_simo_turbine)
    nout = length(turb.nodes.B)
    inprops = calculate_fluid_properties(turb.nodes.A)
    turb.W = 0
    turb.Q = 0

    for i = 1:nout
        temp_nodes = siso_nodes(turb.nodes.A,turb.nodes.B[i])
        ΔT,ΔP,Δu,Δh,Δs      = nodal_changes(temp_nodes)
        turb.ΔT[i] =      ΔT
        turb.ΔP[i] =      ΔP
        turb.ΔU[i] =      Δu * turb.nodes.B[i].mass_flow
        turb.ΔH[i] =      Δh * turb.nodes.B[i].mass_flow
        turb.ΔS[i] =      Δs * turb.nodes.B[i].mass_flow
        turb.W  +=       turb.ΔH[i]
    end
end

function component_info(c::h20_simo_turbine; kw...)
    @printf "\t%s\n\t\t Pin = %.2f Tin = %.2f [K] W = %.2f\n" c.name c.nodes.A.P c.nodes.A.T c.W/10^6
    for (i,nod) in enumerate(c.nodes.B)
        @printf "\t\t Pout = %.2f Tout = %.2f [K] h = %.2f \t ṁ = %.2f \n" nod.P nod.T nod.h/10^6 nod.mass_flow
    end
end
## WATER COMPONENTS
#   Closed feedwater heater
#   Open feedwater heater
#   Mixing Chamer
#   Boiler
#   Condensor

Base.@kwdef mutable struct h20_open_feedwater_heater <: fluid_component
    #note - this parameter will compute based on the decribed tempreatutres of the nodes - they should be fixed
    nodes::miso_nodes                                # inlet and outlet
    yflow::Vector{<:component_nodes}                  # nodes associated with y
    xflow::Vector{<:component_nodes}                  # nodes associated wtih y-1


    #across variables
    ΔT::Vector{Union{Float64,Int64}} = []       # temperature
    ΔP::Vector{Union{Float64,Int64}} = []       # pressure change

    ΔH::Vector{Union{Float64,Int64}} = []       # enthalpy
    ΔU::Vector{Union{Float64,Int64}} = []       # internal energy
    ΔS::Vector{Union{Float64,Int64}} = []       # entropy

    Q::Union{Float64,Int64} = 0.0 
    W::Union{Float64,Int64} = 0.0
    
    #through variables
    mass_flow::Union{Float64,Int64}     =   0.0
    name::String                        =   "OFWH"
end

function advance(ofw::h20_open_feedwater_heater; verbose = false)
    h1 = ofw.nodes.A[1].h # flow y
    h2 = ofw.nodes.A[2].h # flow 1-y
    h3 = ofw.nodes.B.h
    #   computing mass fraction
    dh1     =  h1 - h2
    dh2     =  h3 - h2

    # @assert (dh1 >= 0)
    # @assert (dh2 >= 0)

    y = abs(dh2 / dh1)

    for node in ofw.yflow
        set_mass_flow(node, y * ofw.mass_flow;    use_existing_mflow_as_reference = false)
    end
    for node in ofw.xflow
        set_mass_flow(node, (1-y) * ofw.mass_flow;    use_existing_mflow_as_reference = false)
    end
end

function process(ofw::h20_open_feedwater_heater)
    advance(ofw)
    return
end

function component_info(c::h20_open_feedwater_heater; kw...)
    @printf "\t %s OFW - Pout = %.2f Tout = %.2f [K] hout = %.2f \n" c.name c.nodes.B.P c.nodes.B.T c.nodes.B.h/10^6
    for (i,nod) in enumerate(c.nodes.A)
        @printf "\t\t\t Pin = %.2f Tin = %.2f [K] hin = %.2f \n" nod.P nod.T nod.h
    end
end

###########################
# ============================================================================#
#                                                                             #
#                              NETWORK ANAYLYS
#
# # =============================================================================

# Base.@kwdef mutable struct component_network_params
#     Tmin::Float64 = 650.0     # temperature
#     Pmin::Float64 = 80.0      # pressure
#     Tmax::Float64 = 900.0     # max allowable temperatures
#     working_fluid::Symbol = :helium
# end

# mutable struct component_network
#     params::component_network_params
#     nodes::Vector{<:component_nodes}
#     elements::Vector{<:model_component}
#     mass_flow::Float64
#     name::String
# end

# mutable struct system_network
#     cycles::Vector{component_network}
#     heat_exchangers::Vector{<:heat_exchanger}
#     name::String
#     eval_order
# end

# function show_params(cn::component_network_params)
#     @printf "NETWORK: Tmin = %.2f, \t Pmin = %.2f,\t Tmax = %.2f\t working fluid = %s\n" cn.Tmin cn.Pmin cn.Tmax cn.working_fluid
# end

# function show_params(cn::component_network)
#     show_params(cn.params)
# end

# breeder_network_params()        = component_network_params(700 + 273.15, 32, 1273.15, :pbli)
# blanket_network_params()        = component_network_params(450 + 273.15, 80, 750+273.15, :helium)
# divertor_network_params()       = component_network_params(350 + 273.15, 80, 650+273.15, :helium)

# breeder_network_params2()       = component_network_params(550 + 273.15, 32, 750+273.15, :pbli)
# blanket_network_params2()       = component_network_params(450 + 273.15, 80, 450+273.15, :helium)
# divertor_network_params2()      = component_network_params(350 + 273.15, 80, 350+273.15, :helium)
# extraction_network_params2()    = component_network_params(300 + 273.15, 80, 700+273.15, :helium)

# function set_name(cn::component_network,in_string::String)
#     cn.name = in_string
# end

# component_network(params::component_network_params,nodes::Vector{<:component_nodes},elements::Vector{<:model_component}) = component_network(params,nodes,elements,1.0,"noname")

# function node_count(cn::component_network)
#     return length(cn.nodes)
# end

# function elem_count(cn::component_network)
#     return length(cn.elements)
# end

# function network2graph(cycle_network::component_network; verbose = false)
#     num_edge    = elem_count(cycle_network) # number of elements
#     num_vert    = node_count(cycle_network) # number of nodes (vertices)
#     g           = DiGraph(num_vert)
#     nms =[]
#     edgelabel_dict = Dict()
#     edgelabel_mat = Array{String}(undef, num_vert, num_vert)

#     for el in cycle_network.elements
#         push!(nms,el.name)
#         if typeof(el.nodes) ==  siso_nodes
#             i = el.nodes.A.node_num
#             j = el.nodes.B.node_num
#             add_edge!(g,i,j)
#             edgelabel_mat[i, j] = edgelabel_dict[(i, j)] = el.name
#             verbose ? println("Node added :  $(i) to $(j)  <=> $(el.name) ") : nothing
#         elseif typeof(el.nodes) == miso_nodes
#             #   mult input, single output
#             j = el.nodes.B.node_num
#             for in_node in el.nodes.A
#                 i = in_node.node_num
#                 add_edge!(g,i,j)
#                 edgelabel_mat[i, j] = edgelabel_dict[(i, j)] = el.name
#                 verbose ? println("Node added :  $(i) to $(j)  <=> $(el.name) ") : nothing
#             end
#         elseif typeof(el.nodes) == simo_nodes
#             i = el.nodes.A.node_num
#             for OUT_node in el.nodes.B
#                 j = OUT_node.node_num
#                 add_edge!(g,i,j)
#                 edgelabel_mat[i, j] = edgelabel_dict[(i, j)] = el.name
#                 verbose ? println("Node added :  $(i) to $(j)  <=> $(el.name) ") : nothing
#             end
#         end
#     end
#     return g,nms, edgelabel_mat, edgelabel_dict

# end

# function network2graph2(cycle_network::component_network; verbose = false)
#     num_elem    = elem_count(cycle_network) # number of elements
#     num_vert    = node_count(cycle_network) # number of nodes (vertices)
#     g           = DiGraph(num_elem)
#     nms =[]
#     edgelabel_dict = Dict()
#     edgelabel_mat = Array{String}(undef, num_vert, num_vert)

#     for (idx1,el1) in enumerate(cycle_network.elements)
#         # add_vertex!(g)
#         push!(nms,el1.name)
#         inn = inlet_nodes(el1.nodes)
#         for node_idx in inn
#             for (idx2,el2) in  enumerate(cycle_network.elements)
#                 onn = outlet_nodes(el2.nodes)
#                 for src_node in onn
#                     if src_node == node_idx
#                         add_edge!(g,idx2,idx1)
#                         edgelabel_mat[idx2, idx1] = edgelabel_dict[(idx2, idx1)] = string(node_idx)
#                         verbose ? println("EDGE added :  $(el2.name) to $(el1.name)  <=> $(src_node) ") : nothing
#                     end
#                 end
#             end
#         end
#     end

#     return g,nms, edgelabel_mat, edgelabel_dict

# end

# function component_info(cn::component_network)
#     @printf "CYCLE OBJECT COMPONENT INFO => %s\n" cn.name
#     for el in cn.elements
#         component_info(el)
#     end
# end

# # ============================================================================#
# #                                                                             #
# #                             Saved cases
# #
# # =============================================================================

# function match_compression_ratio(Pmin,Nc,Nt;    ic_dp = 0.2, rh_dp = 0.2,   regen_dp = 0.2,    heating_dp = 0.2,   cooling_dp = 0.2)
#     Pout(pin,ncomp,rp,dp) = pin*rp^ncomp - dp * sum(rp .^ (1:(ncomp-1)))
#     Pout_comp(rpc_x) = Pout(Pmin,Nc,rpc_x,ic_dp)
#     ΔP_1 = regen_dp + heating_dp;
#     ΔP_2 = regen_dp + cooling_dp;

#     fzer(rp1,rp2) = Pout(Pout_comp(rp1)-ΔP_1, Nt, 1/rp2, rh_dp) - ΔP_2 - Pmin
#     if Nc > Nt
#         rpt = 3.75;
#         f_c(rp) = fzer(rp,rpt)
#         rpc = find_zero(f_c,3)
#     else 
#         rpc = 3.75;
#         f_t(rp) = fzer(rpc,rp)
#         rpt = find_zero(f_t,3)
#     end
#     return rpt, rpc
# end

# function brayton_network(Tmax_guess; Tmin = 300.0, Nc=1, Nt=1, Nhx = 3, regen=true, Pmin = 50.0,cycle_fluid = :helium)
#     params = component_network_params(Tmin,Pmin,Tmax_guess,cycle_fluid)
#     # flow source and initial conditinos
#     addl_states = 2;
#     states_in_between_comp_and_turb = 2 + Nhx   # compout, regen_out, 
     
#     totStates   = 1+ 2*Nc + 2*Nt+ 1 + Nhx
#     # totComp     = Nt+Nc+numIC+numRH+addl_comp+2 # +2 for primaryheat and primary cooling

#     heating_dp = 0.2*Nhx;
#     cooling_dp = 0.2;

#     rpt,rpc = match_compression_ratio(Pmin,Nc,Nt)

#     # Initializing
#     nodes = init_gas_nodes(totStates)

#     component_vector = Vector{model_component}(undef,1)
#     component_vector[1] = mass_flow_source(nodes[1:2];params = mass_flow_source_paremeter(1.0,cycle_fluid))
    
#     compIdx = 1:(2*Nc-1);           # component indexes 1 - last compressor
#     idcount = repmat(1:Nc,Nc);      # for labels
#     for i in compIdx
#         idx = i+1
#         node_idx = [idx,idx+1]
#         if iseven(idx)
#             cparams = comp_params(;rpin = rpc)
#             push!(component_vector,gas_compressor(nodes[node_idx];params = cparams, name = "compressor_$(idcount[idx])"))
#         else
#             push!(component_vector, intercooler(nodes[node_idx],Tmin; name_string =  "intercool_$(idcount[idx])"))
#         end
#     end

#     idx += 1
#     node_idx = [idx,idx+1]
#     push!(component_vector, init_regen_cold(nodes[node_idx]))

#     regen_cold_idx = length(component_vector)

#     for i =1:Nhx
#         idx += 1
#         node_idx = [idx,idx+1]
#         push!(component_vector, heat_exchanger_heater(nodes[node_idx]; name_string = "Heat_EX_$(i)"))
#     end
#     hx_idx = (regen_cold_idx+1:length(component_vector))

#     turbIdx = 1:(2*Nt-1);           # component indexes 1 - last compressor
#     idcount = repmat(1:Nt,Nt);      # for labels
#     for i =1:(2*Nt-1)
#         idx += 1
#         node_idx = [idx,idx+1]
#         if isodd(i)
#             tparams = turb_params(;rpin = rpt)
#             push!(component_vector,gas_turbine(nodes[node_idx];params = tparams, name = "turbine_$(idcount[i])"))
#         else
#             push!(component_vector, simple_reheater(nodes[node_idx],Tmin; name_string =  "reheat_$(idcount[i])"))
#         end
#     end

    
#     idx += 1
#     node_idx = [idx,idx+1]
#     push!(component_vector, init_regen_hot(nodes[node_idx]))
#     regen_hot_idx = length(component_vector)
#     idx += 1
#     node_idx = [idx,1]
#     push!(component_vector, heat_sink(nodes[node_idx],Tmin,Pmin))

#     return component_network(params,nodes,component_vector,1.0,"BRAYTON CYCLE"), [regen_cold_idx, regen_hot_idx], hx_idx

# end

# function cooling_network(para::component_network_params; Qin_s = 0.0, verbose = false, dynamic = false)
#     numNodes = 5
#     fluid_type          = gastype(para.working_fluid)
#     component_vector    = Vector{model_component}(undef,1)
#     if fluid_type == :liq
#         nodes                   = init_liq_nodes(numNodes,para.Pmin)
#         component_vector[1]     = mass_flow_source(nodes[1:2];  params = mass_flow_source_paremeter(1.0,  para.working_fluid))
#         idx = [2,3]
#         nxt_nodes = nodes[2:3]
#         fix_pressure(nodes[1],para.Pmin)
#         fix_pressure(nodes[2],para.Pmin)
#         push!(component_vector, init_simple_pump(nxt_nodes)); 
#         nxt_nodes = nodes[3:4]
#         push!(component_vector, ideal_heat_source(nxt_nodes; Qin = 0.0, ΔP = 6.0, name_string = "Fusion Heat"));

#         nxt_nodes =  nodes[4:5]
#         push!(component_vector, heat_exchanger_heater(nxt_nodes; Qin = 0.0, ΔP = 6.0, name_string = "To_cycle"))
#         idx = [5,1]
#         nxt_nodes = nodes[idx]
#         if dynamic == false
#             push!(component_vector, heat_sink(nxt_nodes, para.Tmin, para.Pmin; hmode = :Tset, name_string = "Heat_Sink"))
#         else
#             push!(component_vector, heat_sink(nxt_nodes, para.Tmin, para.Pmin; hmode = :cont, name_string = "Heat_Sink"))
#         end
#     else
#         nodes = init_gas_nodes(numNodes,para.Pmin)
#         fix_pressure(nodes[1],para.Pmin)
#         fix_pressure(nodes[2],para.Pmin)
#         rpc = (para.Pmin+3.0) / para.Pmin
#         component_vector[1] = mass_flow_source(nodes[1:2];params = mass_flow_source_paremeter(1.0,  para.working_fluid))
#         nxt_nodes = nodes[2:3]
#         push!(component_vector, gas_compressor(nxt_nodes; params = comp_params(;rpin = rpc), name = "Circulator"))
#         nxt_nodes = nodes[3:4]
#         push!(component_vector, ideal_heat_source(nxt_nodes; Qin =  0.0, name_string = "Fusion Heat"))
#         nxt_nodes =  nodes[4:5]
#         push!(component_vector, heat_exchanger_heater(nxt_nodes; Qin = 0.0, name_string = "To_Cycle"))
#         idx = [5,1]
#         nxt_nodes = nodes[idx]

#         if dynamic == false
#             push!(component_vector, heat_sink(nxt_nodes, para.Tmin, para.Pmin; hmode = :Tset, name_string = "Heat_Sink"))
#         else
#             push!(component_vector, heat_sink(nxt_nodes, para.Tmin, para.Pmin; hmode = :cont, name_string = "Heat_Sink"))
#         end
#     end
#     set_Q(component_vector[3],Qin_s; mode = :in)
#     for n in nodes
#         overide_set_temperature(n,para.Tmin)        # initalizing tempreatures
#     end
#     return component_network(para,nodes,component_vector)
# end

# function extraction_loop(para::component_network_params; dynamic = false)
#     numNodes = 8
#     component_vector = Vector{model_component}(undef,1)

#     nodes = init_gas_nodes(numNodes,para.Pmin)
#     fix_pressure(nodes[1],para.Pmin)
#     fix_pressure(nodes[2],para.Pmin)
#     rpc = (para.Pmin+3.0) / para.Pmin
#     component_vector[1] = mass_flow_source(nodes[1:2];params = mass_flow_source_paremeter(1.0,  para.working_fluid))
#     nxt_nodes = nodes[2:3]
#     push!(component_vector, gas_compressor(nxt_nodes; params = comp_params(;rpin = rpc), name = "Circulator"))
#     nxt_nodes = nodes[3:4]
#     push!(component_vector, heat_exchanger_heater(nxt_nodes; Qin = 0.0, name_string = "HX1"))
#     nxt_nodes =  nodes[4:5]
#     push!(component_vector, heat_exchanger_heater(nxt_nodes; Qin = 0.0, name_string = "HX2"))
#     nxt_nodes =  nodes[5:6]
#     push!(component_vector, heat_exchanger_heater(nxt_nodes; Qin = 0.0, name_string = "HX3"))
#     nxt_nodes =  nodes[6:7]
#     push!(component_vector, heat_exchanger_heater(nxt_nodes; Qin = 0.0, name_string = "Boiler_hx"))
#     idx = [7,1]

#     if dynamic == false
#         push!(component_vector, heat_sink(nodes[idx], para.Tmin, para.Pmin; name_string = "Heat_Sink"))
#     else
#         push!(component_vector, heat_sink(nxt_nodes, para.Tmin, para.Pmin; hmode = :cont, name_string = "Heat_Sink"))
#     end

#     for n in nodes
#         overide_set_temperature(n,para.Tmin)        # initalizing tempreatures
#     end
#     return component_network(para,nodes,component_vector)
# end

# function default_coolant_networks()
#     breeder_network     = cooling_network(breeder_network_params();Qin_s =  740.0e6)
#     divertor_network    = cooling_network(divertor_network_params();Qin_s =  140.0e6)
#     blanket_network     = cooling_network(blanket_network_params();Qin_s =  100.0e6)
#     breeder_network.name =  "breeder_network" 
#     divertor_network.name = "divertor_network" 
#     blanket_network.name =  "blanket_network" 
#     return breeder_network, divertor_network,blanket_network
# end

# function simple_rankine(;Pmax,Pmid,Pmin,TmaxK)
#     nodes = init_water_nodes(8)
#     component_vector = Vector{model_component}(undef,7)

#     mflow_nodes     = siso_nodes(nodes[1:2])
#     turbine_nodes   = simo_nodes(nodes[2],nodes[3:4])
#     condensor_nodes = [nodes[4],nodes[5]]
#     pump1_nodes     = siso_nodes(nodes[5:6])
#     ofw_nodes       = miso_nodes([nodes[3],nodes[6]],nodes[7])
#     pump2_nodes     = siso_nodes(nodes[7:8])

#     boiler_nodes = [nodes[8],nodes[1]]


#     set_temperature(nodes[2],   TmaxK)
#     set_temperature(nodes[1],  TmaxK)
#     # refMax = initSuperHeatedWater(TmaxK,Pmax)
#     # set_h(nodes[1],refMax.h)

#     yflow = [nodes[3]]
#     xflow = nodes[4:6]

#     component_vector[1] = mass_flow_source(nodes = mflow_nodes; params = mass_flow_source_paremeter(working_fluid = :water, mass_flow = 100.0))
#     component_vector[5] = h20_simo_turbine(nodes = turbine_nodes; params = h20_turbine_parameters(1.0))
#     component_vector[6] = simple_condensor(condensor_nodes)

#     component_vector[3] = h20_pump(nodes = pump1_nodes; params = h20_pump_parameters(1.0), name = "Pump1")
#     component_vector[4] = h20_pump(nodes = pump2_nodes; params = h20_pump_parameters(1.0), name = "Pump2")

#     component_vector[7] = h20_open_feedwater_heater(nodes = ofw_nodes, xflow = xflow, yflow = yflow)

#     component_vector[2] = heat_exchanger_heater(boiler_nodes; name_string = "boiler_qin", fixed_T = true, To = TmaxK)

#     #initializing
#     initialize(component_vector[3], Pmin,Pmid)
#     initialize(component_vector[4], Pmid,Pmax)
#     initialize(component_vector[5], Pmax,[Pmid,Pmin])

#     para = component_network_params(Tmin = 273, Pmin = Pmin, Tmax = TmaxK, working_fluid = :water)
#     cm = component_network(para,nodes,component_vector,1.0,"Rankine")

#     advance(cm)
#     process(cm)
#     advance(cm)
#     process(cm)
#     return cm
# end

# function rankine_with_loops()
#     rankine = simple_rankine(; Pmax = 80.0, Pmid = 12.0, Pmin = 0.1, TmaxK = 600+273.15)
#     breeder_circuit    =    cooling_network(breeder_network_params2(); Qin_s =  740.0e6, dynamic = false)
#     divertor_circuit    =   cooling_network(divertor_network_params2();Qin_s =  140.0e6, dynamic = false)
#     blanket_circuit    =    cooling_network(blanket_network_params2(); Qin_s =  100.0e6, dynamic = false)

#     breeder_circuit.name =  "breeder_network" 
#     divertor_circuit.name = "divertor_network" 
#     blanket_circuit.name =  "blanket_network" 

#     ensure_pressure_match(breeder_circuit)
#     ensure_pressure_match(divertor_circuit)
#     ensure_pressure_match(blanket_circuit)

#     change_and_update_mass_flow(breeder_circuit, 20000.0)
#     change_and_update_mass_flow(divertor_circuit,1000.0)
#     change_and_update_mass_flow(blanket_circuit,1000.0)

#     extraction_network = extraction_loop(extraction_network_params2(), dynamic = false)
#     extraction_network.name = "Intermediate_loop"

#     ensure_pressure_match(extraction_network)
#     change_and_update_mass_flow(extraction_network,1000.0)
    
#     change_and_update_mass_flow(rankine,300.0)
#     #hx in 3 4 5
#     hx = [3,4,5]

#     # fix_outlet_temp(breeder_circuit.elements[5],breeder_circuit.params.Tmin)

#     divertor_heat_exchanger = heat_exchanger(extraction_network.elements[hx[1]],divertor_circuit.elements[4]; name = "div_hx")
#     blanket_heat_exchanger  = heat_exchanger(extraction_network.elements[hx[2]],blanket_circuit.elements[4]; name = "blk_hx")
#     breeder_heat_exchanger  = heat_exchanger(extraction_network.elements[hx[3]],breeder_circuit.elements[4]; name = "breed_hx")
#     boiler_heat_exchanger   = heat_exchanger(rankine.elements[2],   extraction_network.elements[6]; name = "boiler_hx")

#     cycle_vec   = [divertor_circuit,    blanket_circuit, breeder_circuit, extraction_network, rankine]
#     hx_vec      = [divertor_heat_exchanger, blanket_heat_exchanger, breeder_heat_exchanger, boiler_heat_exchanger]

#     eval_order = Dict()
#     eval_order[1] = divertor_circuit
#     eval_order[2] = divertor_heat_exchanger
#     eval_order[3] = extraction_network
#     eval_order[4] = blanket_circuit
#     eval_order[5] = blanket_heat_exchanger
#     eval_order[6] = extraction_network
#     eval_order[7] = breeder_circuit
#     eval_order[8] = breeder_heat_exchanger
#     eval_order[9] = extraction_network
#     eval_order[10] = boiler_heat_exchanger
#     eval_order[11] = rankine

#     return system_network(cycle_vec, hx_vec, "Rankine with Intermediate loops", eval_order)
# end


# # http://docs.juliaplots.org/latest/graphrecipes/examples/
# # ============================================================================#
# #                                                                             #
# #                              Printing
# #
# # =============================================================================


# # utilities for setup
# function ensure_pressure_match(cn::component_network)
#     # only use this function for the cooling loops
#     Pmin = cn.params.Pmin
#     ΔP_0 = 0
#     for c in cn.elements
#         if typeof(c) != mass_flow_source
#             if :default_ΔP ∈ fieldnames(typeof(c.params))
#                 ΔP_0 += c.params.default_ΔP
#             end
#         end
#     end

#     fluid_type  = gastype(cn.params.working_fluid)
#     if fluid_type == :gas
#         cn.elements[2].params.rp = (Pmin + ΔP_0)/Pmin
#     else
#         fix_pressure(cn.elements[2].nodes.B,(Pmin + ΔP_0))
#     end
# end

# function change_source_flow(cn::component_network,new_flow::Float64)
#     for c in cn.elements
#         if typeof(c) == mass_flow_source
#             c.mass_flow = new_flow
#             cn.mass_flow = new_flow
#             break
#         end
#     end
# end

# function update_mass_flow(cn::component_network)
#     for n in cn.nodes
#         n.mass_flow = cn.mass_flow
#     end
#     for c in cn.elements
#         if typeof(c) != mass_flow_source
#             c.mass_flow = cn.mass_flow
#         end
#     end
# end

# function change_and_update_mass_flow(cn::component_network,mass_flow::Float64)
#     change_source_flow(cn,mass_flow)
#     update_mass_flow(cn)
# end

# function set_initial_temperatures(cn::component_network,Tmin::Float64)
#     for node in cn.nodes
#         set_temperature(node,Tmin)
#     end
# end

# function set_initial_temperatures(cn::component_network)
#     Tmin = cn.params.Tmin
#     for node in cn.nodes
#         set_temperature(node,Tmin)
#     end
# end

# function set_initial_pressure(cn::component_network,Pmin::Float64)
#     for node in cn.nodes
#         set_pressure(node,Pmin)
#     end
# end

# function set_initial_pressure(cn::component_network)
#     Pmin = cn.params.Pmin
#     for node in cn.nodes
#         set_pressure(node,Pmin)
#     end
# end

# function process(cn::component_network; kw...)

#     for comp in cn.elements
#         process(comp; kw...)
#     end
# end

# function advance(cn::component_network; kw...)
#     for comp in cn.elements
#         advance(comp; kw...)
#     end
# end

# #####################################PRINTING AND CHECK
# # CHECK FUNCTIONS
# function node_enthalpies(node::component_nodes)
#     props = calculate_fluid_properties(node)
#     @printf "\t %s h = %.2f \n" node.node_num node.h/1e3
# end

# function node_enthalpies(cm::component_network)
#     for n in cm.nodes
#         node_enthalpies(n)
#     end
# end

# function show_node_simple(n::node) where node <: component_nodes
#     @printf "\t# %i, \t T = %.2f, \t P = %.2f, \t ṁ = %.2f, \t fluid = %s \n"    n.node_num n.T n.P n.mass_flow n.working_fluid
# end

# function show_node_simple(cn::component_network)
#     @printf "NETWORK NODES - %s \n Working fluid: %s\n Tmin = %.2f\n Pmin: %.2f \n" cn.name cn.params.working_fluid cn.params.Tmin cn.params.Pmin
#     for n in cn.nodes
#         show_node_simple(n)
#     end
# end

# function show_component(c::comp) where comp<:fluid_component
#     str = String(nameof(typeof(c)))
#     @printf "%s => ΔT = %.2f ,\t ΔP = %.2f ,\t ΔH = %.2f ,\t ΔQ = %.2f,\t ṁ = %.2f  \n" str c.ΔT c.ΔP c.ΔH c.Q c.mass_flow
# end

# function show_comp_details(cn::component_network)
#     for c in cn.elements
#         if typeof(c) != mass_flow_source
#             show_component(c)
#         end
#     end
# end

# function check_energy_balance(cn::component_network; sum_only = false, show_math = false)
#     Hsum = 0
#     @printf "%s \n" cn.name
#     for c in cn.elements
#         if typeof(c) != mass_flow_source
#             if !sum_only
#                 if show_math
#                     @printf "\t %s  \t ΔH = H_total + H_element = %.2f + %.2f = %.2f \n" c.name Hsum/10^6 c.ΔH/10^6 (Hsum + c.ΔH)/10^6
#                 else
#                     @printf "\t %s\t ΔH = %.2f\n" c.name c.ΔH
#                 end
#             end
#             Hsum += c.ΔH
#         end
#     end
#     @printf "\t TOTAL: %s SUM OF ENTHALPY = %.2f\n" cn.name Hsum
# end

# function check_entropy_balance(cn::component_network; sum_only = false)
#     Ssum = 0
#     # @printf "%s \n" cn.name
#     for c in cn.elements
#         if typeof(c) != mass_flow_source
#             if !sum_only
#                 @printf "\t %s\t ΔS = %.2f\n" c.name c.ΔS
#             end
#             Ssum += c.ΔS
#         end
#     end
#     @printf "\t TOTAL: %s SUM OF ENTROPY = %.2f\n" cn.name Ssum
# end

# function show_node(n::node) where node <: component_nodes
#     ptyp,ttyp = get_node_fix_type(n)
#     @printf "Node %i \t T [K] = %.2f %s , P [Bar] = %.2f %s \n" n.node_num (n.T) ttyp (n.P) (ptyp)
# end

# function show_node(cvec::Vector{<:component_nodes})
#     for nod in cvec
#         show_node(nod)
#     end
# end

# function show_comp_nodes(c::comp) where comp<:fluid_component
#     inP, inT, oP,oT = get_node_fix_type(c.nodes)
#     f1str = (inP,inT)
#     f2str = (oP,oT)
#     str = c.name
#     @printf "%s \t \tNode = %s %i => %i %s  \n" str f1str c.nodes.A.node_num c.nodes.B.node_num f2str
# end

# function show_comp_nodes(cvec::Vector{comp}) where comp<: fluid_component
#     for c in cvec
#         show_comp_nodes(c)
#     end
# end

# function show_comp_nodes(cn::component_network)
#     for c in cn.elements
#         show_comp_nodes(c)
#     end
# end

# function component_dp(c::comp) where comp<: fluid_component
#     if typeof(c) != gas_regenerator                  
#         return (c.nodes.B.P - c.nodes.A.P)
#     else
#         #evaluating regenerator
#         dp1 = (c.hot_nodes.B.P - c.hot_nodes.A.P)
#         dp2 =(c.cold_nodes.B.P - c.cold_nodes.A.P)
#         return dp1 + dp2
#     end
# end

# function show_node(cn::component_network)
#     println(cn.name)
#     for nod in cn.nodes
#         show_node(nod)
#     end
# end

# function nodal_vs_element_H(cn::component_network)
#     Hsum = 0
#     @printf "%s \n" cn.name
#     for c in cn.elements
#         if typeof(c) != mass_flow_source
#             ΔT,ΔP,Δu,Δh,Δs = nodal_changes(c.nodes)
#             @printf "\t%s\t (Element) ΔH = %.2f,\t(Nodal) ΔH = %.2f\n" c.name c.ΔH/10^6 Δh*c.mass_flow/10^6
#             # dh_simple = dh_nodes(c.nodes)
#             # @printf "\t%s\t (Element) ΔH = %.2f,\t(Nodal) ΔH = %.2f,\t (Nodal 2) ΔH = %.2f\n" c.name c.ΔH/10^6 Δh*c.mass_flow/10^6 dh_simple*c.mass_flow/10^6
#         end
#     end
# end

# function element_io_temp(cn::component_network)
#     @printf "%s \n" cn.name
#     for c in cn.elements
#         @printf "\t %s, \t Tin = %.2f, Tout = %.2f \t Pin = %.2f, Pout = %.2f\n" c.name get_inlet_temp(c.nodes) get_outlet_temp(c.nodes) get_inlet_press(c.nodes) get_outlet_press(c.nodes)
#     end
# end  

# function element_io_any(cn::component_network,io_var::Symbol)
#     @printf "%s \n" cn.name
#     for c in cn.elements
#         inp,outp = calculate_fluid_properties(c.nodes);
#         @printf "\t %s, \t %s in = %.2f, %s out = %.2f,\t Difference = %.2f \n" c.name String(io_var) getproperty(inp,io_var) String(io_var)  getproperty(outp,io_var) getproperty(outp,io_var)-getproperty(inp,io_var)
#     end
# end

# function element_QW_simple(cn::component_network; show_math = false)
# for c in cn.elements
#     @printf "%s \t\t OUTPUT\t Qin = %.2f \t  Win = %.2f \n" c.name c.Q/10^6 c.W/10^6
# end
# end

# function element_QW(cn::component_network; show_math = false)
#     @printf "%s \n" cn.name
#     sumH = 0
#     sumW = 0
#     sumQ = 0
#     Qin = 0; Qout =
#     Wout = 0; Win = 0;
#     for c in cn.elements
#         if show_math == true
#             @printf "\t %s \t\t Q_net + W_net = %.3f + %.3f = %.3f \t ΔH = %.3f \n " c.name c.Q/10^6 c.W/10^6 (c.Q/10^6 + c.W/10^6) c.ΔH/10^6
#         else
#             @printf "\t %s \t\t Q_net + W_net = %e \t ΔH = %e \n " c.name (c.Q + c.W) c.ΔH
#         end
#         sumQ += c.Q
#         sumW += c.W
#         sumH += c.ΔH
#         c.Q > 0 ? Qin += abs(c.Q) : Qout += abs(c.Q)
#         c.W < 0 ? Wout += abs(c.W) : Win += abs(c.W)
#     end
#     @printf "\t TOTAL \tQ = %.2f \t  W = %.2f \t ΔH = %.2f \n " sumQ/10^6 sumW/10^6 sumH/10^6
#     @printf "\t INPUT \tQin = %.2f \t  Win = %.2f \n" Qin/10^6 Win/10^6
#     @printf "\t OUTPUT\tQout = %.2f \t  Wout = %.2f \n" Qout/10^6 Wout/10^6
#     @printf "\t\t THERMAL EFFECIENCY η = %.2f " sumW/Qin 
# end

# function calculate_circuit_net(cn::component_network;verbose = false)
#     Qin = 0
#     Win = 0
#     Wout = 0
#     Qout = 0
#     # Q_hx = ["", 0]
#     Q_hx = []
#     for c in cn.elements
#         if typeof(c) == heat_source         #   heat source
#             mode = c.params.source_mode     #   external heating
#             if mode == :env 
#                 c.Q > 0 ? Qin += abs(c.Q) : Qout += abs(c.Q)
#             else    # heat exchanger
#                 # Q_hx =vcat(Q_hx, [c.name,c.Q])
#                 push!(Q_hx,c.Q)
#             end 
#         else
#             c.W < 0 ? Wout += abs(c.W) : Win += abs(c.W)
#         end
#     end

#     if verbose
#         println(cn.name)
#         QQ = (Qin/10^6 - Qout/10^6 + sum(Q_hx)/10^6 + Win/10^6 - Wout/10^6)
#         @printf "\t Qin = %.2f MW \t Qout = %.2f MW \t Qhx = %.2f \t Win = %.2f MW \t Wout = %.2f MW \n " Qin/10^6 Qout/10^6 sum(Q_hx)/10^6 Win/10^6 Wout/10^6
#         @printf "\t\t TOTAL = %.2f \n" QQ
#     end
#     return Qin, Qout, Win, Wout, Q_hx
# end

