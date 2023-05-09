#
#
#
# const fluid_opts = [:pbli, :water, :helium]
using Printf, XSteam
mutable struct gas_prop_obj 
    cp::Union{Float64,Int64}    # specific heat
    cv::Union{Float64,Int64}    # specific heat constant volume
    gam::Union{Float64,Int64}   # ratio of specific heats
    v::Union{Float64,Int64}     # specific volume
    R::Union{Float64,Int64}     # gas constant
    s::Union{Float64,Int64}
    h::Union{Float64,Int64}
end

mutable struct liq_prop_obj
    P::Real     #   
    T::Real     # 
    cp::Real    # specific heat
    v::Real     # specific volume
    h::Real     # enthalpy (kJ/kg)
    u::Real     # specific energy
    s::Real     # entropy (kJ/kg)
    x::Real     # vapor fraction
end

mutable struct fluid_prop_obj
    P::Real     # Pressure Bar
    T::Real     # Temperature (K)
    cp::Real    # specific heat (J/kg*K)
    cv::Real    # cv (J/kg*K)
    gam::Real   # Ratio of specific heats (unitless)
    v::Real     # specific volume (m3/kg)
    h::Real     # Specific enthalpy (J/kg)
    u::Real     # specific energy
    s::Real     # entropy (J/kg)
    x::Real     # vapor fraction = 0 for gasses
end

# pbli_fluid_obj(Tin::Union{Float64,Int64},Pin::Union{Float64,Int64},cpin::Union{Float64,Int64},vin::Union{Float64,Int64}) = liq_prop_obj(P=Pin,T=Tin,cp=cpin,v=vin,h = -1,s = -1, x = 0.0)

function gas_properties(gas::Symbol;P=1,T=300)
    #default pressure and temperature (1 bar, 300K)
    gas     = Symbol(lowercase(String(gas)))
    R_dict  = Dict(:helium => 2076.9,   :air => 287.0)  #"J/kg/K"
    cp_dict = Dict(:helium => 5192.6,   :air => 1005.0) #"J/kg/K"
    Pref = 1; #atm
    Tref = 300;
    try
        cp = cp_dict[gas]    # constant pressure specific heat
        R  = R_dict[gas]     # Gas Constant
        cv = cp-R            # constant volume specific heat
        v = R*T/(P*1e5)
        k = cp/cv
        sref = cp * log(T/Tref) - R*log(P/Pref)  
        return gas_prop_obj(cp, cv, k, v, R, sref,cp * T)
    catch e
        println("uknown fluid type in dependent_gas_props()")
        display(e)
    end
end
# liq_prop_obj(P=Pin,T=Tin,cp=cpin,v=vin,h = -1,s = -1, x = 0.0)

function PbLi_properties(Tin::Union{Float64,Int64},Pin::Union{Float64,Int64}; verbose = false)
    if Tin>1615.15
        println("Exceeded lead lithium boiling point $(Tin)")
    elseif Tin<453 && verbose
        println("Warning - lead lithium temperature below melting point $(Tin)")
    end

    # data from [2]
    cpin   = (0.195 - 9.116e-6 .* Tin) .* 1000  # J/kgK
    vin    = 1/(10520.35 - 1.19051 .* Tin)      # m3/kg
    
    cp_ref = (0.195 - 9.116e-6 .* Tin) .* 1000
    # specific entropy of lead lithium, with reference to 508 K
    a = 0.195*1e3
    b = 9.11e-6*1000
    s_pbli      = a*log(Tin) - b*Tin - 1210.316002277804
    u(T)        = 195 * T - .0009116/2 * T^2                # integral of cp
    return liq_prop_obj(Pin, Tin, cpin, vin, -1.0, u(Tin),  s_pbli, 0.0)
end

# superheated water (above dome)
function initSuperHeatedWater(P::Real,T::Real)
    #   assumes vapor fraction = 1, i.e. fluid is 100% superheated 
    cp = Cp_pT(P,T-273.15);
    h = h_pT(P,T-273.15);
    v = v_ph(P,h)
    s = s_ph(P,h)
    u = u_ph(P,h)
    if isnan(s)
        @printf "\t Nan XSTEAM query h2s for T = %.2f P = %.2f h = %.2f\n" T-273.15 P h
    end
    return liq_prop_obj(P,T,cp*1e3,v,h*1e3,u*1e3,s*1e3,1);
end

# any water
function initWater_PT(p::Real,T::Real,x::Real)
    #Tin in Kelvin, Pin in Bar
    Tref = T - 273.15;
    cp = Cp_pT(p,Tref)
    h = h_px(p,x)
    v = v_ph(p,h)
    s = s_ph(p,h)
    u = u_pT(p,Tref)
    x = x_ph(p,h)
    return liq_prop_obj(p,T,cp*1e3,v,h*1e3,u*1e3,s*1e3,x)
end

"""
    initLiquidWater(P::Real)

    Initializes saturated water properties
"""
function initLiquidWater(P::Real)
    cp = CpV_p(P);
    h = hL_p(P);
    v = vL_p(P);
    s = s_ph(P,h);
    x = 0.0
    T = Tsat_p(P);
    u = u_pT(P,T)
    T = T+273.15
    return liq_prop_obj(P,T,cp*1e3,v,h*1e3,u*1e3,s*1e3,x);
end

function liq_properties(liq::Symbol; P = 1, T = 500, x=-1.0)
    liq = enforce_lowercase(liq)
    if liq == :pbli
        return PbLi_properties(T,P)
    elseif liq == :water
        if x != 1.0 && x != 0.0
            return initWater_PT(P,T,x)
        elseif x == 1
            return initSuperHeatedWater(P,T)
        elseif x == 0
            return initLiquidWater(P)
        end
    end

end

function gastype(working_fluid::Symbol)

    if working_fluid in [:air, :helium] 
        return :gas
    else
        return :liq
    end
end


