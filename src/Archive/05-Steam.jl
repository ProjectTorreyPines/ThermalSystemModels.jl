# module FluidProperties
    using Unitful, XSteam, Printf

    const HE = "helium"
    const AIR = "air"
    const PBLI = "pbli"
    const WATER = "water"
    function fluid_properties(name::AbstractString; kw...)
        if name == HE
            return helium_properties(kw...)
        elseif name == AIR
            return air_properties(kw...)
        elseif name == PBLI
            return pbli_properties(kw...)
        elseif name == WATER
            return water_properties(kw...)
        else
            error("Unsupported fluid name")
        end
    end


    # object to hold fluid properties
    struct FluidProp{T<:Real}
        cp::T      # constant pressure specific heat [J/kg/K]
        cv::T      # constant volume specific heat [J/kg/K]
        k::T       # specific heat ratio [-]
        v::T       # specific volume [m³/kg]
        R::T       # gas constant [J/kg/K]
        s::T       # specific entropy [J/kg/K]
        h::T       # specific enthalpy [J/kg/K]
        x::T       # quality [-]
    end

    """
        gas_properties(gas::Symbol; P = 1u"bar", T = 300u"K")

        Return the fluid properties of a gas at a given pressure and temperature.

        Args:
        - `gas`: Symbol representing the gas (e.g., :air, :helium)
        - `P`: pressure [Pa] (default: 1 bar)
        - `T`: temperature [K] (default: 300 K)
    """
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
    air_properties(kw...)       = gas_properties(:helium;kw...)
    helium_properties(kw...)    = gas_properties(:helium;kw...)
    
    """
        PbLi_properties(Tin::Union{Float64,Int64}, Pin::Union{Float64,Int64};verbose = false)

        Return the fluid properties of lead lithium at a given temperature and pressure.

        Args:
        - `Tin`: temperature [K]
        - `Pin`: pressure [Pa]
        - `verbose`: boolean to print warning when Tin is below melting point (default: false)
    """
    function PbLi_properties(Tin::Union{Float64, Int64}, Pin::Union{Float64, Int64};
                            verbose = false)
        if Tin > 1615.15
            println("Exceeded lead lithium boiling point $(Tin)")
        elseif Tin < 453 && verbose
            println("Warning - lead lithium temperature below melting point $(Tin)")
        end
        # Pin = bar
        # data from [2]
        cpin = (0.195 - 9.116e-6 .* Tin) .* 1000  # J/kgK
        vin = 1 / (10520.35 - 1.19051 .* Tin)      # m³/kg

        # specific entropy of lead lithium, with reference to 508 K
        a = 0.195*1e3
        b = 9.11e-6*1000
        s_pbli      = a*log(Tin) - b*Tin - 1210.316002277804
        u(T)        = 195 * T - .0009116/2 * T^2                # integral of cp
        #                                                   h ~ cpT + vP
        return FluidProp(cpin, cpin, 1.0, vin, 1.0, s_pbli, cp * T + vin * (Pin * 1e5) , 0.0)
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
    export fluid_properties, FluidProp
# end

# mutable struct steamObj
#     P::Real # Pressure (Bar)
#     T::Real # Temperature (C)
#     h::Real # enthalpy (kJ/kg)
#     s::Real # entropy (kJ/kg)
#     x::Real # vapor fraction
# end


# """
#     initSuperHeatedWater(P::Real)

#     Initializes superheatedwater properties
# """
# function initSuperHeatedWater(P::Real,T::Real)
#   #  assumes vapor fraction = 1, i.e. fluid is 100% superheated 
#     h = h_pT(P,T);
#     s = s_pT(P,T)
#     return steamObj(P,T,h,s,1);
# end

# """
#     initLiquidWater(P::Real)

#     Initializes saturated water properties
# """
# function initLiquidWater(P::Real)
#     h = h_px(P,0);
#     s = s_ph(P,h);
#     x = x_ph(P,h);
#     T = Tsat_p(P);
#     return steamObj(P,T,h,s,x);
# end

# """
#     steam_initPrePump(P::Vector{Real}) or (P::Real)

#     Initializes saturated water properties for a given pressure at the inlet to a pump, if there are multiple pumps supply the pressrues as a
# """
# function h2o_initPrePump(P::Vector{<:Real})
#     #assumes vapor fraction = 0, i.e. fluid is 100% liquid
#     npts = length(P)
#     if npts==0
#         throw(error("Cant do that"))
#     else
#         i = 1;
#         steams=Vector{steamObj}(undef,npts)
#         for p in P
#             steams[i] = initLiquidWater(P[i])
#             i = i+1;
#         end
#         return steams
#     end
# end


# function h2o_initPrePump(P::Real)
#     #assumes vapor fraction = 0, i.e. fluid is 100% liquid
#     return initLiquidWater(P)
# end


# """
#     steam_initPreTurb(P::Vector{Real}) or (P::Real)

#     Initializes saturated water properties for a given pressure at the inlet to a pump, if there are multiple pumps supply the pressrues as a
# """
# function h2o_initPreTurb(P::Real,T::Real)
#     #assumes vapor fraction = 1, i.e. fluid is 100% gas
#     return initSuperHeatedWater(P,T)
# end

# function h2o_initPreTurb(P::Vector{<:Real},T::Vector{<:Real})
#     npts = length(P)
#     if npts==0
#         throw(error("Cant do that"))
#     else
#         i = 1;
#         steams=Vector{steamObj}(undef,npts)
#         for p in P
#             steams[i] = initSuperHeatedWater(P[i],T[i])
#             i = i+1;
#         end
#         return steams
#     end
# end

# """
#     convert_pressure(val::Real,from_unit::Symbol,to_unit::Symbol)

#     convert pressure units, optional to specify to_unit, default = bar because that is what XSTeam uses
# """
# function convert_pressure(val::Real,from_unit::Symbol,to_unit::Symbol)
#     # default = SI
#     # converting to lowercase to avoid errors
#     from_unit = Symbol(lowercase(string(from_unit)))
#     to_unit = Symbol(lowercase(string(to_unit)))

#     # unit vector (pun)
#     units = (pa = 1.0, kpa = 1.0e3, mpa = 1.0e6, gpa = 1.0e9, psi = 6.8947572932e3, ksi = 6.8947572932e3 * 1000 ,bar = 100000.0, atm = 101325.0, inh2o = 249.082);

#     inUnit = units[from_unit]
#     outUnit = units[to_unit]
#     return val/(outUnit/inUnit)
# end
# function convert_pressure(val::Real,from_unit::Symbol)
#     # default = SI
#     # converting to lowercase to avoid errors
#     from_unit = Symbol(lowercase(string(from_unit)))
#     to_unit = :bar

#     # unit vector (pun)
#     units = (pa = 1.0, kpa = 1.0e3, mpa = 1.0e6, gpa = 1.0e9, psi = 6.8947572932e3, ksi = 6.8947572932e3 * 1000 ,bar = 100000.0, atm = 101325.0, inh2o = 249.082);

#     inUnit = units[from_unit]
#     outUnit = units[to_unit]
#     return val/(outUnit/inUnit)
# end

# """
#     steam_closed_feedwater_heater(steamKnown::steamObj,P_unkown::Real)

#     Analyses closed feedwater heater
# """
# function h2o_closed_feedwater_heater(steamKnown::steamObj,P_unknown::Real)
#     #   Steam_OFW CLOSED FEEDWATER HEATER
#     #                _______________
#     #               |               |
#     #           -> ___/\/\/\/\/\/\______ Steam_known = outlet
#     #               |               |
#     #               |_______________|
#     t_unkown = steamKnown.T;
#     h = h_pT(P_unknown,t_unkown);
#     s = s_pT(P_unknown,t_unkown);
#     x =  x_ph(P_unknown,h);
#     return steamObj(P_unknown,t_unkown,h,s,x);   #outlet steam
# end

# """
#     h2o_MixingChamber(y, steam::Vector{steamObj})

#     Analyses Mixing chamber
# """
# function h2o_MixingChamber(y, steam::Vector{steamObj})
#     # steam is a variable number of arguments, steam[1] corresponds to mass fraction y
#     # and steam[2] corresponds to (1-y)
    
#     h_out = (1-y)*steam[2].h + y*steam[1].h
#     Pout = steam[1].P
#     s = s_ph(Pout, h_out)
#     T = T_ph(Pout, h_out)
#     x = x_ps(Pout, s)
#     steamOut = steamObj(Pout, T,  h_out, s, x)
#     return steamOut
# end

# """
#     h2o_pump(Pou::Real, steam_in::steamObj, etaP)

#     Analyses pump
# """

# function h2o_pump(Pout::Real, steam_in::steamObj, etaP::Float64)
#     Pout = Pout
#     h2s = h_ps(Pout, steam_in.s)
#     w = (h2s - steam_in.h) / etaP
#     h2 = w + steam_in.h
#     s2 = s_ph(Pout, h2)
#     T2 = T_ps(Pout, s2)
#     steam = steamObj(Pout,T2,h2,s2,0.0)
#     return steam, w
# end

# function h2o_pump(Pout_pump::Vector{<:Real}, steam_in::Vector{steamObj}, etaP::Vector{Float64})
#     if length(Pout_pump) != length(steam_in)
#         error("WHAAT")
#     end
#     Ninst = length(Pout_pump)
#     steams = Array{steamObj}(undef,Ninst)
#     wout = Vector{Real}(undef,Ninst)
#     for i = 1:Ninst
#         steams[i],wout[i] = h2o_pump(Pout_pump[i],steam_in[i],etaP[i])
#     end
#     return steams, wout
# end
# h2o_pump(Pout_pump::Vector{<:Real}, steam_in::Vector{steamObj}, etaP::Float64) = h2o_pump(Pout_pump,steam_in,ones(length(Pout_pump)).*etaP)

# """
#     h2o_turbine(Pou::Real, steam_in::steamObj, etaP)
#     Analyses steam turbine
# """
# function h2o_turbine(Pout::Real, steam_in::steamObj, etaT)
#     # Pin is in kPa, Pout is converted to bar
#     h2s = h_ps(Pout, steam_in.s)
#     w   = etaT * (steam_in.h - h2s)
#     h   = steam_in.h - w
#     s2  = s_ph(Pout, h)
#     x2  = x_ps(Pout, s2)
#     T   = T_ps(Pout, s2)
#     steam = steamObj(Pout,T,h,s2,x2)
#     return steam, w
# end

# function h2o_turbine(Pout::Vector{<:Real}, steam_in::Vector{steamObj}, etaT::Vector{Float64})
#     if length(Pout) != length(steam_in)
#         error("WHAAT")
#     end
#     Ninst = length(Pout)
#     steam_out = Array{steamObj}(undef,Ninst)
#     w = Vector{Real}(undef,Ninst)
#     for i = 1:length(Pout)
#         steam_out[i],w[i] = h2o_turbine(Pout[i],steam_in[i],etaT[i])
#     end
#     return steam_out,w
# end
# h2o_turbine(Pout::Vector{<:Real},steam_in::Vector{steamObj},etaT::Float64) = h2o_turbine(Pout,steam_in,ones(size(Pout)).*etaT)

# #need to make this more robust to handle multiple turbine inputs with different number of outputs
# function h2o_simo_turbine(Pout::Vector{<:Real},steam_in::steamObj,etaT::Float64)
#     Ninst = length(Pout)
#     steam_out = Array{steamObj}(undef,Ninst)
#     w = Vector{Real}(undef,Ninst)
#     for i = 1:length(Pout)
#         steam_out[i],w[i] = h2o_turbine(Pout[i],steam_in,etaT)
#     end
#     return steam_out,w
# end

# h2o_simo_turbine(Pout::Vector{<:Real},steam_in::Vector{steamObj},etaT::Float64) = h2o_simo_turbine(Pout,steam_in[1],etaT)

# function steamTable(steamStruct)
#     # Detailed explanation goes here
#     sz = length(steamStruct)
#     tarr = zeros(sz, 5)
#     for i = 1:sz
#         try
#             tarr[i,:] = [steamStruct[i].P, steamStruct[i].T, steamStruct[i].h, steamStruct[i].s, steamStruct[i].x]
#         catch e
#             tarr[i,:] = [-1 -1 -1 -1 -1]
#         end
#     end
#     steamTab = DataFrame(tarr, :auto)
#     rename!(steamTab, [:P_bar, :T_degC, :h_kJ_kg, :s_kJ_kg, :x])
#     return steamTab
# end

# function h2o_massFraction(type::String,steamin::Vector{steamObj},steamout::Vector{steamObj})
#     #H2O_MASSFRACTION Summary of this function goes here
#     #   id = "CFWH" = closed feedwater heater
#     #           Then the first steamin[1],steamout[1] = flow y inlet and outlet
#     #           The second 2 inputs = flow 1-y in, out
#     #   id = "OFWH"
#     #           Then the first input is the outlet
#     #           The second 2 inputs = flow y and 1-y respectively

#         if occursin("ofw",type)
#             h1 = steamout[1].h
#             h2 = steamin[1].h
#             h3 = steamin[2].h
#             y = (h1-h3)/(h2-h3);
#             return y

#         elseif occursin("cfw",type)
#             h1 = steamin[1].h
#             h2 = steamout[1].h
#             h3 = steamin[2].h
#             h4 = steamout[2].h
#             y = (h4-h3)/((h1-h2)+(h4-h3));
#             return y
#         end
# end
