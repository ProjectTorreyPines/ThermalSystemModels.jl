module fluid_properties

    using Printf, XSteam

    function fprops(name::Symbol; kw...)
            if name == :helium
                return helium_properties(kw...)
            elseif name == :air
                return air_properties(kw...)
            elseif name == :pbli
                return pbli_properties(kw...)
            elseif name == :water
                return water_properties(kw...)
            else
                error("Unsupported fluid name")
            end
    end

    # object to hold fluid properties
    mutable struct FluidProp{Ty <: Union{Float64, Int64}}
        T::Ty       # Temperautre K
        P::Ty       # Pressure Bar
        cp::Ty      # constant pressure specific heat [J/kg/K]
        cv::Ty      # constant volume specific heat [J/kg/K]
        k::Ty       # specific heat ratio [-]
        v::Ty       # specific volume [m³/kg]
        R::Ty       # gas constant [J/kg/K]
        s::Ty       # specific entropy [J/kg/K]
        h::Ty       # specific enthalpy [J/kg/K]
        x::Ty       # quality [-]
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
    function pbli_properties(; T::Union{Float64, Int64} = 508.0, P::Union{Float64, Int64} = 10.0, verbose = false)
        if T > 1615.15
            println("Exceeded lead lithium boiling point $(Tin)")
        elseif T < 453 && verbose
            println("Warning - lead lithium temperature below melting point $(Tin)")
        end
        # Pin = bar
        # data from [2]
        cpin = (0.195 - 9.116e-6 .* T) .* 1000  # J/kgK
        vin = 1 / (10520.35 - 1.19051 .* T)      # m³/kg

        # specific entropy of lead lithium, with reference to 508 K
        a = 0.195*1e3
        b = 9.11e-6*1000
        s_pbli      = a*log(T) - b*T - 1210.316002277804
        u(T)        = 195 * T - .0009116/2 * T^2                # integral of cp
        #                                                   h ~ cpT + vP
        return FluidProp(T, P, cpin, cpin, 1.0, vin, 1.0, s_pbli, cpin * T + vin * (P * 1e5) , 0.0)
    end

    function superheatedwater(P::Real,T::Real)
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

    export fprops, FluidProp, gastype
end