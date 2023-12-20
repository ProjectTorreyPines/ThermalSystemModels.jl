module FastBrayton
    using LinearAlgebra

    struct BraytonOutput{T<:Real}
        η_th::T
        w_net::T
        w_in::T
        w_out::T
        r_bw::T
        q_H::T
        q_L::T
        T_HX::T
    end


    """

    braytonCycle(rp::Union{Float64,Int}, 
                        Pmax::Union{Float64,Int}, 
                        Tmin::Union{Float64,Int}, 
                        Tmax::Union{Float64,Int}, 
                        Nt::Int, 
                        Nc::Int, 
                        ηt::T=0.93, 
                        ηc::T=0.89, 
                        ϵr::T=0.0;
                        cp::Union{Float64,Int} = 5.1926e3,
                        cv::Union{Float64,Int} = 3.1156e3) where {T<:Real}
                        Brayton evaluates the thermal performance for a specified brayton cycle
        
        ===================================================================
        INPUTS        UNIT        DESCRIPTION
        rp          Scalar      Total Compression Ratio
        Pmax        [Pa]        Specified pressure, maximum
        Tmin        [K]         Lowest Cycle Temp
        Tmax        [K]         Highest Cycle Temp
        Nc,Nt       Scalars     Number of compression (Nc) and Turbine (Nt) stages
        ηc          Frac        Fractional compressor isentropic effeciency
        ηt          Frac        Fractional Turbine isentropic effeciency
        ϵr          Frac        Fractional Regenerator effectiveness

    """
    function braytonCycle(rp::Union{Float64,Int}, 
                            Pmax::Union{Float64,Int}, 
                            Tmin::Union{Float64,Int}, 
                            Tmax::Union{Float64,Int}, 
                            Nt::Int, 
                            Nc::Int; 
                            ηt::T=0.93, 
                            ηc::T=0.89, 
                            ϵr::T=0.0,
                            cp::Union{Float64,Int} = 5.1926e3,
                            cv::Union{Float64,Int} = 3.1156e3) where {T<:Real}
        # cp = 5.1926e3
        # cv = 3.1156e3

        #SUMMARIZED COEFFECICENTS
        a1t = (ηt * (rp^(1.0 / Nt))^(cv / cp) - ηt * rp^(1.0 / Nt) + rp^(1.0 / Nt)) / rp^(1.0 / Nt)
        a1c = (ηc * (rp^(1.0 / Nc))^(cv / cp) - (rp^(1.0 / Nc))^(cv / cp) + rp^(1.0 / Nc)) / (ηc * (rp^(1.0 / Nc))^(cv / cp))
        Pmin = Pmax / rp

        #STATE VARIABLES
        phi = [Tmin; Tmax]
        θ_L = [Tmin; Pmin]
        θ_H = [Tmax; Pmax]

        #MATRICES, DEFINED IN WRITE UP
        A_t = Diagonal([a1t, rp^(-1.0 / Nt)])
        B_t = Diagonal([1 / a1t, 1.0])           #turbine
        A_c = Diagonal([a1c, rp^(1.0 / Nc)])
        B_c = Diagonal([1 / a1c, 1.0])#Comp
        C = Diagonal([cp * (a1c - 1.0), cp * (a1t - 1.0)])
        N = [Nc 0; (1-Nc) 0; 0 Nt; 0 (1-Nt)]
        E = [(1.0-ϵr) ϵr; ϵr (1.0-ϵr)]

        #OUTPUTS OF THE COMPRESSOR AND TURBINE CIRCUITS
        θ_ci = (A_c * B_c)^(Nc - 1) * A_c * θ_L    #AFTER COMP
        θ_ti = (A_t * B_t)^(Nt - 1) * A_t * θ_H    #AFTER TURB

        θ_1 = E * [θ_ci[1]; θ_ti[1]]  #OUTLET TEMPERATURES OF REGENERATOR
        θ_2 = [Tmax; Tmin]            #FOR USE IN THE FOLLOWING

        #RESULTS
        wc, q_intercool, wt, q_reheat = N * C * phi           #[COMP WORK, INTERCOOL HEAT LOSS, TURBINE WORK, REHEAT HEAD ADDITION]
        q_primary, q_waste = cp * (θ_2 - θ_1)                 #FUSION HEAT ADDITION, HEAT WASTED
        wnet = abs(wt) - abs(wc)                              #NET WORK OUT 
        q_in = abs(q_primary) + abs(q_reheat)                 #NET HEAT ADDED TO THE SYSTEM
        q_out = abs(q_waste) + abs(q_intercool)               #NET HEAT REMOVED FROM SYSTEM
        η_thermal = wnet / q_in                               #THERMAL EFFECIENCY
        Primary_inlet_temp = θ_1[1]                           #INLET HELIUM TEMPERATURE AT REACTOR HEAT EXCHANGERS
        return BraytonOutput(η_thermal, abs(wnet), abs(wc), abs(wt), abs(wc / wt), abs(q_in), abs(q_out), Primary_inlet_temp)
    end
end