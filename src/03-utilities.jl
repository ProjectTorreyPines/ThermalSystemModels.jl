
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

function hxeff(Thi::Real, Ch::Real, Tci::Real, Cc::Real, eff::Real; verbose = false)
    # if Ch and Cc are just cp, Qact = Qspecific
    Cmin = min(Ch, Cc)
    dt_mx = Thi - Tci
    Q_act = Cmin .* dt_mx .* eff
    dtH = Q_act ./ Ch
    dtC = Q_act ./ Cc
    Tho = Thi - dtH
    Tco = Tci + dtC
    if verbose == true
        @printf "HX eff = %.2f, Tci = %.2f, Thi = %.2f \n" eff Tci Thi
        @printf "\t ΔT max = %.2f \n\t ΔT cold = %.2f => Tco = %.2f \n\t ΔT hot = %.2f  => Tho = %.2f\n" dt_mx dtC Tco dtH Tho
    end
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
    ΔT_1 = abs(Thi - Tco)
    ΔT_2 = abs(Tho - Tci)

    ΔT_logmean = (ΔT_1-ΔT_2)/log(ΔT_1/ΔT_2);
    return ΔT_1,ΔT_2,ΔT_logmean
end

repmat(v,n)  =  repeat(v;inner=n)
repelem(v,n) =  repeat(v;outer=n)

bar2pascal(P_bar) = P_bar * 100000