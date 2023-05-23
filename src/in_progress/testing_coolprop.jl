using ModelingToolkit, Revise, Unitful, CoolProp,  Logging,  Printf, XSteam, PyCall
using DifferentialEquations, Printf
Logging.disable_logging(Logging.Warn)
using Plots
using PlotlyJS
include("03-MTK_UTILS.jl")
module Gas
    include("01-ThermoGas.jl")
end
# http://www.coolprop.org/coolprop/HighLevelAPI.html
#=============================================================#
#                   CURRENT METHODS     
#=============================================================#
    @variables t
    Pmax    = 30
    Pmin    = Pmax/3
    Tin     = 800
    η_isen  = 0.9
    begin
        @named hotres   = Gas.SinglePortReservoir(P=30, T=800)
        @named hmsrc    = Gas.GasFlowSource(Ṁ = 1.0)
        @named turbine  = Gas.PassiveThermoTurbine(η = .9)
        @named ep       = Gas.SetPressure(P=30/3)
        connections = vcat(gas_connect(hotres.n,hmsrc.p),
                        gas_connect(hmsrc.n,turbine.p), 
                        gas_connect(turbine.n,ep.p))
        systemNames =[hotres,hmsrc,turbine,ep]
        @named odemodel = ODESystem(connections,t; systems = systemNames)
        smodel = substitute(odemodel,Gas.propDict)
        simple_mod = structural_simplify(smodel)
        prob = ODAEProblem(simple_mod,Pair[],(0.0,1.0))
        sol=solve(prob)
        showsol(systemNames,sol)
    end

#=============================================================#
#                   Using COOL PROPS                            
#=============================================================#
    # isenteropic -> constant entropy

    S_initial   = PropsSI("SMASS", "T", Tin*u"K", "P", Pmax*u"bar","He")
    H_initial   = PropsSI("HMASS", "T", Tin*u"K", "P", Pmax*u"bar","He")

    #   h2 - h1 = integral of cp

    H2s = PropsSI("HMASS", "S", S_initial, "P", Pmin*u"bar","He")
    ΔHs = H2s - H_initial
    ΔHa = ΔHs * η_isen
    H2a = H_initial +ΔHa
    To = PropsSI("T", "Hmass", H2a, "P", Pmin*u"bar","He") 

    T2s = PropsSI("T","SMASS", S_initial, "P", Pmin*u"bar","He").val
    ΔTs = T2s - Tin
    ΔTa = ΔTs * η_isen
    To = Tin + ΔTa
    Ho =  PropsSI("Hmass", "T", To * u"K", "P", Pmin*u"bar","He") 


    @connector function Pin(; name, Pdef = 10.0, Tdef = 300, ṁdef = 0.0)
        ext_var  =  @variables    P(t)=Pdef     T(t)=190    s(t)=0.0    h(t)=191 
        int_var =   @variables    x(t)=0.0      v(t)=.001   cp(t)=0.0
        thru_var =  @variables    ṁ(t)=0.0      Φ(t)=0.0                     # mass flow and energy flow
        ODESystem(Equation[], t, [thru_var..., ext_var..., int_var...], ps; name = name)
    end


#=============================================================#
#                   Graphical Approach                       
#=============================================================#
# State 1 (Pre turbine) with compression ration = 3.5
P_i = 80.0u"bar"
T_i = 750u"°C"
s_i   = PropsSI("SMASS", "T", T_i, "P", P_i,"He");
v_i   = PropsSI("D", "T", T_i, "P", P_i,"He")^(-1);
u_i   = PropsSI("UMASS", "T", T_i, "P", P_i,"He");
h_i   = PropsSI("HMASS", "T", T_i, "P", P_i,"He");



@printf "State 1: P = %.2f %s,\tT = %.2f %s,\ts = %.2e %s, u = %.2e %s,\tv = %.2f %s \n" P_i.val unit(P_i) T_i.val unit(T_i) s_i.val unit(s_i) u_i.val unit(u_i) v_i.val unit(v_i)
#compression ration = 3.5
r_p = 3.5

#outlet pressure
P_o = P_i/r_p

# P array = P_i -> P_o
P_arr = LinRange(P_i,P_o,100)

# Isentropic Process
vfcn = PropsSI.("D","SMASS",s_i,"P", P_arr,"HE")

using PyCall, PyPlot, CoolProp
begin
    pyimport("sys").executable
    pyimport_conda("CoolProp.CoolProp", "CoolProp")
    pyimport_conda("CoolProp.Plots", "CoolProp")
    @pyimport CoolProp.CoolProp as CP
    @pyimport CoolProp.Plots as CPlot
end
p = CPlot.PropertyPlot("Helium", "TS")
p.calc_isolines()
p.show();
fig1 = gcf()
display(fig1)
# CCP = CoolProp.CoolProp

# CP.PropertyPlot("He","ph")

