using ModelingToolkit, Plots, Revise, Unitful
using Symbolics, Logging
using XSteam
using Printf
using Graphs, Plots, GraphRecipes
@variables t
Logging.disable_logging(Logging.Warn)
MTK = ModelingToolkit

##
#       Need todo: Coolprops, validation
#               Connect work nodes - master energy balance?
##

#   PROPERTIES DICT
begin
    @register_symbolic cpfunc(x,y)
    @register_symbolic kfunc(x,y)
    @register_symbolic vphfunc(x,y)

    @register_symbolic Tpsfunc(x,y)
    @register_symbolic Tphfunc(x,y)

    @register_symbolic xphfunc(x,y)
    @register_symbolic xptfunc(x,y)

    @register_symbolic sptfunc(x,y)
    @register_symbolic sphfunc(x,y)

    @register_symbolic hptfunc(x,y)
    @register_symbolic hpsfunc(x,y)
    @register_symbolic hsatfunc(x)

    hydro_prop_dict = Dict(Tpsfunc => T_ps,
                            Tphfunc => T_ph,
                            xphfunc => x_ph,
                            sptfunc => s_pT,
                            sphfunc => s_ph,
                            hptfunc => h_pT,
                            hpsfunc => h_ps,
                            hsatfunc => hL_p,
                            vphfunc => v_ph);

end

#   STEAM PINS
"""
    BasicSteamPin()
    Self computes T,s,x,V
    Must have methods for ṁ,Φ,P,h
"""
@connector function BasicSteamPin(; name, Pdef=0.1)
    across_var  = @variables  P(t)=Pdef T(t)=190 s(t)=0.0 h(t)=191 x(t)=0.0 v(t)=.001
    thru_var    = @variables  ṁ(t)=0.0 Φ(t)=0.0                     # mass flow and energy flow
    
    ps = []
    
    eqs = Equation[T ~ Tphfunc(P,h)
                    s ~ sphfunc(P,h)
                    x ~ xphfunc(P,h)
                    v ~ vphfunc(P,h)]

    ODESystem(eqs, t, [thru_var..., across_var...], ps; name = name, defaults = [P => Pdef, ṁ => 0, Φ => 0, h=>191])
end

@connector function WorkPin(; name)
    sts = @variables Ẇ(t) = 0.0 
    ODESystem(Equation[], t, sts, []; name = name)
end

@connector function HeatTransferPin(; name)  #input in W
    sts = @variables Q̇(t) = 0.0
    ODESystem(Equation[], t, sts, []; name = name, defaults = [Q̇ => 0.0])
end

@connector function FixedHeatFlowPin(; name, Qin = 1e3)  #input in W
    sts = @variables Q̇(t)=Qin 
    ps = @parameters Qin = Qin
    eqs=[Q̇ ~ Qin]
    ODESystem(eqs,t, sts, ps; name = name, defaults = [Q̇ => Qin])
end
# @named bsp_test = BasicSteamPin()
# p = WorkPin(name = :testwork)
# q = FixedHeatFlowPin(name = :testfixedheatflow)
# ht = HeatTransferPin(name = :testHeatTransferPin)

"""
    hydro_connect(pins...)
    Adds mass and energy balance eqs to all nodes

"""
function hydro_connect(pins...)

    eqs = [
        sum(pin -> pin.ṁ, pins) ~ 0.0, # Mass
        sum(pin -> pin.Φ, pins) ~ 0.0, # Energy
    ]

    for pin in pins[2:end]
        push!(eqs, pins[1].P ~ pin.P)   # p_1 = p_2, p_1 = p_3,...
        push!(eqs, pins[1].h ~ pin.h)   # p_1 = p_2, p_1 = p_3,...
    end
    return eqs
end

function hydro_basic_connect(n,p)
    eqs = [
        0 ~ p.ṁ + n.ṁ
        0 ~ p.Φ + n.Φ
        p.h ~ n.h
        p.P ~ n.P
    ]
    return eqs
end

function hydro_series_connect(comp ;  returnmode = :eq)
    eqs =  hydro_basic_connect(comp[1].n,comp[2].p)   # n -> p

    if length(comp)>2
        for i = 2:length(comp)-1
            eqs = vcat(eqs,  hydro_basic_connect(comp[i].n,comp[i+1].p))
        end
    end

    if returnmode == :ode
        return ODESystem(eqs,t; name = :connections)
    end
    
    return eqs
end

function extenda(odevec::Vector{ODESystem})
    if length(odevec)==1
        return odevec[1]
    else
        return extend(odevec[1],extenda(odevec[2:end]))
    end
end

begin
    @component function HeatTransferPort(; name)
        @named p = HeatTransferPin()
        @named n = HeatTransferPin()
        eqs = [0 ~p.Q̇ + n.Q̇]
        @variables Q̇(t) 
        ODESystem(Equation[],t, [Q̇], []; name = name, systems = [p,n])
    end

    @component function Reservoir(;name, P = 0.1)
        @named n = BasicSteamPin(Pdef = P)
        ps = @parameters P=P 
        eqs = [
            n.P ~ P
            n.h  ~ hsatfunc(P)
            n.Φ ~ 0
        ]
        ODESystem(eqs, t,[], ps; name = name, systems = [n], defaults = [P => 0.1])
    end

    @component function SuperHeatedReservoir(;name, P = 150, T = 600)
        @named n = BasicSteamPin(Pdef = P)
        ps = @parameters P=P T=T 
        eqs = [
            n.P ~ P
            0 ~ hptfunc(P,T) + n.h
            n.Φ ~ -n.ṁ * n.h
        ]
        ODESystem(eqs, t,[], ps; name = name, systems = [n], defaults = [P => 0.1, T => T])
    end

    @component function ioReservoir(;name, P = 0.1)
        @named p = BasicSteamPin()
        @named n = BasicSteamPin(Pdef = P)
        ps = @parameters P=P
        eqs = [
            n.P ~ P
            n.h  ~ hsatfunc(P)
            n.Φ ~ 0
            0~n.ṁ + p.ṁ
        ]
        ODESystem(eqs, t,[], ps; name = name, systems = [n,p], defaults = [P => 0.1 ])
    end

    @component function SetPressure(;name, P = 0.1)
        @named p = BasicSteamPin(Pdef = P)
        @named n = BasicSteamPin(Pdef = P)
        ps = @parameters P(t)=P 
        eqs = [
            p.P ~ P
            n.P ~ p.P
            n.h ~ p.h
            0 ~ p.Φ + n.Φ             # conservation of energy
            0 ~ p.ṁ + n.ṁ
        ]
        ODESystem(eqs, t,[], ps; name = name, systems = [p,n], defaults = [P => 0.1])
    end
end

begin 
    @component function SteamFlowSource(;name, ṁ = 1.0)
        @named p = BasicSteamPin()
        @named n = BasicSteamPin()
        ps = @parameters Ṁ = ṁ
        eqs = [
                0 ~ p.Φ + n.Φ   # conservation of energy
                0 ~ p.ṁ + n.ṁ
                p.ṁ ~ Ṁ
                p.h ~ n.h 
                n.P ~ p.P
            ]
        ODESystem(eqs, t, [], ps; name = name, systems = [n,p], defaults = [Ṁ => ṁ])
        # compose(ODESystem(eqs, t, [], ps; name = name, defaults = [Ṁ => 1.0]),p,n)
    end

    @component function AdiabaticPump(;name, η = 0.6, setpressure = false ,Pout = 10)
        #Pout in bar
        @named p = BasicSteamPin()
        @named n = BasicSteamPin()
        @named w = WorkPin()
        
        ps = @parameters η = η P=Pout

        eqs = [
            w.Ẇ ~ p.Φ + n.Φ                         # conservation of energy
            0 ~ p.ṁ + n.ṁ
            n.h  ~ p.h + p.v *100* (n.P - p.P) / η  # work, multiply by 100 to get to kPa
            w.Ẇ  ~ p.ṁ * (n.h - p.h)
        ]

        if setpressure
            eqs = vcat(eqs,n.P  ~ P)
        end
        # extenda([ODESystem(eqs, t,[], ps; name = name, defaults = [η => 0.6, P => Pout]), p,n,w])
        ODESystem(eqs, t,[], ps; name = name, systems = [p,n,w], defaults = [η => .6, P => 10])
    end

    @component function Splitter(;name)
        @named p = BasicSteamPin()
        @named y = BasicSteamPin()
        @named z = BasicSteamPin()

        #  0 ~ pMax₊P - trbsimo₊split₊z₊P(t)
        eqs = [
            p.h ~ y.h               # enthalpy
            0 ~ y.h - z.h
            p.P ~ y.P               # pressure
            z.P ~ y.P 
            0 ~ p.ṁ + y.ṁ + z.ṁ
            0 ~ p.Φ + y.Φ + z.Φ
        ]

        ODESystem(eqs, t,[], []; name = name, systems = [p,y,z])
    end

    @component function AdiabaticTurbine(;name, η = 1.0, setpressure = false, Pout = 0.1)
        #Pout in bar
        @named p = BasicSteamPin()
        @named n = BasicSteamPin()
        @named w = WorkPin()
        ps = @parameters P=Pout η = η

        eqs = Equation[
            0 ~ p.ṁ + n.ṁ      
            n.h  ~ p.h - (p.h-hpsfunc(n.P,p.s))*η
            w.Ẇ  ~ p.ṁ * (n.h - p.h)
            w.Ẇ ~ p.Φ + n.Φ                         # conservation of energy
        ]

        if setpressure
            println("we be setting pressure in this bitch")
            eqs = vcat(eqs, n.P ~ P)
        end

        # extenda([ODESystem(eqs, t,[], ps; name = name, defaults = [η => 0.6, P => Pout]), p,n,w])
        ODESystem(eqs, t,[], ps; name = name, systems = [p,n,w], defaults = [η => 1.0, P => 10])
    end   

    #  0 ~ pMax₊P - trbsimo₊split₊z₊P(t)
    @component function SIMO_AdiabaticTurbine(;name, η = 1.0, setpressure = false, Pout = [12,0.1])
        #Pout in bar
        ps = @parameters Py=Pout[1] Pz=Pout[2] η = η

        sp = setpressure

        @named p = BasicSteamPin()
        @named split = Splitter()
        @named hp = AdiabaticTurbine(η = η, Pout = Py, setpressure = sp)
        @named lp = AdiabaticTurbine(η = η, Pout = Pz, setpressure = sp)
        
        sts = @variables Whp(t)=0.0 Wlp(t)=0.0

        eqs = Equation[
            Whp ~ hp.w.Ẇ
            Wlp ~ lp.w.Ẇ
            ]

        # connections

        incon = hydro_basic_connect(p,split.p)
        ycon = hydro_basic_connect(split.y,hp.p)
        zcon = hydro_basic_connect(split.z,lp.p)

        eqs = vcat(eqs,incon,ycon,zcon)

        ODESystem(eqs, t,sts, ps; name = name, systems = [p,split,hp,lp], defaults = [η => 1.0 Py => 10 Pz => 0.1])
        # extend(ODESystem(eqs, t,sts, ps; name = name, systems = [hp,lp], defaults = [η => 1.0 Py => 10 Pz => 0.1]),split)
    end

    @component function SIMOAdiabaticTurbine(;name, ηin = 1.0, setpressure = false, Pyin = 10, Pzin  = 0.1)
        #Pout in bar
        ps = @parameters Py=Pyin Pz=Pzin η = ηin
        sp = setpressure

        # EXTERNAL NODES FOR INTERFACING
        @named p = BasicSteamPin()  # inlet node
        @named w = WorkPin()
        @named hp = AdiabaticTurbine(η = ηin, Pout = Pyin, setpressure = sp)
        @named lp = AdiabaticTurbine(η = ηin, Pout = Pzin, setpressure = sp)

        #
        #                yn(-) -------> y
        # IN --> p (+) --|
        #                zn(-) -------> z
        #
        # INTERNAL NODES
        @named yn = BasicSteamPin()
        @named zn = BasicSteamPin()

        split_connect = hydro_connect(p, yn, zn)    # mflow, p.ṁ => positive
        hp_connect = hydro_connect(yn,hp.p)
        lp_connect = hydro_connect(zn,lp.p)


        eqs = vcat(split_connect,hp_connect,lp_connect)


        # if setpressure == true
        #     display("SIMO TURBINE PRESSURE MODE => SET")
        #     eqs = vcat(eqs,[y.P ~ Py; z.P ~ Pz])
        # end

        compose(ODESystem(eqs, t,[], ps; name = name, systems = [yn,zn,p], defaults = [η => 1.0, Py => 10, Pz => 0.1]), hp,lp)
        # ODESystem(eqs, t,[], ps; name = name, systems = [p,hp,lp], defaults = [η => 1.0, Py => 10, Pz => 0.1])
        # extend(ODESystem(eqs, t,sts, ps; name = name, systems = [hp,lp], defaults = [η => 1.0 Py => 10 Pz => 0.1]),split)
    end

    @component function IdealBoiler(;name, Tout = 350)
        #Pout in bar
        @named p = BasicSteamPin()
        @named n = BasicSteamPin()
        @named q = HeatTransferPin()
        
        ps = @parameters T = Tout

        eqs = Equation[
            q.Q̇ ~ p.Φ + n.Φ             # conservation of energy
            0 ~ p.ṁ + n.ṁ
            n.P  ~  p.P                 # no pressure
            n.h  ~ hptfunc(p.P,T)       # work, multiply by 100 to get to kPa
            q.Q̇ ~ p.ṁ * (n.h - p.h)
        ]
        ODESystem(eqs, t,[], ps; name = name, systems = [p,n,q], defaults = [ T => Tout])
    end

    @component function IdealCondensor(;name)
        #Pout in bar
        @named p = BasicSteamPin()
        @named n = BasicSteamPin()
        @named q = HeatTransferPin()
        
        eqs = Equation[
            q.Q̇ ~ p.Φ + n.Φ             # conservation of energy
            0 ~ p.ṁ + n.ṁ
            n.P  ~  p.P                 # no pressure
            n.h  ~ hsatfunc(p.P)        # work, multiply by 100 to get to kPa
            q.Q̇ ~ p.ṁ * (n.h - p.h)
        ]
        ODESystem(eqs, t,[], []; name = name, systems = [p,n,q])
    end

    @component function OpenFeedwaterHeater(;name)
        # flows x and y are the inlets
        @named y = BasicSteamPin()
        @named z = BasicSteamPin()
        @named n = BasicSteamPin()  

        sts=@variables yfrac(t) = 0.5
        # n.Φ ~ yfrac * y.Φ + (1-yfrac) * z.Φ
        # n.P ~ y.P
        eqs  =[
            n.P ~ z.P
            n.h ~ hsatfunc(n.P) # sat liquid
            0 ~ n.h - (yfrac * y.h + (1-yfrac) * z.h)
            0 ~ n.ṁ*yfrac + y.ṁ
            0 ~ n.ṁ + y.ṁ + z.ṁ
            0 ~ n.Φ + y.Φ + z.Φ
        ]
        ODESystem(eqs, t,sts, []; name = name, systems = [y,z,n], defaults = [yfrac => 0.5])
    end
end

begin
    @named res          = Reservoir(P = 150)
    @named lowres       = Reservoir(P=0.1)
    @named iores        = ioReservoir(P=0.1)
    @named ioresB       = ioReservoir(P=10)
    @named freenodeA    = BasicSteamPin(Pdef = 150)

    @named valve        = SteamFlowSource(ṁ = 1)
    @named pumpB         = AdiabaticPump(Pout = 150, setpressure = true)
    @named boil         = IdealBoiler(Tout = 600)
    @named turbine      = SIMOAdiabaticTurbine(setpressure = true, Pyin = 10, Pzin = 0.1,ηin = 1.0) 
    @named pumpA         = AdiabaticPump(Pout = 10, setpressure = true)
    @named cond         = IdealCondensor()
    @named openfw       = OpenFeedwaterHeater()

    connections = vcat(hydro_connect(pumpB.n,valve.p),
                        hydro_connect(valve.n, boil.p),
                        hydro_connect(boil.n, turbine.p),
                        hydro_connect(turbine.hp.n, openfw.y),
                        hydro_connect(turbine.lp.n, cond.p),
                        hydro_connect(cond.n,iores.p),
                        hydro_connect(iores.n,pumpA.p),
                        hydro_connect(pumpA.n,openfw.z),
                        hydro_connect(openfw.n,pumpB.p))

    systemNames =[valve,
                    boil,
                    turbine,
                    pumpB,
                    res,
                    iores,
                    pumpA,
                    cond,
                    openfw]
                    
    @named odemodel = ODESystem(connections,t; systems = systemNames)

    smodel = substitute(odemodel,hydro_prop_dict)
    simple_mod = structural_simplify(smodel)
    prob = ODAEProblem(simple_mod,Pair[],(0.0,1.0))

    plotmod(smodel)
    prob.f.observed(openfw.n.P,prob.u0,prob.p,1.0)

    pf(xin) = prob.f.observed(xin,prob.u0,prob.p,1.0)
    nodes = [turbine.p, turbine.yn, turbine.zn, turbine.hp.p, turbine.lp.p,turbine.hp.n, turbine.lp.n,  cond.p, cond.n, pumpA.n,pumpB.p, openfw.z, openfw.y]
    for n in nodes
        @printf "\t %s ṁ = %.2f \t P = % .2f \t h = %.2f\n" n.name pf(n.ṁ) pf.(n.P) pf(n.h)
    end
end



for s in odemodel.systems
@show propertynames(s)
end














#   WORKING
begin
    @named res          = Reservoir(P = 10)
    @named iores        = ioReservoir(P=10)
    @named freenodeA    = BasicSteamPin(Pdef = 150)

    @named valve        = SteamFlowSource(ṁ = 1)
    @named pumpB         = AdiabaticPump(Pout = 150, setpressure = true)
    @named boil         = IdealBoiler(Tout = 600)
    @named turbine      = SIMOAdiabaticTurbine(setpressure = true, Pyin = 10, Pzin = 0.1,ηin = 1.0) 
    @named pump         = AdiabaticPump(Pout = 10, setpressure = true)
    @named cond         = IdealCondensor()
    @named openfw       = OpenFeedwaterHeater()

    connections = vcat(hydro_connect(res.n, valve.p),
                        hydro_connect(valve.n,boil.p),
                        hydro_connect(boil.n, turbine.p),
                        hydro_connect(turbine.hp.n, openfw.y),
                        hydro_connect(turbine.lp.n, cond.p),
                        hydro_connect(cond.n,pump.p),
                        hydro_connect(pump.n,openfw.z),
                        hydro_connect(openfw.n,freenodeB),
                        hydro_connect(freenodeB,pumpB.p),
                        hydro_connect(pumpB.n,iores.p))

    systemNames =[res,valve,boil,turbine,pump,cond,openfw,freenodeB]

    @named odemodel = ODESystem(connections,t; systems = systemNames)

    smodel = substitute(odemodel,hydro_prop_dict)
    simple_mod = structural_simplify(smodel)
    prob = ODAEProblem(simple_mod,Pair[],(0.0,1.0))
end

begin
    println("Turbine")
    nodes = [turbine.p, turbine.yn, turbine.zn, turbine.hp.p, turbine.lp.p, cond.p, cond.n, openfw.z, openfw.y]
    for n in nodes
        @printf "\t %s ṁ = %.2f \t P = % .2f\n" n.name pf(n.ṁ) pf.(n.P)
    end
end

function plotmod(smodel)
    # structural_simplify(model)

    var2idx = Dict(s => i for (i, s) in enumerate(states(smodel)))
    idx2var = Dict(i => s for (i, s) in enumerate(states(smodel)))

    varvar = varvar_dependencies(asgraph(smodel), variable_dependencies(smodel))
        
    eqgraph = asgraph(equation_dependencies(smodel),var2idx)
    equation_dep_digraph = asdigraph(eqgraph,smodel)
    equation_dep_digraph.fadjlist
    varvar = varvar_dependencies(asgraph(smodel), variable_dependencies(smodel))
    varvar.badjlist
    # # asd = asdigraph(digr,odesys)
    nm = [idx2var[idx] for idx in 1:nv(varvar)]
    p=plot(varvar,
            size = (1500,1500),
            names = nm,
            method = :sfpd)
    

    return p,var2idx,idx2var,varvar,eqgraph
        
end

p,var2idx,idx2var,varvar,eqgraph=plotmod(smodel)
    

begin
    println("pump inlet")
    node = pump.n
    @show pf(node.h)
    @show pf(node.P)
    @show pf(node.s)
    @show pf(node.v)
    @show pf(node.x)
    @show pf(node.Φ)
    @show pf(node.ṁ)
    @show pf(node.T)
end

# begin
#     vlv         = SteamFlowSource(name = :vlv, ṁ = 1)
#     pmpa        = AdiabaticPump(name = :pmpa, Pout = 12,setpressure = true)
#     pmpb        = AdiabaticPump(name = :pmpb, Pout = 150,setpressure = true)
#     boil        = IdealBoiler(name = :boil, Tout = 600)
#     trbsimo     = SIMO_AdiabaticTurbine(name = :trbsimo, setpressure = false) 
#     openfw      = OpenFeedwaterHeater(name = :openfw)
#     cond        = IdealCondensor(name = :cond)
#     ioRes       = ioReservoir(name = :ioRes, P = 0.1)
# end


# begin
#     con2 = hydro_basic_connect(pmpa.n,openfw.z)
#     con3 = hydro_series_connect([openfw,vlv,pmpb,boil,trbsimo])
#     con4 = hydro_basic_connect(trbsimo.hp.n,openfw.y)
#     con5 = hydro_series_connect([trbsimo.lp,cond,ioRes,pmpa])
#     connections = vcat(con2,con3,con4,con5)
# end


# begin
# odesys = ODESystem(connections,t; name = :odesys, systems = [vlv,pmpa,pmpb,boil,trbsimo,openfw,cond,ioRes])
# model = substitute(odesys,hydro_prop_dict)
# end




# @named sftest = SteamFlowSource()
# # sftest.defaults
# @named pftest = AdiabaticPump()
# # pftest.defaults
# # compose(sftest,pftest)
# @named sisotest = SIMO_AdiabaticTurbine()
# mode = ODEProblem(model,model.defaults,parameters(model))
# parameters(model)
# parameters(model)
# parameters(model)
# probfun = ODEFunction(model)
# MTK.moderntoolkitize(model)
# model.defaults
# model.states
# states(model)
# equations(model)
# dump(MTK.get_u0_p(model, model.defaults, parameters(model)))
# s[P]
# states(model)