using ModelingToolkit, Plots, DifferentialEquations, Revise
@variables t
# C:\Users\harvey\.julia\dev\test\ThermalSystem_Models.jl\src\in_progress\01-FoundationBlocks.jl
#
# Parameters T, P, mflow, θ where mflow * theta = E

@connector function FlowPin(; name)
    sts = @variables P=1.0 ṁ=1.0 [connect = Flow]
    # ODESystem(Equation[], t, sts, []; name = name)
    NonlinearSystem(Equation[],  sts, []; name = name)
end

@connector function HeatPin(; name)
    sts = @variables T=300 Q̇=0.0 [connect = Flow]
    # ODESystem(Equation[], t,sts, [], name = name)
    # NonlinearSystem(Equation[], t, sts, []; name = name)
    NonlinearSystem(Equation[],  sts, []; name = name)
end

# convention = +Q = Qin
@connector function HeatFlowPin(; name)
    sts = @variables Q̇=0.0 [connect = Flow]
    # ODESystem(Equation[], t,sts, [], name = name)
    NonlinearSystem(Equation[],  sts, []; name = name)
end

# convention = +W = Win
@connector function WorkPin(; name)
    sts = @variables Ẇ = 0.0 [connect = Flow]
    # ODESystem(Equation[], t, sts, []; name = name)
    NonlinearSystem(Equation[],  sts, []; name = name)
end

@connector function ThermoPin(; name, cp = 5192.0)
    sts = @variables P=1.0 ṁ=1.0 [connect = Flow] T = 300 θ = 0.0 [connect = Flow]
    ps = @parameters cp = 5192.0 Tref = 293.15
    sts = [P,T,ṁ,θ]
    eqs = [0 ~ θ - ṁ * cp * (T-Tref)]
    # ODESystem(Equation[], t,sts, [], name = name)
    NonlinearSystem(eqs,  sts, ps; name = name)
end

@component function ThermoGround(; name)
    @named gnd = ThermoPin()
    eqs = [gnd.P ~ 0
            gnd.T ~ 0]
    compose(NonlinearSystem(eqs,  [], []; name = name),gnd)
end

@component function ThermoPort(; name)
    @named inlet = ThermoPin()
    @named outlet = ThermoPin()
    @named q = HeatFlowPin()
    @named w = WorkPin()

    sts = @variables P=0.0 T=0.0 ṁ=1.0

    eqs = [ T ~ outlet.T - inlet.T
            P ~ outlet.P - inlet.P
            0 ~ inlet.ṁ - outlet.ṁ
            0 ~ inlet.ṁ - ṁ
            0 ~ q.Q̇ + w.Ẇ + inlet.ṁ * inlet.θ - outlet.ṁ * outlet.θ]

    # compose(ODESystem(eqs, t, sts, []; name = name), inlet, outlet, q, w)
    compose(NonlinearSystem(eqs,  sts, []; name = name),inlet,outlet,q,w)
end

@component function PassiveThermoPort(; name)
    @named inlet = ThermoPin()
    @named outlet = ThermoGround()

    sts = @variables P=0.0 T=0.0 ṁ=1.0

    eqs = [
        0 ~ inlet.T - T
        0 ~ inlet.P - P
        0 ~ inlet.ṁ - outlet.gnd.ṁ
        0 ~ inlet.ṁ - ṁ
    ]
    # extend(ODESystem(eqs, t, [], []; name = name), tp)
    compose(NonlinearSystem(eqs, sts, []; name = name),inlet,outlet)
end

# compression ratio is always A --> B
@component function gas_compressor(;name)
    @named tp = ThermoPort()
    ps = @parameters rp η cp k
    kcoeff = (k-1)/k
    eqs = [
        0 ~ rp* tp.inlet.P - tp.P
        0 ~ tp.inlet.T*(1 - rp^kcoeff)/(-η) - tp.T
        0 ~ tp.w.Ẇ - tp.ṁ * cp * tp.T
        0 ~ tp.q.Q̇
        0 ~ tp.ṁ * tp.outlet.θ - (tp.inlet.θ * tp.ṁ + tp.w.Ẇ)
    ]
    # extend(ODESystem(eqs, t,[], ps; name = name), tp)
    extend(NonlinearSystem(eqs,  [], ps; name = name),tp)
end

# compression ratio is always A --> B
@component function gas_turbine(;name)
    @named tp = ThermoPort()
    ps = @parameters rp η cp k
    kcoeff = (k-1)/k
    eqs = [
        0 ~ rp* tp.inlet.P - tp.P
        0 ~ η*(rp^kcoeff-1)*tp.inlet.T - tp.T 
        0 ~ tp.w.Ẇ - tp.ṁ * cp * T
        0 ~ tp.q.Q̇
        0 ~ tp.ṁ * tp.outlet.θ - (tp.inlet.theta * tp.ṁ + tp.w.Ẇ)
    ]
    # extend(ODESystem(eqs, t,[], ps; name = name), tp)
    extend(NonlinearSystem(eqs,  [], ps; name = name),tp)
end

@component function MassFlowSource(; name,  Ṁ = 1.0 , Pin=10 , Tin = 300 )
    @named tp = ThermoPort()
    ps = @parameters Ṁ = Ṁ Pin = Pin Tin = Tin
    eqs = [
            0 ~ tp.ṁ - Ṁ
            0 ~ tp.P - Pin
            0 ~ tp.T - Tin
            0 ~ tp.q.Q̇
            0 ~ tp.w.Ẇ]

    # extend(ODESystem(eqs, t, [], ps; name = name), tp)
    extend(NonlinearSystem(eqs,  [], ps; name = name),tp)
end


function testcase()
    @named g1 = ThermoGround()
    @named mf = MassFlowSource()
    # @named cs = gas_compressor()
    # @named ptm = PassiveThermoPort()
    println("Here")

    # con = [connect(g1.gnd,mf.inlet)]
    #     connect(mf.outlet,cs.inlet)
    #     connect(cs.outlet,ptm.inlet)]



    @named odesys = NonlinearSystem(con,[],[]; systems = [g1,mf,cs,ptm])

    @named nlsys = NonlinearSystem(equations(expand_connections(odesys)), states(odesys),parameters(odesys))
end
