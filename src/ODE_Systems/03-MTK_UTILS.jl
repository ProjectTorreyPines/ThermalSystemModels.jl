using ModelingToolkit, Symbolics, Printf, Graphs, MetaGraphs
MTK = ModelingToolkit


function gas_connect(pins...)

    eqs = [
        sum(pin -> pin.ṁ, pins) ~ 0.0, # Mass
        sum(pin -> pin.Φ, pins) ~ 0.0, # Energy
    ]

    for pin in pins[2:end]
        push!(eqs, pins[1].P ~ pin.P)   # p_1 = p_2, p_1 = p_3,...
        push!(eqs, pins[1].T ~ pin.T)   # p_1 = p_2, p_1 = p_3,...
    end
    return eqs
end

function work_connect(pins...)
    eqs = [
        sum(pin -> pin.Ẇ, pins) ~ 0.0, # Mass
    ]
    return eqs
end

function heat_connect(pins...)
    eqs = [
        sum(pin -> pin.Q̇, pins) ~ 0.0, # Mass
    ]
    return eqs
end

function showsol(c, sol::ODESolution)
    estorage = []
    for cel in c
        cvec = []
        nodematch = [:p, :n, :z, :y, :n1, :n2, :p1, :p2]

        for prop in nodematch
            if hasproperty(cel, prop)
                push!(cvec, getproperty(cel, prop))
            end
        end
        if (hasproperty(cel, :Ẇ) || hasproperty(cel, :Q̇)) &&
           (!hasproperty(cel, :p) && !hasproperty(cel, :n))
            push!(estorage, cel)
        end
        for n in cvec
            @printf "%-13s  P[bar]= %-8.1f T[K]= %-8.1f ṁ[kg/s]= %+-7.2f Φ=%+-7.2f" n.name sol[n.P][end] sol[n.T][end] sol[n.ṁ][end] sol[n.Φ][end] /
                                                                                                                                       10^6
            if hasproperty(n, :h)
                @printf "h[kJ/kg]= %-8.1f" sol[n.h][end] / 1e3
            end
            if hasproperty(n, :x)
                @printf "x= %-5.2f" sol[n.x][end]
            end

            if hasproperty(cel, :Q̇)
                @printf "Q̇[kW]= %-8.2f" sol[cel.Q̇][end] / 10^3
                # @printf "Q̇[k]= %-8.2i" sol[cel.Q̇][end]
            elseif hasproperty(cel, :q)
                @printf "Q̇[kW]= %-8.2f" sol[cel.q.Q̇][end] / 10^3

                # @printf "Q̇[k]= %-8.2s" sol[cel.q.Q̇][end]
            end
            if hasproperty(cel, :yfrac)
                @printf "y= %-8.2f" sol[cel.yfrac][end]
            end

            if hasproperty(cel, :C)
                @printf "C = %-8.2f" sol[cel.C][end]
            end

            if hasproperty(cel, :Ẇ)
                @printf "Ẇ[kW]= %-8.2f" sol[cel.Ẇ][end] / 1e3
            elseif hasproperty(cel, :w)
                @printf "Ẇ[kW]= %-8.2f" sol[cel.w.Ẇ][end] / 1e3
            end
            @printf "\n"
        end
    end
    @printf "\n"
    sumE = 0.0
    for es in estorage
        if hasproperty(es, :Ẇ)
            @printf "%-8s\t Ẇ [kW] = %.2f" es.name sol[es.Ẇ][end] / 10^3
            sumE += sol[es.Ẇ][end]
        elseif hasproperty(es, :Q̇)
            @printf "%-8s\t Q̇ [kW] = %.2f" es.name sol[es.Q̇][end] / 10^3
            sumE += sol[es.Q̇][end]
        end
        @printf "\n"
    end
    @printf "\n"
    @printf "Sum of Energy Storage: [kW] = %.2f\n" sumE / 10^3
end
function variable2symbol(var)
    sm = Symbol[]
    for i = 1:length(var)
        push!(sm, Symbolics.tosymbol(var[i]; escape = false))
    end
    return sm
end

function showsol(c, prob::ODEProblem; prob_attr = nothing)
    estorage = []
    evars = []
    if isnothing(prob_attr)
        prob_attr = (prob.u0, prob.p, 0.0)
    end

    sol(varr) = prob.f.observed(varr, prob_attr...)
    for cel in c
        cvec = []
        nodematch = [:p, :n, :z, :y, :n1, :n2, :p1, :p2]

        for prop in nodematch
            if hasproperty(cel, prop)
                push!(cvec, getproperty(cel, prop))
            end
        end
        if (hasproperty(cel, :Ẇ) || hasproperty(cel, :Q̇)) &&
           (!hasproperty(cel, :p) && !hasproperty(cel, :n))
            push!(estorage, cel)
        end
        for n in cvec
            sname = String(n.name)
            if length(sname) > 13
                sname = sname[1:3] * "..." * sname[end-8:end]
            end
            push!(evars, n.Φ)
            push!(evars, n.ṁ)
            @printf "%-13s  P[bar]= %-8.1f T[K]= %-8.1f ṁ[kg/s]= %+-7.2f Φ=%+-8.2f" sname sol(
                n.P,
            ) sol(n.T) sol(n.ṁ) sol(n.Φ) / 10^6
            if hasproperty(n, :h)
                @printf "h[kJ/kg]= %-8.1f" sol(n.h) / 1e3
            end
            if hasproperty(n, :x)
                @printf "x= %-5.2f" sol(n.x)
            end

            # if hasproperty(n,:v)
            #     @printf "v[g/cm3]= %-5.2f" sol(n.v)*1e3
            # end

            # if hasproperty(n,:cp)
            #     @printf "cp[g/cm3]= %-5.2f" sol(n.cp)
            # end

            if hasproperty(cel, :Q̇)
                @printf "Q̇[kW]= %-8.2f" sol(cel.Q̇) / 10^3
                push!(evars, cel.Q̇)
                hasproperty(cel, :q) ? push!(evars, cel.q.Q̇) : nothing
            elseif hasproperty(cel, :q)
                @printf "Q̇[kW]= %-8.2f" sol(cel.q.Q̇) / 10^3
                push!(evars, cel.q.Q̇)
            end
            if hasproperty(cel, :yfrac)
                @printf "y= %-8.2f" sol(cel.yfrac)
            end

            if hasproperty(cel, :C)
                @printf " C= %-8.2f" sol(cel.C)
            end

            if hasproperty(cel, :Ẇ)
                @printf "Ẇ[kW]= %-8.2f" sol(cel.Ẇ) / 1e3
                push!(evars, cel.Ẇ)
                hasproperty(cel, :w) ? push!(evars, cel.w.Ẇ) : nothing
            elseif hasproperty(cel, :w)
                @printf "Ẇ[kW]= %-8.2f" sol(cel.w.Ẇ) / 1e3
                push!(evars, cel.w.Ẇ)
            end
            @printf "\n"
        end
    end
    @printf "\n"
    sumE = 0.0
    for es in estorage
        if hasproperty(es, :Ẇ)
            @printf "%-8s\t Ẇ [kW] = %.2f" es.name sol(es.Ẇ) / 10^3
            sumE += sol(es.Ẇ)
            push!(evars, es.Ẇ)
        elseif hasproperty(es, :Q̇)
            @printf "%-8s\t Q̇ [kW] = %.2f" es.name sol(es.Q̇) / 10^3
            sumE += sol(es.Q̇)
            push!(evars, es.Q̇)
        end
        @printf "\n"
    end
    @printf "\n"
    @printf "Sum of Energy Storage: [kW] = %.2f\n" sumE / 10^3

    evar_sym = variable2symbol(evars)
    evar_dict = Dict(evar_sym[i] => evars[i] for i = 1:length(evar_sym))
    evar_string_dict = Dict(String(evar_sym[i]) => evars[i] for i = 1:length(evar_sym))
    return evar_dict, evar_string_dict
end

function showSpecificSol(c, sol)
    estorage = []
    for cel in c
        cvec = []
        if hasproperty(cel, :p) && hasproperty(cel, :w)
            @printf "%-13s  Ẇ[kW]=%+-8.2f ṁ[kg/s]= %+-7.2f ẇ[kJ/kg] = %+-7.2f\n" cel.name sol[cel.w.Ẇ][end] /
                                                                                             1e3 sol[cel.p.ṁ][end] sol[cel.w.Ẇ][end] /
                                                                                                                    sol[cel.p.ṁ][end] /
                                                                                                                    1e3
        end
    end

end

function showsys(sys)
    display([s.name for s in sys.systems])
    return [s for s in sys.systems]
end

function sys2dict(sys::Vector{ODESystem})
    d = Dict()
    [d[s.name] = s for s in sys]
    return d
end

function sys2dict(sys::ODESystem)
    d = Dict()
    [d[s.name] = s for s in sys.systems]
    return d
end

function plantsys(Plant; compoenent_level = 2)
    d = Dict()
    for s in Plant.systems
        for ss in s.systems
            id = String(s.name) * "." * String(ss.name)
            d[id] = ss
        end
    end
    return d
end

function systems(p::ODESystem)
    return p.systems
end

#   returns all connection equations
function connection_equations(BaseSys)
    equat = equations(BaseSys)
    sysEq = Equation[]
    for sys in BaseSys.systems
        push!(sysEq, equations(sys)...)
    end
    return setdiff(equat, sysEq)
end


"""
    check_for_var(var; tofind = :ṁ)

    checks if a an array of variables are "tofind" variables 
    returns false is any of the input variables are not the "tofind" variable
    Useful for identifying balance equations within a larger system equations vector
        example mass flow continuity:  ṁ_in + ṁ_out = 0
        only variables in that equations are ṁ
"""
function check_for_var(var; tofind = :ṁ)

    strnam = [String(s) for s in var]

    tfStr = String(tofind)
    flowvar = collect(tfStr)[1]
    #   finds all system variables with ṁ
    for (i, str) in enumerate(strnam)
        chr = collect(str)[end]
        if chr != flowvar
            return false
        end
    end
    return true
end
check_for_ṁ(var) = check_for_var(var; tofind = :ṁ)
check_for_Ẇ(var) = check_for_var(var; tofind = :Ẇ)
check_for_Q(var) = check_for_var(var; tofind = Symbol("̇"))


"""
    find_flow_vars(BaseSys; tofind = :ṁ)

    this returns all system variables with mass flow

"""
function find_flow_vars(BaseSys; tofind = :ṁ)
    v2n = BaseSys.var_to_name
    symnam = collect(keys(BaseSys.var_to_name))
    strnam = [String(s) for s in symnam]
    tfStr = String(tofind)
    flowvar = collect(tfStr)[1]
    flowkeys = Symbol[]
    flowvals = SymbolicUtils.BasicSymbolic{Real}[]
    #   finds all system variables with ṁ
    for (i, str) in enumerate(strnam)
        chr = collect(str)[end]
        if chr == flowvar
            push!(flowkeys, symnam[i])
            push!(flowvals, v2n[symnam[i]])
        end
    end
    return v2n, flowkeys, flowvals
end
mass_flow_vars(x) = find_flow_vars(x; tofind = :ṁ)
work_flow_vars(x) = find_flow_vars(x; tofind = :Ẇ)
heat_flow_vars(x) = find_flow_vars(x; tofind = "̇") # this is odd, took forever to figure this out



"""
    find_flow_vars(BaseSys; tofind = :ṁ)
    :divertor_circulator₊p₊ṁ => ["divertor_circulator", "p", "ṁ"]
"""
function split_variable_sym(var)
    return split(String(Symbolics.tosymbol(var; escape = false)), "₊")
end

function join_variable_str(splitvar, num2keep::Int; joinstr = "₊")
    if num2keep > length(splitvar)
        return join(splitvar, joinstr)
    else
        return join(splitvar[1:num2keep], joinstr)
    end
end

function component_names(eqs; unique_only = true, num2keep = 1)
    compnames = Symbol[]
    for eq in eqs
        vars = ModelingToolkit.get_variables(eq)
        symvar = variable2symbol(vars) # variables in the form :variablenmame
        bol = check_for_ṁ(symvar)
        if bol == true
            # cname = [split_variable_sym(s)[1] for s in symvar]
            cname = [join_variable_str(split_variable_sym(s), num2keep) for s in symvar]
            for cn in cname
                push!(compnames, Symbol(cn))
            end
        end
    end
    unique_only ? compnames = unique(compnames) : nothing
    return compnames
end

function component_names(eqs, tofind; unique_only = true)
    compnames = Symbol[]
    for eq in eqs
        vars = ModelingToolkit.get_variables(eq)
        symvar = variable2symbol(vars) # variables in the form :variablenmame
        bol = check_for_var(symvar; tofind)
        if bol == true
            cname = [split_variable_sym(s)[1] for s in symvar]
            for cn in cname
                push!(compnames, Symbol(cn))
            end
        end
    end
    unique_only ? compnames = unique(compnames) : nothing
    return compnames
end

function all_component_names(sys; num2keep = 1, unique_only = true)
    long_names = collect(keys(sys.var_to_name))
    compnames = Symbol[]
    for ln in long_names
        symvar = variable2symbol([ln]) # variables in the form :variablenmame
        cname = join_variable_str(split_variable_sym(symvar[1]), num2keep)
        push!(compnames, Symbol(cname))
    end
    unique_only ? compnames = unique(compnames) : nothing
    return compnames
end

function system_details(sys::ODESystem; alias_elimate = true)
    @printf "************ SYSTEM: %s ************\n" sys.name
    @printf "SYSTEM: %s \n" sys.name
    @printf "\t # of equations = %i\n" length(equations(sys))
    @printf "\t # of states = %i\n" length(states(sys))
    @printf "\t # of parameters = %i\n" length(parameters(sys))
    if alias_elimate
        @printf "After Alias Elemination:\n"
        al_sys = alias_elimination(sys)
        @printf "\t # of equations = %i\n" length(equations(al_sys))
        @printf "\t # of states = %i\n" length(states(al_sys))
        @printf "\t # of parameters = %i\n" length(parameters(al_sys))
    end
end

function connection_equations2(BaseSys; flagvec = [])
    equat = equations(BaseSys)
    sysEq = Equation[]
    compnames = Symbol[]
    for eq in equat
        vars = MTK.get_variables(eq)
        symvar = variable2symbol(vars) # variables in the form :variablenmame
        cname = [split_variable_sym(s) for s in symvar]

        # checking n vs p node
        # portid = [String(c[end-1]) for c in cname]
        # system name
        portNames = [String(c[1]) for c in cname]
        portSym = [Symbol(c[1]) for c in cname]
        # end variable
        endvars = [String(c[end]) for c in cname]

        for (i, ps) in enumerate(portSym)
            if portNames[i] != endvars[i]
                push!(compnames, ps)
            end
        end

        if length(portNames) == length(unique(portNames)) && length(unique(endvars)) == 1
            push!(sysEq, eq)
        elseif length(unique(endvars)) == 1
            if length(intersect(flagvec, portSym)) > 0
                push!(sysEq, eq)
            end
        end
    end
    return sysEq, unique(compnames)
end

function system2graph(ODESYS::ODESystem; verbose = false)
    eqs = connection_equations2(ODESYS)
    v2n, flowkeys, flowvals = mass_flow_vars(ODESYS)
    compnames = component_names(eqs; unique_only = true)

    # verbose ? println(compnames) : nothing
    component_name_dict = Dict(cc => ii for (ii, cc) in enumerate(compnames))

    verbose ? display(collect(keys(component_name_dict))) : nothing
    num_to_name_dict = Dict(ii => cc for (ii, cc) in enumerate(compnames))
    g = DiGraph()
    add_vertices!(g, length(compnames))

    # dict of variables 
    edge_power_dict = Dict()

    # loop through connection equations
    for eq in eqs
        vars = MTK.get_variables(eq)
        symvar = variable2symbol(vars) # variables in the form :variablenmame
        bol = check_for_ṁ(symvar)

        # bol means that mass flow is a used variable (flow connected)
        if bol == true
            # Splitting name in system
            # example ["reservoir", "n", "ṁ"], ["valve", "p", "ṁ"]
            cname = [split_variable_sym(s) for s in symvar]

            # checking n vs p node
            portid = [String(c[end-1]) for c in cname]

            # system name
            portNames = [String(c[1]) for c in cname]

            verbose ? println(portNames) : nothing

            src = findall(x -> contains(x, "n"), portid)    # Source nodes
            dst = findall(x -> contains(x, "p"), portid)    # Destination nodes

            if !isempty(src) && !isnothing(dst) && !(portNames[src] == portNames[dst])
                # adding edge from n -> p
                if length(src) > 1
                    println(
                        "Warning - may incorrectly create graph @ src = $(portNames[src]) to dst = $(portNames[dst])",
                    )
                end

                for s in src
                    for d in dst
                        verbose ? println(src) : nothing
                        if portNames[s] != portNames[d]
                            add_edge!(
                                g,
                                component_name_dict[Symbol(portNames[s])],
                                component_name_dict[Symbol(portNames[d])],
                            )
                            cp = cname[d]
                            cp[end] = "Φ"
                            pwr_id = Symbolics.tosymbol(join(cp, "₊"); escape = false)
                            edge_power_dict[(
                                component_name_dict[Symbol(portNames[s])],
                                component_name_dict[Symbol(portNames[d])],
                            )] = Symbol(pwr_id)
                            # Symbolics.tosymbol(replace(String(symvar[d]), "ṁ" => "Φ"); escape = false)
                        end
                    end
                end
            end
        end
    end
    return g, component_name_dict, num_to_name_dict, edge_power_dict
end

function system2graph2(
    ODESYS::ODESystem,
    utility_vector::Vector{Symbol};
    verbose = false,
    soln = nothing,
)
    eqs, compnames = connection_equations2(ODESYS)

    verbose ? println(compnames) : nothing

    component_name_dict = Dict(cc => ii for (ii, cc) in enumerate(compnames))
    num_to_name_dict = Dict(ii => cc for (ii, cc) in enumerate(compnames))

    num_to_sys_dict = Dict()
    for sys in ODESYS.systems
        num = component_name_dict[sys.name]
        num_to_sys_dict[num] = sys
    end


    g = DiGraph()
    add_vertices!(g, length(compnames))

    # dict of variables 
    edge_power_dict = Dict()

    # loop through connection equations
    for eq in eqs
        vars = MTK.get_variables(eq)
        symvar = variable2symbol(vars) # variables in the form :variablenmame
        # bol means that mass flow is a used variable (flow connected)
        if check_for_ṁ(symvar)
            # Splitting name in system
            # example ["reservoir", "n", "ṁ"], ["valve", "p", "ṁ"]
            cname = [split_variable_sym(s) for s in symvar]

            # checking n vs p node
            portid = [String(c[end-1]) for c in cname]

            # system name
            portNames = [String(c[1]) for c in cname]

            verbose ? println(portNames) : nothing

            src = findall(x -> contains(x, "n"), portid)    # Source nodes
            dst = findall(x -> contains(x, "p"), portid)    # Destination nodes

            if !isempty(src) && !isnothing(dst) && !(portNames[src] == portNames[dst])
                # adding edge from n -> p
                if length(src) > 1
                    println(
                        "Warning - may incorrectly create graph @ src = $(portNames[src]) to dst = $(portNames[dst])",
                    )
                end

                for s in src
                    for d in dst
                        verbose ? println(src) : nothing
                        if portNames[s] != portNames[d]
                            add_edge!(
                                g,
                                component_name_dict[Symbol(portNames[s])],
                                component_name_dict[Symbol(portNames[d])],
                            )
                            cp = cname[s]
                            cp[end] = "Φ"
                            pwr_id = Symbolics.tosymbol(join(cp, "₊"); escape = false)
                            edge_power_dict[(
                                component_name_dict[Symbol(portNames[s])],
                                component_name_dict[Symbol(portNames[d])],
                            )] = Symbol(pwr_id)
                            # Symbolics.tosymbol(replace(String(symvar[d]), "ṁ" => "Φ"); escape = false)
                        end
                    end
                end
            end
        elseif check_for_Q(symvar)
            #consumer flags
            pflags = ["cold", "env", "waste"]

            cname = [split_variable_sym(s) for s in symvar]
            nconnected = length(symvar)
            portNames = [Symbol(c[1]) for c in cname]
            non_util = setdiff(portNames, utility_vector)

            if length(non_util) != length(portNames)
                # means there is a utility in the equation
                util_component = intersect(portNames, utility_vector)[1]
                verbose ? println(util_component) : nothing
                util_src = true
                for flg in pflags
                    verbose ?
                    println(
                        "  $(flg) => $(util_component) $( occursin(lowercase(String(flg)),lowercase(String(util_component))))",
                    ) : nothing
                    if occursin(lowercase(String(flg)), lowercase(String(util_component)))
                        util_src = false
                    end
                end

                for (i, port) in enumerate(non_util)
                    src = util_component
                    dst = port

                    idx = findfirst(x -> x == port, portNames)
                    # checking if direction needs to be switched
                    if util_src == false
                        verbose ? println(util_component) : nothing
                        src = port
                        dst = util_component
                    end
                    svar = symvar[idx]
                    add_edge!(g, component_name_dict[src], component_name_dict[dst])
                    edge_power_dict[(component_name_dict[src], component_name_dict[dst])] =
                        Symbolics.tosymbol(svar)
                end


            elseif length(portNames) == 2
                # direct heat exchange

                dst = portNames[2]
                src = portNames[1]
                svar = symvar[2]
                if !isnothing(soln)
                    if soln(svar) < 0
                        dst = portNames[1]
                        svar = symvar[1]
                        src = portNames[2]
                    end
                end

                add_edge!(g, component_name_dict[src], component_name_dict[dst])
                edge_power_dict[(component_name_dict[src], component_name_dict[dst])] =
                    Symbolics.tosymbol(svar)
            end

        elseif check_for_Ẇ(symvar)
            #consumer flags
            pflags = ["pump", "circulator", "compr"]
            # producer flags
            tflags = ["turb"]

            cname = [split_variable_sym(s) for s in symvar]
            nconnected = length(symvar)
            portNames = [Symbol(c[1]) for c in cname]

            non_util = setdiff(portNames, utility_vector)

            if length(non_util) != length(portNames)
                # means there is a utility in the equation
                util_component = intersect(portNames, utility_vector)

                for (i, port) in enumerate(non_util)
                    src = util_component[1]
                    dst = port

                    idx = findfirst(x -> x == port, portNames)
                    # checking if direction needs to be switched
                    for flg in tflags
                        if occursin(String(flg), String(port))
                            src = port
                            dst = util_component[1]
                        end
                    end
                    verbose ? display(src) : nothing
                    add_edge!(g, component_name_dict[src], component_name_dict[dst])
                    edge_power_dict[(component_name_dict[src], component_name_dict[dst])] =
                        Symbolics.tosymbol(symvar[idx])
                end
            end
        end
    end

    # edge_pow_verts = collect(keys(edge_power_dict))
    # edge_pow_vars = collect(values(edge_power_dict));
    # edge_pow_var_sym = [Symbol(join(split_variable_sym(var),".")) for var in edge_pow_vars]

    # revised_edge_power_dict = Dict(edge_pow_verts[i] => edge_pow_var_sym[i] for i = 1:length(edge_pow_vars))


    return g, component_name_dict, num_to_name_dict, edge_power_dict, num_to_sys_dict
end

function add_graph_connection!(
    sys,
    g,
    component_dict,
    vertex_dict,
    Component_Variable,
    edge_power_dict;
    flip_dir = false,
    verbose = false,
)
    comp_names = collect(keys(component_dict))
    eqs = equations(sys)
    parsed_name = split_variable_sym(Component_Variable)
    element_name = parsed_name[1]
    element_var = parsed_name[end]
    str_el = String(element_name)

    if element_var == "Q̇"
        element_var = "̇"
    end

    deps = MTK.variable_dependencies(sys; variables = Component_Variable)
    eqs_to_search = eqs[deps.fadjlist...]
    eq_names = component_names(eqs_to_search, element_var; unique_only = true)
    new_work_names = setdiff(eq_names, comp_names)

    for n in new_work_names
        add_vertex!(g)
        component_dict[n] = nv(g)
        vertex_dict[nv(g)] = n
    end


    for eq in eqs_to_search
        vars = MTK.get_variables(eq)
        symvar = variable2symbol(vars)
        parsed_elem = [split_variable_sym(s) for s in symvar]
        parsed_names = [String(s[1]) for s in parsed_elem]
        non_input = setdiff(parsed_names, str_el)

        for (i, ni) in enumerate(non_input)
            # default = reservoir -> component
            src = Symbol(element_name)
            dst = Symbol(ni)

            #component -> reservoir
            if occursin("turb", ni) ||
               occursin("conde", ni) ||
               occursin("Rej", ni) ||
               flip_dir == true
                # # always from the hot utility
                # @show "switch on $(ni) -> $(element_name)"
                src = Symbol(ni)
                dst = Symbol(element_name)
            end
            if src != dst
                # display(eq)
                display(vars)
                if verbose
                    println("Added $(ni) -> $(element_name)")
                    println("    Idx $(component_dict[src]) -> $(component_dict[dst])")
                    println(
                        "    var $(vertex_dict[component_dict[src]]) -> $(vertex_dict[component_dict[dst]])",
                    )
                end
                add_edge!(g, component_dict[src], component_dict[dst])
                edge_power_dict[(component_dict[src], component_dict[dst])] =
                    Symbolics.tosymbol(symvar[i])
            end
        end
    end
end

function system2metagraph(
    sys::ODESystem,
    utility_vector::Vector{Symbol};
    verbose = false,
    soln = nothing,
)
    # initialize component names and relation equations
    eqs, compnames = connection_equations2(sys; flagvec = utility_vector)

    G = MetaDiGraph()
    add_vertices!(G, length(compnames))

    # Adding system property
    set_prop!(G, :sys, sys)
    set_prop!(G, :utility_vector, utility_vector)
    set_prop!(G, :soln, soln)

    component_name_dict = Dict(cc => ii for (ii, cc) in enumerate(compnames))
    num_to_name_dict = Dict(ii => cc for (ii, cc) in enumerate(compnames))
    num_to_sys_dict = Dict()
    for s in sys.systems
        num = component_name_dict[s.name]
        num_to_sys_dict[num] = s
        set_prop!(G, num, :name, s.name)# 
        set_prop!(G, num, :sys, s)
    end
    set_indexing_prop!(G, :name)
    # loop through connection equations
    for eq in eqs
        vars = MTK.get_variables(eq)
        symvar = variable2symbol(vars) # variables in the form :variablenmame
        # bol means that mass flow is a used variable (flow connected)
        if check_for_ṁ(symvar)
            # Splitting name in system
            # example ["reservoir", "n", "ṁ"], ["valve", "p", "ṁ"]
            cname = [split_variable_sym(s) for s in symvar]

            # checking n vs p node
            portid = [String(c[end-1]) for c in cname]

            # system name
            portNames = [String(c[1]) for c in cname]

            # verbose ? println(portNames) : nothing

            src = findall(x -> contains(x, "n"), portid)    # Source nodes
            dst = findall(x -> contains(x, "p"), portid)    # Destination nodes

            if !isempty(src) && !isnothing(dst) && !(portNames[src] == portNames[dst])
                # adding edge from n -> p
                if length(src) > 1
                    println(
                        "Warning - may incorrectly create graph @ src = $(portNames[src]) to dst = $(portNames[dst])",
                    )
                end

                for s in src
                    for d in dst
                        # verbose ? println(src) : nothing
                        if portNames[s] != portNames[d]
                            add_edge!(
                                G,
                                component_name_dict[Symbol(portNames[s])],
                                component_name_dict[Symbol(portNames[d])],
                            )
                            cp = cname[s]
                            cp[end] = "Φ"
                            pwr_id = Symbolics.tosymbol(join(cp, "₊"); escape = false)
                            # verbose ? println(cp) : nothing
                            dsys =
                                get_prop(G, component_name_dict[Symbol(portNames[d])], :sys)
                            verbose ? println("$(Symbol(portNames[s])) => $(dsys.name)") :
                            nothing
                            # verbose ? println("$(dsys.name) $(typeof(portid[d]))") : nothing
                            set_prop!(
                                G,
                                Edge(
                                    component_name_dict[Symbol(portNames[s])],
                                    component_name_dict[Symbol(portNames[d])],
                                ),
                                :etype,
                                :flow,
                            )
                            if portid[d] == "p"
                                set_prop!(
                                    G,
                                    Edge(
                                        component_name_dict[Symbol(portNames[s])],
                                        component_name_dict[Symbol(portNames[d])],
                                    ),
                                    :evar,
                                    dsys.p.Φ,
                                )
                                set_prop!(
                                    G,
                                    Edge(
                                        component_name_dict[Symbol(portNames[s])],
                                        component_name_dict[Symbol(portNames[d])],
                                    ),
                                    :mflow,
                                    dsys.p.ṁ,
                                )
                            elseif portid[d] == "p1"
                                set_prop!(
                                    G,
                                    Edge(
                                        component_name_dict[Symbol(portNames[s])],
                                        component_name_dict[Symbol(portNames[d])],
                                    ),
                                    :evar,
                                    dsys.p1.Φ,
                                )
                                set_prop!(
                                    G,
                                    Edge(
                                        component_name_dict[Symbol(portNames[s])],
                                        component_name_dict[Symbol(portNames[d])],
                                    ),
                                    :mflow,
                                    dsys.p1.ṁ,
                                )
                            elseif portid[d] == "p2"
                                set_prop!(
                                    G,
                                    Edge(
                                        component_name_dict[Symbol(portNames[s])],
                                        component_name_dict[Symbol(portNames[d])],
                                    ),
                                    :evar,
                                    dsys.p2.Φ,
                                )
                                set_prop!(
                                    G,
                                    Edge(
                                        component_name_dict[Symbol(portNames[s])],
                                        component_name_dict[Symbol(portNames[d])],
                                    ),
                                    :mflow,
                                    dsys.p2.ṁ,
                                )
                            else
                                println(
                                    "Failed to find destination node name for $(Symbol(portNames[s])) => $(dsys.name) , add key $(portid[d])",
                                )
                            end
                        end
                    end
                end
            end
        elseif check_for_Q(symvar)
            verbose ? println(" Heat Transfer variables FLAG") : nothing
            #consumer flags
            pflags = ["cold", "env", "waste"]

            cname = [split_variable_sym(s) for s in symvar]
            nconnected = length(symvar)
            portNames = [Symbol(c[1]) for c in cname]
            non_util = setdiff(portNames, utility_vector)

            if length(non_util) != length(portNames)
                # means there is a utility in the equation
                util_component = intersect(portNames, utility_vector)[1]
                verbose ? println(util_component) : nothing
                util_src = true
                for flg in pflags
                    verbose ?
                    println(
                        "  $(flg) => $(util_component) $( occursin(lowercase(String(flg)),lowercase(String(util_component))))",
                    ) : nothing
                    if occursin(lowercase(String(flg)), lowercase(String(util_component)))
                        util_src = false
                    end
                end

                for (i, port) in enumerate(non_util)
                    src = util_component
                    dst = port
                    etyp = :externalheat

                    idx = findfirst(x -> x == port, portNames)
                    # checking if direction needs to be switched
                    if util_src == false
                        verbose ? println(util_component) : nothing
                        src = port
                        dst = util_component
                        etyp = :externalcool
                    end
                    svar = vars[idx]
                    add_edge!(G, component_name_dict[src], component_name_dict[dst])
                    set_prop!(
                        G,
                        Edge(component_name_dict[src], component_name_dict[dst]),
                        :etype,
                        etyp,
                    )
                    set_prop!(
                        G,
                        Edge(component_name_dict[src], component_name_dict[dst]),
                        :evar,
                        svar,
                    )
                end


            elseif length(portNames) == 2
                # direct heat exchange

                dst = portNames[2]
                src = portNames[1]
                svar = vars[2]
                if !isnothing(soln)
                    if soln(vars[2]) < 0
                        dst = portNames[1]
                        svar = vars[1]
                        src = portNames[2]
                    end
                end
                add_edge!(G, component_name_dict[src], component_name_dict[dst])
                set_prop!(
                    G,
                    Edge(component_name_dict[src], component_name_dict[dst]),
                    :etype,
                    :transferheat,
                )
                set_prop!(
                    G,
                    Edge(component_name_dict[src], component_name_dict[dst]),
                    :evar,
                    svar,
                )
            end

        elseif check_for_Ẇ(symvar)
            verbose ? println(" Working variables FLAG") : nothing
            #consumer flags
            pflags = ["pump", "circulator", "compr"]
            # producer flags
            tflags = ["turb"]

            cname = [split_variable_sym(s) for s in symvar]
            nconnected = length(symvar)
            portNames = [Symbol(c[1]) for c in cname]

            non_util = setdiff(portNames, utility_vector)

            if length(non_util) != length(portNames)
                # means there is a utility in the equation
                util_component = intersect(portNames, utility_vector)

                for (i, port) in enumerate(non_util)
                    src = util_component[1]
                    dst = port

                    idx = findfirst(x -> x == port, portNames)
                    # checking if direction needs to be switched
                    for flg in tflags
                        if occursin(String(flg), String(port))
                            src = port
                            dst = util_component[1]
                        end
                    end
                    # verbose ? display(src) : nothing
                    # verbose ? println("$(Symbol(portNames[src])) => $(portNames[dst])") : nothing
                    add_edge!(G, component_name_dict[src], component_name_dict[dst])
                    set_prop!(
                        G,
                        Edge(component_name_dict[src], component_name_dict[dst]),
                        :etype,
                        :electric,
                    )
                    set_prop!(
                        G,
                        Edge(component_name_dict[src], component_name_dict[dst]),
                        :evar,
                        vars[idx],
                    )
                end
            end
        end
    end
    return G
end

