using ModelingToolkit, Symbolics, Printf, Graphs
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

function showsol(c,sol)
    estorage = []
    for cel in c
        cvec = []
        nodematch = [:p, :n, :z, :y, :n1, :n2, :p1, :p2]

        for prop in nodematch
            if hasproperty(cel,prop)
                push!(cvec,getproperty(cel,prop))
            end
        end
        if (hasproperty(cel,:Ẇ) || hasproperty(cel,:Q̇)) && (!hasproperty(cel,:p) && !hasproperty(cel,:n))
            push!(estorage,cel)
        end
        for n in cvec
            @printf "%-13s  P[bar]= %-8.1f T[K]= %-8.1fṁ[kg/s]= %+-7.2fΦ=%+-7.2f" n.name sol[n.P][end] sol[n.T][end] sol[n.ṁ][end] sol[n.Φ][end]/10^6
            if hasproperty(n,:h)
                @printf "h[kJ/kg]= %-8.1f" sol[n.h][end]/1e3
            end
            if hasproperty(n,:x)
                @printf "x= %-5.2f" sol[n.x][end] 
            end

            if hasproperty(cel,:Q̇)
                @printf "Q̇[kW]= %-8.2f" sol[cel.Q̇][end]/10^3
            elseif hasproperty(cel,:q)
                @printf "Q̇[kW]= %-8.2f" sol[cel.q.Q̇][end]/10^3
            end
            if hasproperty(cel,:yfrac)
                @printf "y= %-8.2f" sol[cel.yfrac][end] 
            end

            if hasproperty(cel,:C)
                @printf "C= %.2e " sol[cel.C][end] 
            end

            if hasproperty(cel, :Ẇ)
                @printf "Ẇ[kW]= %-8.2f"  sol[cel.Ẇ][end]/1e3
            elseif hasproperty(cel,:w)
                @printf "Ẇ[kW]= %-8.2f" sol[cel.w.Ẇ][end]/1e3
            end
            @printf "\n"
        end
    end
    @printf "\n"
    for es in estorage
        if hasproperty(es, :Ẇ)
            @printf "%-8s\t Ẇ [kW] = %.2f" es.name sol[es.Ẇ][end]/10^3
        elseif hasproperty(es,:Q̇)
            @printf "%-8s\t Q̇ [kW] = %.2f" es.name sol[es.Q̇][end]/10^3
        end
        @printf "\n"
    end
    @printf "\n"
end

function showSpecificSol(c,sol)
    estorage = []
    for cel in c
        cvec = []
        if hasproperty(cel,:p) && hasproperty(cel,:w)
            @printf "%-13s  Ẇ[kW]=%+-8.2f ṁ[kg/s]= %+-7.2f ẇ[kJ/kg] = %+-7.2f\n" cel.name  sol[cel.w.Ẇ][end]/1e3 sol[cel.p.ṁ][end] sol[cel.w.Ẇ][end]/sol[cel.p.ṁ][end]/1e3
        end
    end

end

function showsys(sys)
    display([s.name for s in sys.systems])
    return [s for s in sys.systems];
end

function sys2dict(sys)
    d = Dict()
    [d[s.name] = s for s in sys.systems];
    return d
end

function plantsys(Plant; compoenent_level = 2)
    d = Dict()
    for s in Plant.systems
        for ss in s.systems
            id = String(s.name)*"."*String(ss.name);
            d[id] = ss;
        end
    end
    return d
end

function systems(p::ODESystem)  
    return p.systems
end

function substitute_props(sys::ODESystem,subdicts)
    systems = ODESystem[]
end

#   returns all connection equations
function connection_equations(BaseSys)
    equat = equations(BaseSys)
    sysEq = Equation[]
    for sys in BaseSys.systems
        push!(sysEq,equations(sys)...)
    end
    return setdiff(equat,sysEq)
end

function check_for_var(var; tofind = :ṁ)

    strnam = [String(s) for s in var]

    tfStr   = String(tofind)
    flowvar = collect(tfStr)[1]
    #   finds all system variables with ṁ
    for (i,str) in enumerate(strnam)
        chr = collect(str)[end]
        if chr != flowvar
            return false
        end
    end
    return true  
end

check_for_ṁ(var) = check_for_var(var;tofind = :ṁ)
check_for_Ẇ(var) = check_for_var(var; tofind = :Ẇ)
check_for_Q(var) = check_for_var(var; tofind = Symbol("̇"))
# this returns all system variables with mass flow
function find_flow_vars(BaseSys; tofind = :ṁ)
    v2n = BaseSys.var_to_name
    symnam = collect(keys(BaseSys.var_to_name))
    strnam = [String(s) for s in symnam]
    tfStr = String(tofind)
    flowvar = collect(tfStr)[1]
    flowkeys = Symbol[]
    flowvals = SymbolicUtils.BasicSymbolic{Real}[]
    #   finds all system variables with ṁ
    for (i,str) in enumerate(strnam)
        chr = collect(str)[end]
        if chr == flowvar
            push!(flowkeys,symnam[i])
            push!(flowvals,v2n[symnam[i]])
        end
    end
    return v2n,flowkeys, flowvals
end
mass_flow_vars(x) = find_flow_vars(x; tofind = :ṁ)
work_flow_vars(x) = find_flow_vars(x; tofind = :Ẇ )
heat_flow_vars(x) = find_flow_vars(x; tofind = "̇" ) # this is odd, took forever to figure this out

function variable2symbol(var)
    sm = Symbol[]
    for i = 1:length(var)
        push!(sm,Symbolics.tosymbol(var[i]; escape = false))
    end
    return sm

end

function split_variable_sym(var)
    return split(String(Symbolics.tosymbol(var; escape = false)), "₊")
end

function join_variable_str(splitvar, num2keep::Int)
    if num2keep > length(splitvar)
        return join(splitvar,"₊")
    else
        return join(splitvar[1:num2keep],"₊")
    end
end

function component_names(eqs; unique_only = true, num2keep = 1)
    compnames = Symbol[]
    for eq in eqs
        vars = ModelingToolkit.get_variables(eq)
        symvar = variable2symbol(vars); # variables in the form :variablenmame
        bol = check_for_ṁ(symvar)
        if bol == true
            # cname = [split_variable_sym(s)[1] for s in symvar]
            cname = [join_variable_str(split_variable_sym(s),num2keep) for s in symvar]
            for cn in cname
                push!(compnames,Symbol(cn))
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
        symvar = variable2symbol(vars); # variables in the form :variablenmame
        bol = check_for_var(symvar; tofind)
        if bol == true
            cname = [split_variable_sym(s)[1] for s in symvar]
            for cn in cname
                push!(compnames,Symbol(cn))
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
        symvar = variable2symbol([ln]); # variables in the form :variablenmame
        cname = join_variable_str(split_variable_sym(symvar[1]),num2keep)
        push!(compnames,Symbol(cname))
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

function system2graph(ODESYS::ODESystem; verbose = false)
    eqs = connection_equations(ODESYS);
    v2n,flowkeys, flowvals = mass_flow_vars(ODESYS)
    compnames = component_names(eqs; unique_only = true)

    verbose ? println(compnames) : nothing
    component_name_dict = Dict(cc => ii for (ii,cc) in enumerate(compnames))
    num_to_name_dict = Dict(ii => cc for (ii,cc) in enumerate(compnames))
    g = DiGraph()
    add_vertices!(g,length(compnames))

    # dict of variables 
    edge_power_dict = Dict()

    for eq in eqs
        vars    = MTK.get_variables(eq)
        symvar  = variable2symbol(vars); # variables in the form :variablenmame
        bol     = check_for_ṁ(symvar)
        if bol == true
            # Splitting name in system
            # example ["reservoir", "n", "ṁ"], ["valve", "p", "ṁ"]
            cname = [split_variable_sym(s) for s in symvar]

            # checking n vs p node
            portid = [String(c[end-1]) for c in cname]

            # system name
            portNames = [String(c[1]) for c in cname]

            verbose ? println(portNames) : nothing

            src = findall(x -> contains(x, "n"), portid)
            dst = findall(x -> contains(x, "p"), portid)

            if !isempty(src) && !isnothing(dst) && !(portNames[src] == portNames[dst])
                # adding edge from n -> p
                if length(src) > 1
                    println("Warning - may incorrectly create graph @ src = $(portNames[src]) to dst = $(portNames[dst])")
                end

                for d in dst
                    verbose ? println(src) : nothing
                    if portNames[src[1]] != portNames[d]
                        add_edge!(g,component_name_dict[Symbol(portNames[src[1]])], component_name_dict[Symbol(portNames[d])])
                        cp = cname[d]
                        cp[end] = "Φ"
                        pwr_id = Symbolics.tosymbol(join(cp,"₊"); escape = false)
                        edge_power_dict[(component_name_dict[Symbol(portNames[src[1]])], component_name_dict[Symbol(portNames[d])])] = Symbol(pwr_id)
                        # Symbolics.tosymbol(replace(String(symvar[d]), "ṁ" => "Φ"); escape = false)
                    end
                end
            end
        end
    end
    return g, component_name_dict, num_to_name_dict, edge_power_dict
end

function add_graph_connection!(sys, g, component_dict, vertex_dict, Component_Variable,edge_power_dict)
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
        vars    = MTK.get_variables(eq)
        symvar  = variable2symbol(vars);
        parsed_elem = [split_variable_sym(s) for s in symvar]
        parsed_names = [String(s[1]) for s in parsed_elem]
        non_input = setdiff(parsed_names,str_el)

        for (i,ni) in enumerate(non_input)
            # default = reservoir -> component
            src = Symbol(element_name)
            dst = Symbol(ni)

            #component -> reservoir
            if occursin("turb", ni) || occursin("conde", ni) || occursin("Rej", ni)
                # # always from the hot utility
                # @show "switch on $(ni) -> $(element_name)"
                src = Symbol(ni)
                dst = Symbol(element_name)
            end
            if src !=  dst 
                add_edge!(g,component_dict[src], component_dict[dst]) 
                edge_power_dict[(component_dict[src], component_dict[dst])] = Symbolics.tosymbol(symvar[i])
            end
        end
    end
end

function meta_graph(Plant_system; level::Int = 2)
    flow_comps = component_names(Plant_system; num2keep = level)

    component_name_dict = Dict(cc => ii for (ii,cc) in enumerate(flow_comps))
    num_to_name_dict = Dict(ii => cc for (ii,cc) in enumerate(flow_comps))

    con_eq = connection_equations(Plant_system)

    v2n, mflow_syms, mflow_vars = mass_flow_vars(Plant_system)
    v2n, heat_syms, heat_vars = heat_flow_vars(Plant_system)
    v2n, work_syms, work_vars = work_flow_vars(Plant_system)

    all_eqs = equations(Plant_system)

end
