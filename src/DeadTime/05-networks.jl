##########################
# ============================================================================#
#                                                                             #
#                              NETWORK ANAYLYS
#
# =============================================================================
using DataFrames, SparseArrays, Printf
using NetworkLayout
using Graphs

# include("03-utilities.jl")
Base.@kwdef mutable struct component_network_params
    Tmin::Float64 = 650.0     # temperature
    Pmin::Float64 = 80.0      # pressure
    Tmax::Float64 = 900.0     # max allowable temperatures
    working_fluid::Symbol = :helium
end

mutable struct component_network
    params::component_network_params
    nodes::Vector{<:component_nodes}
    elements::Vector{<:model_component}
    mass_flow::Float64
    name::String
end

mutable struct system_network
    cycles::Vector{component_network}
    heat_exchangers::Vector{<:heat_exchanger}
    name::String
    eval_order
end

function show_params(cn::component_network_params)
    @printf "NETWORK: Tmin = %.2f, \t Pmin = %.2f,\t Tmax = %.2f\t working fluid = %s\n" cn.Tmin cn.Pmin cn.Tmax cn.working_fluid
end

function show_params(cn::component_network)
    show_params(cn.params)
end

breeder_network_params()        = component_network_params(700 + 273.15, 32, 1273.15, :pbli)
blanket_network_params()        = component_network_params(450 + 273.15, 80, 750+273.15, :helium)
divertor_network_params()       = component_network_params(350 + 273.15, 80, 650+273.15, :helium)

breeder_network_params2()       = component_network_params(550 + 273.15, 32, 750+273.15, :pbli)
blanket_network_params2()       = component_network_params(450 + 273.15, 80, 450+273.15, :helium)
divertor_network_params2()      = component_network_params(350 + 273.15, 80, 350+273.15, :helium)
extraction_network_params2()    = component_network_params(300 + 273.15, 80, 700+273.15, :helium)

function set_name(cn::component_network,in_string::String)
    cn.name = in_string
end

component_network(params::component_network_params,nodes::Vector{<:component_nodes},elements::Vector{<:model_component}) = component_network(params,nodes,elements,1.0,"noname")

function node_count(cn::component_network)
    return length(cn.nodes)
end

function elem_count(cn::component_network)
    return length(cn.elements)
end

function network2graph(cycle_network::component_network; verbose = false)
    num_edge    = elem_count(cycle_network) # number of elements
    num_vert    = node_count(cycle_network) # number of nodes (vertices)
    g           = DiGraph(num_vert)
    nms =[]
    edgelabel_dict = Dict()
    edgelabel_mat = Array{String}(undef, num_vert, num_vert)

    for el in cycle_network.elements
        push!(nms,el.name)
        if typeof(el.nodes) ==  siso_nodes
            i = el.nodes.A.node_num
            j = el.nodes.B.node_num
            add_edge!(g,i,j)
            edgelabel_mat[i, j] = edgelabel_dict[(i, j)] = el.name
            verbose ? println("Node added :  $(i) to $(j)  <=> $(el.name) ") : nothing
        elseif typeof(el.nodes) == miso_nodes
            #   mult input, single output
            j = el.nodes.B.node_num
            for in_node in el.nodes.A
                i = in_node.node_num
                add_edge!(g,i,j)
                edgelabel_mat[i, j] = edgelabel_dict[(i, j)] = el.name
                verbose ? println("Node added :  $(i) to $(j)  <=> $(el.name) ") : nothing
            end
        elseif typeof(el.nodes) == simo_nodes
            i = el.nodes.A.node_num
            for OUT_node in el.nodes.B
                j = OUT_node.node_num
                add_edge!(g,i,j)
                edgelabel_mat[i, j] = edgelabel_dict[(i, j)] = el.name
                verbose ? println("Node added :  $(i) to $(j)  <=> $(el.name) ") : nothing
            end
        end
    end
    return g,nms, edgelabel_mat, edgelabel_dict

end

# function network2graph2(cycle_network::component_network; verbose = false, component_dict = Dict())
#     #   Input component Dict to directly add to exisiting Dict
#     num_elem    = elem_count(cycle_network) # number of elements
#     num_vert    = node_count(cycle_network) # number of nodes (vertices)
#     g           = DiGraph(num_elem)
#     nms =[]
#     edgelabel_dict = Dict()
#     edgelabel_mat = Array{String}(undef, num_vert, num_vert)
#     comp_offset = length(component_dict)
#     for (idx1,el1) in enumerate(cycle_network.elements)
#         component_dict[idx1+comp_offset] = el1
#         # add_vertex!(g)
#         push!(nms,el1.name)
#         inn, Pin, Tin, work_fluid = inlet_nodes(el1.nodes)
#         for (idxa,node_idx) in enumerate(inn)
#             press = Pin[idxa]
#             temp = Tin[idxa]
#             edg_str = @sprintf "%s \n P = %.2f \n T = %.2f" String(work_fluid) press temp-273.15
#             for (idx2,el2) in  enumerate(cycle_network.elements)
#                 onn,p2,t2 = outlet_nodes(el2.nodes)
#                 for src_node in onn
#                     if src_node == node_idx
#                         add_edge!(g,idx2,idx1)
#                         edgelabel_mat[idx2, idx1] = edgelabel_dict[(idx2, idx1)] = edg_str
#                         verbose ? println("EDGE added :  $(el2.name) to $(el1.name)  <=> $(src_node) ") : nothing
#                     end
#                 end
#             end
#         end
#     end
#     return g,nms, edgelabel_mat, edgelabel_dict, component_dict
# end

function component_info(cn::component_network)
    @printf "CYCLE OBJECT COMPONENT INFO => %s\n" cn.name
    for el in cn.elements
        component_info(el)
    end
end

"""
    network2graph(sys::system_network; verbose = false)
"""
function network2graph2(sys::system_network; do_plot = false, verbose = false)
    println("waddup")
    @info "Running network2gragh2 of system network"
    gt = [network2graph2(x) for x in sys.cycles]
    gdawg = gt[1][1]
    name_list = gt[1][2]
    edgelabel_mat = gt[1][3]
    Edict = gt[1][4]
    Cdict = gt[1][5]

    @info "Running network2gragh2 of system network"
    
    for i = 2:length(sys.cycles)
        gdawg =  blockdiag(gdawg,gt[i][1])
        vcat_dict!(Edict,gt[i][4])
        name_list = [name_list; gt[i][2]]
        vcat_dic!(Cdict, gt[i][5])
    end

    # graph, Edge dict, name list, vcat dic created
    for hx in sys.heat_exchangers
        hot_element     = hx.hot_stream
        cold_element    = hx.cold_stream
    
        #initializing src trg
        src = -1
        trg = -1
        Qmag = hot_element.Q
        hot_element.Q > 0 ? dir = 1 : dir = -1
    
        #first find by name
        hname = hx.hot_stream.name  # strings
        bname = hx.cold_stream.name
    
        hot_opt = findall(x -> x == hname,name_list)   # node indexes
    
        # if more than 1 index, i.e. hot_elements have the same name
        # check the comp_dict to find matches
        if length(hot_opt) > 1
            for opt in hot_opt
                actual_element = Cdict[opt]
                if get_inlet_temp(actual_element) == get_inlet_temp(hot_element)
                    # actual element found
                    # println("Found element $(get_inlet_temp(actual_element)) and $(get_inlet_temp(hot_element))")
                    src=opt
                end
            end
        else
            src = hot_opt[1]
        end
        cld_opt = findall(x-> x == bname,name_list)
        if length(cld_opt) > 1
            for opt in cld_opt
                actual_element = Cdict[opt]
                if get_inlet_temp(actual_element) == get_inlet_temp(cold_element)
                    # actual element found
                    # println("Found element $(get_inlet_temp(actual_element)) and $(get_inlet_temp(cold_element))")
                    src=opt
                end
            end
        else
            trg = cld_opt[1]
        end
    
        if dir == -1
            tmp = trg
            trg = src
            src = tmp
            Qmag = -Qmag
        end
        add_edge!(gdawg,src,trg)
        Edict[(src,trg)] = @sprintf "Heat Exchanger\n Q = %.2f MW" (Qmag/10^6)
    end
    
    name_list
    qc_idx = findall(x -> x=="Heat_Sink",name_list)
    qh_idx = findall(x -> x=="Fusion_Heat",name_list)
    mflow_idx = findall(x -> x=="M-FLOW SRC",name_list)

    name_list[mflow_idx] .= "ṁ"
    push!(name_list,"Cold Utility")
    push!(name_list,"Fusion Core")


    add_vertex!(gdawg)
    for id in qc_idx
        add_edge!(gdawg,id,nv(gdawg))
        Edict[(id,nv(gdawg))] = "Qc = $(round(Cdict[id].Q/10^6))"
    end

    add_vertex!(gdawg)
    for id in qh_idx
        add_edge!(gdawg,nv(gdawg),id)
        Edict[(nv(gdawg),id)] = "Qfus = $(round(Cdict[id].Q/10^6))"
    end

    if do_plot == true
        plot(gdawg)
    end
    #     positions = Spring(Ptype = Float64, C =3.0, iterations=100)(adjacency_matrix(gdawg))
    #     length(positions)
    #     x = zeros(length(positions),1)
    #     y = zeros(length(positions),1)
    #     for (idx,pos) in enumerate(positions)
    #         x[idx]=pos[1] * 1
    #         y[idx]=pos[2] * 1
    #     end
    #     fig2 = graphplot(gdawg, names = name_list,
    #         curves=false,
    #         nodeshape   =  :rect,
    #         edgelabel   = Edict,
    #         nodecolor = :lightblue,
    #         nodesize    = 0.1,
    #         edge_label_box = true,
    #         axis_buffer = 0.0)

    #         display(fig2)
    #         return gdawg, Edict, Cdict, name_list, fig2
    # end

    return gdawg, Edict, Cdict, name_list

end


function graph2sankey()
    #                  0             1          2          3                   4               5               6           7
    Energy_Sys = ["Fusion Core", "Divertor" , "blanket" , "breeder" ,"Intermediate_loop", "Thermal Cycle", "Cold Utility", "Electric"];
    
    Edict = Dict()  #   Dict([source,targetr]) => value
    FUSION_IDX      = 0
    INTER_IDX       = 4
    CYCLE_IDX       = 5
    ELECTRIC_IDX    = 7
    UTILITY_IDX     = 6
    
    # imediate loops
    for i = 1:3
        Edict[(FUSION_IDX,i)]     = df.Qh[i]/10^6
        Edict[(ELECTRIC_IDX,i)]   = df.Win[i]/10^6
        Edict[(i,INTER_IDX)]      = -df.Qtx[i]/10^6
        Edict[(i,UTILITY_IDX)]    = df.Qc[i]/10^6
    end
    
    Edict[(INTER_IDX,UTILITY_IDX)] = df.Qc[INTER_IDX] / 10^6
    
    # primary cycles
    Edict[(INTER_IDX,CYCLE_IDX)]    = df.Qtx[CYCLE_IDX] / 10^6
    Edict[(CYCLE_IDX,ELECTRIC_IDX)] = df.Wout[CYCLE_IDX] / 10^6
    Edict[(CYCLE_IDX,UTILITY_IDX)]  = df.Qc[CYCLE_IDX] / 10^6
    Edict[(CYCLE_IDX,CYCLE_IDX)]    = df.Win[CYCLE_IDX] /10^6
     
    
    dict_k = keys(Edict)
    source_t = [keyval[1] for keyval ∈ dict_k]
    target_t = [keyval[2] for keyval ∈ dict_k]
    vals_t    = [v for v in values(Edict)]

    p2 = PlotlyJS.plot(sankey(
        node = attr(
          pad = 15,
          thickness = 20,
          label = Energy_Sys,
        ),
        link = attr(
          source    = source_t, # indices correspond to labels, eg A1, A2, A1, B1, ...
          target    = target_t,
          value     = vals_t,
          label     = vals_t,
        )), 
      Layout(title_text="Basic Sankey Diagram", font_size=10)
    )
    display(p2)
end
# ============================================================================#
#                                                                             #
#                             Saved cases
#
# =============================================================================

function match_compression_ratio(Pmin,Nc,Nt;    ic_dp = 0.2, rh_dp = 0.2,   regen_dp = 0.2,    heating_dp = 0.2,   cooling_dp = 0.2)
    Pout(pin,ncomp,rp,dp) = pin*rp^ncomp - dp * sum(rp .^ (1:(ncomp-1)))
    Pout_comp(rpc_x) = Pout(Pmin,Nc,rpc_x,ic_dp)
    ΔP_1 = regen_dp + heating_dp;
    ΔP_2 = regen_dp + cooling_dp;

    fzer(rp1,rp2) = Pout(Pout_comp(rp1)-ΔP_1, Nt, 1/rp2, rh_dp) - ΔP_2 - Pmin
    if Nc > Nt
        rpt = 3.75;
        f_c(rp) = fzer(rp,rpt)
        rpc = find_zero(f_c,3)
    else 
        rpc = 3.75;
        f_t(rp) = fzer(rpc,rp)
        rpt = find_zero(f_t,3)
    end
    return rpt, rpc
end

function brayton_network(Tmax_guess; Tmin = 300.0, Nc=1, Nt=1, Nhx = 3, regen=true, Pmin = 50.0,cycle_fluid = :helium)
    params = component_network_params(Tmin,Pmin,Tmax_guess,cycle_fluid)
    # flow source and initial conditinos
    addl_states = 2;
    states_in_between_comp_and_turb = 2 + Nhx   # compout, regen_out, 
     
    totStates   = 1+ 2*Nc + 2*Nt+ 1 + Nhx
    # totComp     = Nt+Nc+numIC+numRH+addl_comp+2 # +2 for primaryheat and primary cooling

    heating_dp = 0.2*Nhx;
    cooling_dp = 0.2;

    rpt,rpc = match_compression_ratio(Pmin,Nc,Nt)

    # Initializing
    nodes = init_gas_nodes(totStates)

    component_vector = Vector{model_component}(undef,1)
    component_vector[1] = mass_flow_source(nodes[1:2];params = mass_flow_source_paremeter(1.0,cycle_fluid))
    
    compIdx = 1:(2*Nc-1);           # component indexes 1 - last compressor
    idcount = repmat(1:Nc,Nc);      # for labels
    for i in compIdx
        idx = i+1
        node_idx = [idx,idx+1]
        if iseven(idx)
            cparams = comp_params(;rpin = rpc)
            push!(component_vector,gas_compressor(nodes[node_idx];params = cparams, name = "compressor_$(idcount[idx])"))
        else
            push!(component_vector, intercooler(nodes[node_idx],Tmin; name_string =  "intercool_$(idcount[idx])"))
        end
    end

    idx += 1
    node_idx = [idx,idx+1]
    push!(component_vector, init_regen_cold(nodes[node_idx]))

    regen_cold_idx = length(component_vector)

    for i =1:Nhx
        idx += 1
        node_idx = [idx,idx+1]
        push!(component_vector, heat_exchanger_heater(nodes[node_idx]; name_string = "Heat_EX_$(i)"))
    end
    hx_idx = (regen_cold_idx+1:length(component_vector))

    turbIdx = 1:(2*Nt-1);           # component indexes 1 - last compressor
    idcount = repmat(1:Nt,Nt);      # for labels
    for i =1:(2*Nt-1)
        idx += 1
        node_idx = [idx,idx+1]
        if isodd(i)
            tparams = turb_params(;rpin = rpt)
            push!(component_vector,gas_turbine(nodes[node_idx];params = tparams, name = "turbine_$(idcount[i])"))
        else
            push!(component_vector, simple_reheater(nodes[node_idx],Tmin; name_string =  "reheat_$(idcount[i])"))
        end
    end

    
    idx += 1
    node_idx = [idx,idx+1]
    push!(component_vector, init_regen_hot(nodes[node_idx]))
    regen_hot_idx = length(component_vector)
    idx += 1
    node_idx = [idx,1]
    push!(component_vector, heat_sink(nodes[node_idx],Tmin,Pmin))

    return component_network(params,nodes,component_vector,1.0,"BRAYTON CYCLE"), [regen_cold_idx, regen_hot_idx], hx_idx

end

function cooling_network(para::component_network_params; Qin_s = 0.0, verbose = false, dynamic = false)
    numNodes = 5
    fluid_type          = gastype(para.working_fluid)
    component_vector    = Vector{model_component}(undef,1)
    if fluid_type == :liq
        nodes                   = init_liq_nodes(numNodes,para.Pmin)
        component_vector[1]     = mass_flow_source(nodes[1:2];  params = mass_flow_source_paremeter(1.0,  para.working_fluid))
        idx = [2,3]
        nxt_nodes = nodes[2:3]
        fix_pressure(nodes[1],para.Pmin)
        fix_pressure(nodes[2],para.Pmin)
        push!(component_vector, init_simple_pump(nxt_nodes)); 
        nxt_nodes = nodes[3:4]
        push!(component_vector, ideal_heat_source(nxt_nodes; Qin = 0.0, ΔP = 6.0, name_string = "Fusion Heat"));

        nxt_nodes =  nodes[4:5]
        push!(component_vector, heat_exchanger_heater(nxt_nodes; Qin = 0.0, ΔP = 6.0, name_string = "To_cycle"))
        idx = [5,1]
        nxt_nodes = nodes[idx]
        if dynamic == false
            push!(component_vector, heat_sink(nxt_nodes, para.Tmin, para.Pmin; hmode = :Tset, name_string = "Heat_Sink"))
        else
            push!(component_vector, heat_sink(nxt_nodes, para.Tmin, para.Pmin; hmode = :cont, name_string = "Heat_Sink"))
        end
    else
        nodes = init_gas_nodes(numNodes,para.Pmin)
        fix_pressure(nodes[1],para.Pmin)
        fix_pressure(nodes[2],para.Pmin)
        rpc = (para.Pmin+3.0) / para.Pmin
        component_vector[1] = mass_flow_source(nodes[1:2];params = mass_flow_source_paremeter(1.0,  para.working_fluid))
        nxt_nodes = nodes[2:3]
        push!(component_vector, gas_compressor(nxt_nodes; params = comp_params(;rpin = rpc), name = "Circulator"))
        nxt_nodes = nodes[3:4]
        push!(component_vector, ideal_heat_source(nxt_nodes; Qin =  0.0, name_string = "Fusion Heat"))
        nxt_nodes =  nodes[4:5]
        push!(component_vector, heat_exchanger_heater(nxt_nodes; Qin = 0.0, name_string = "To_Cycle"))
        idx = [5,1]
        nxt_nodes = nodes[idx]

        if dynamic == false
            push!(component_vector, heat_sink(nxt_nodes, para.Tmin, para.Pmin; hmode = :Tset, name_string = "Heat_Sink"))
        else
            push!(component_vector, heat_sink(nxt_nodes, para.Tmin, para.Pmin; hmode = :cont, name_string = "Heat_Sink"))
        end
    end
    set_Q(component_vector[3],Qin_s; mode = :in)
    for n in nodes
        overide_set_temperature(n,para.Tmin)        # initalizing tempreatures
    end
    return component_network(para,nodes,component_vector)
end

function extraction_loop(para::component_network_params; dynamic = false)
    numNodes = 8
    component_vector = Vector{model_component}(undef,1)

    nodes = init_gas_nodes(numNodes,para.Pmin)
    fix_pressure(nodes[1],para.Pmin)
    fix_pressure(nodes[2],para.Pmin)
    rpc = (para.Pmin+3.0) / para.Pmin
    component_vector[1] = mass_flow_source(nodes[1:2];params = mass_flow_source_paremeter(1.0,  para.working_fluid))
    nxt_nodes = nodes[2:3]
    push!(component_vector, gas_compressor(nxt_nodes; params = comp_params(;rpin = rpc), name = "Circulator"))
    nxt_nodes = nodes[3:4]
    push!(component_vector, heat_exchanger_heater(nxt_nodes; Qin = 0.0, name_string = "HX1"))
    nxt_nodes =  nodes[4:5]
    push!(component_vector, heat_exchanger_heater(nxt_nodes; Qin = 0.0, name_string = "HX2"))
    nxt_nodes =  nodes[5:6]
    push!(component_vector, heat_exchanger_heater(nxt_nodes; Qin = 0.0, name_string = "HX3"))
    nxt_nodes =  nodes[6:7]
    push!(component_vector, heat_exchanger_heater(nxt_nodes; Qin = 0.0, name_string = "Boiler_hx"))
    idx = [7,1]

    if dynamic == false
        push!(component_vector, heat_sink(nodes[idx], para.Tmin, para.Pmin; name_string = "Heat_Sink"))
    else
        push!(component_vector, heat_sink(nxt_nodes, para.Tmin, para.Pmin; hmode = :cont, name_string = "Heat_Sink"))
    end

    for n in nodes
        overide_set_temperature(n,para.Tmin)        # initalizing tempreatures
    end
    return component_network(para,nodes,component_vector)
end

function default_coolant_networks()
    breeder_network     = cooling_network(breeder_network_params();Qin_s =  740.0e6)
    divertor_network    = cooling_network(divertor_network_params();Qin_s =  140.0e6)
    blanket_network     = cooling_network(blanket_network_params();Qin_s =  100.0e6)
    breeder_network.name =  "breeder_network" 
    divertor_network.name = "divertor_network" 
    blanket_network.name =  "blanket_network" 
    return breeder_network, divertor_network,blanket_network
end

function brayton_system()
    # initialize cycle
    cycle_network, reg, hx = brayton_network(500; Nc = 2,Pmin = 15.0, Nt = 1);
    set_initial_pressure(cycle_network,15.0)
    change_and_update_mass_flow(cycle_network,460.0)

    # create regen heat exchanger (part of cycle)
    regen   = heat_exchanger(cycle_network.elements[reg]; name = "REGEN")

    # Setup cycle
    advance(regen)          # initializes heat exchanger inlet temperatuers
    advance(cycle_network)  # initializes pressures / initial temperatures
    process(regen)          # calculating
    process(cycle_network)  # calcualating

    # inialize coolant circuits
    breeder_circuit,  divertor_circuit,  blanket_circuit = default_coolant_networks()
    @show hx
    # Self explanitory
    ensure_pressure_match(breeder_circuit)
    ensure_pressure_match(divertor_circuit)
    ensure_pressure_match(blanket_circuit)

    # set mass flow rates
    change_and_update_mass_flow(breeder_circuit,17723.22)
    change_and_update_mass_flow(divertor_circuit,760.0)
    change_and_update_mass_flow(blanket_circuit,570.0)

    hx_div = heat_exchanger(cycle_network.elements[hx[1]],divertor_circuit.elements[4]; name = "DIVERTOR HX")
    hx_blk = heat_exchanger(cycle_network.elements[hx[2]],blanket_circuit.elements[4]; name = "BLANKET HX")
    hx_brd = heat_exchanger(cycle_network.elements[hx[3]],breeder_circuit.elements[4]; name = "BREEDER HX")

    cycle_vec   = [divertor_circuit,    blanket_circuit, breeder_circuit, cycle_network]
    hx_vec      = [regen,hx_div, hx_blk, hx_brd]

    Edict = Dict() 
    return system_network(cycle_vec, hx_vec, "Rankine with Intermediate loops", Edict)
end

function simple_rankine(;Pmax,Pmid,Pmin,TmaxK)
    nodes = init_water_nodes(8)
    component_vector = Vector{model_component}(undef,7)

    mflow_nodes     = siso_nodes(nodes[1:2])
    turbine_nodes   = simo_nodes(nodes[2],nodes[3:4])
    condensor_nodes = [nodes[4],nodes[5]]
    pump1_nodes     = siso_nodes(nodes[5:6])
    ofw_nodes       = miso_nodes([nodes[3],nodes[6]],nodes[7])
    pump2_nodes     = siso_nodes(nodes[7:8])

    boiler_nodes = [nodes[8],nodes[1]]


    set_temperature(nodes[2],   TmaxK)
    set_temperature(nodes[1],  TmaxK)
    # refMax = initSuperHeatedWater(TmaxK,Pmax)
    # set_h(nodes[1],refMax.h)

    yflow = [nodes[3]]
    xflow = nodes[4:6]

    component_vector[1] = mass_flow_source(nodes = mflow_nodes; params = mass_flow_source_paremeter(working_fluid = :water, mass_flow = 100.0))
    component_vector[5] = h20_simo_turbine(nodes = turbine_nodes; params = h20_turbine_parameters(1.0))
    component_vector[6] = simple_condensor(condensor_nodes)

    component_vector[3] = h20_pump(nodes = pump1_nodes; params = h20_pump_parameters(1.0), name = "Pump1")
    component_vector[4] = h20_pump(nodes = pump2_nodes; params = h20_pump_parameters(1.0), name = "Pump2")

    component_vector[7] = h20_open_feedwater_heater(nodes = ofw_nodes, xflow = xflow, yflow = yflow)

    component_vector[2] = heat_exchanger_heater(boiler_nodes; name_string = "boiler_qin", fixed_T = true, To = TmaxK)

    #initializing
    initialize(component_vector[3], Pmin,Pmid)
    initialize(component_vector[4], Pmid,Pmax)
    initialize(component_vector[5], Pmax,[Pmid,Pmin])

    para = component_network_params(Tmin = 273, Pmin = Pmin, Tmax = TmaxK, working_fluid = :water)
    cm = component_network(para,nodes,component_vector,1.0,"Rankine")

    advance(cm)
    process(cm)
    advance(cm)
    process(cm)
    return cm
end

function rankine_with_loops()
    rankine = simple_rankine(; Pmax = 80.0, Pmid = 12.0, Pmin = 0.1, TmaxK = 600+273.15)
    breeder_circuit    =    cooling_network(breeder_network_params2(); Qin_s =  740.0e6, dynamic = false)
    divertor_circuit    =   cooling_network(divertor_network_params2();Qin_s =  140.0e6, dynamic = false)
    blanket_circuit    =    cooling_network(blanket_network_params2(); Qin_s =  100.0e6, dynamic = false)

    breeder_circuit.name =  "breeder_network" 
    divertor_circuit.name = "divertor_network" 
    blanket_circuit.name =  "blanket_network" 

    ensure_pressure_match(breeder_circuit)
    ensure_pressure_match(divertor_circuit)
    ensure_pressure_match(blanket_circuit)

    change_and_update_mass_flow(breeder_circuit, 20000.0)
    change_and_update_mass_flow(divertor_circuit,1000.0)
    change_and_update_mass_flow(blanket_circuit,1000.0)

    extraction_network = extraction_loop(extraction_network_params2(), dynamic = false)
    extraction_network.name = "Intermediate_loop"

    ensure_pressure_match(extraction_network)
    change_and_update_mass_flow(extraction_network,1000.0)
    
    change_and_update_mass_flow(rankine,300.0)
    #hx in 3 4 5
    hx = [3,4,5]

    # fix_outlet_temp(breeder_circuit.elements[5],breeder_circuit.params.Tmin)

    divertor_heat_exchanger = heat_exchanger(extraction_network.elements[hx[1]],divertor_circuit.elements[4]; name = "div_hx")
    blanket_heat_exchanger  = heat_exchanger(extraction_network.elements[hx[2]],blanket_circuit.elements[4]; name = "blk_hx")
    breeder_heat_exchanger  = heat_exchanger(extraction_network.elements[hx[3]],breeder_circuit.elements[4]; name = "breed_hx")
    boiler_heat_exchanger   = heat_exchanger(rankine.elements[2],   extraction_network.elements[6]; name = "boiler_hx")

    cycle_vec   = [divertor_circuit,    blanket_circuit, breeder_circuit, extraction_network, rankine]
    hx_vec      = [divertor_heat_exchanger, blanket_heat_exchanger, breeder_heat_exchanger, boiler_heat_exchanger]

    eval_order = Dict()
    eval_order[1] = divertor_circuit
    eval_order[2] = divertor_heat_exchanger
    eval_order[3] = extraction_network
    eval_order[4] = blanket_circuit
    eval_order[5] = blanket_heat_exchanger
    eval_order[6] = extraction_network
    eval_order[7] = breeder_circuit
    eval_order[8] = breeder_heat_exchanger
    eval_order[9] = extraction_network
    eval_order[10] = boiler_heat_exchanger
    eval_order[11] = rankine

    return system_network(cycle_vec, hx_vec, "Rankine with Intermediate loops", eval_order)
end

function run_track_temp(cyc::component_network,niter; do_plot = false, min_conv = .01, min_iter = 3, verbose = false,)
    TEMPS = get_multi_temp(cyc.nodes)
    dt = []
    for i = 1:niter
        advance(cyc)
        process(cyc)
        TEMPS = vcat(TEMPS,get_multi_temp(cyc.nodes))
        t_old = sum(TEMPS[i,:])
        t_new = sum(TEMPS[i+1,:])
        push!(dt,abs((t_old - t_new)))
        if (sum(dt[i]) < min_conv) && i > min_iter
            break
        end
    end
    iter_vec = [1:niter+1]

    if do_plot
        plot(dt)
    end

    if verbose
        @printf "\t Network converged at %i itteration, δT = %.3d\n" length(dt) dt[end]
    end
    return TEMPS,dt,iter_vec
end

function run_sys(sys::system_network; NITER::Int64 = 10)
    for cool in sys.cycles
        advance(cool)
        process(cool)
    end
    for i = 1:NITER
        for hx in sys.heat_exchangers
            advance(hx)
            advance(hx.hot_stream)
            advance(hx.cold_stream)
            process(hx)
        end
    end
    for cool in sys.cycles
        advance(cool)
        process(cool)
        for hx in sys.heat_exchangers
            advance(hx)
            advance(hx.hot_stream)
            advance(hx.cold_stream)
        end
    end
end

function eval_sys(sys::system_network; verbose = false)
    Qin_all = []
    Qout_all = []
    Win_all = []
    Wout_all = []
    Qhx_all = []
    name = []
    for (idx,cyc) in enumerate(sys.cycles)
        Qin, Qout, Win, Wout, Qhx =calculate_circuit_net(cyc; verbose = verbose)
        push!(Qin_all, Qin)
        push!(Qout_all, Qout)
        push!(Win_all, Win)
        push!(Wout_all, Wout)
        push!(Qhx_all, sum(Qhx))
        push!(name,cyc.name)
    end
    for energy in [Qin_all, Qout_all, Win_all, Qhx_all, Wout_all]
        sume = sum(energy[1:end])
        push!(energy,sume)
    end
    push!(name,"Total")

    df = DataFrame(name = name,Qh = Qin_all,Qc = Qout_all,Qtx = Qhx_all, Win = Win_all,Wout = Wout_all)
    return df
end

# http://docs.juliaplots.org/latest/graphrecipes/examples/
# ============================================================================#
#                                                                             #
#                              Printing
#
# =============================================================================
# utilities for setup
function ensure_pressure_match(cn::component_network)
    # only use this function for the cooling loops
    Pmin = cn.params.Pmin
    ΔP_0 = 0
    for c in cn.elements
        if typeof(c) != mass_flow_source
            if :default_ΔP ∈ fieldnames(typeof(c.params))
                ΔP_0 += c.params.default_ΔP
            end
        end
    end

    fluid_type  = gastype(cn.params.working_fluid)
    if fluid_type == :gas
        cn.elements[2].params.rp = (Pmin + ΔP_0)/Pmin
    else
        fix_pressure(cn.elements[2].nodes.B,(Pmin + ΔP_0))
    end
end

function change_source_flow(cn::component_network,new_flow::Float64)
    for c in cn.elements
        if typeof(c) == mass_flow_source
            c.mass_flow = new_flow
            cn.mass_flow = new_flow
            break
        end
    end
end

function update_mass_flow(cn::component_network)
    for n in cn.nodes
        n.mass_flow = cn.mass_flow
    end
    for c in cn.elements
        if typeof(c) != mass_flow_source
            c.mass_flow = cn.mass_flow
        end
    end
end

function change_and_update_mass_flow(cn::component_network,mass_flow::Float64)
    change_source_flow(cn,mass_flow)
    update_mass_flow(cn)
end

function set_initial_temperatures(cn::component_network,Tmin::Float64)
    for node in cn.nodes
        set_temperature(node,Tmin)
    end
end

function set_initial_temperatures(cn::component_network)
    Tmin = cn.params.Tmin
    for node in cn.nodes
        set_temperature(node,Tmin)
    end
end

function set_initial_pressure(cn::component_network,Pmin::Float64)
    for node in cn.nodes
        set_pressure(node,Pmin)
    end
end

function set_initial_pressure(cn::component_network)
    Pmin = cn.params.Pmin
    for node in cn.nodes
        set_pressure(node,Pmin)
    end
end

function process(cn::component_network; kw...)

    for comp in cn.elements
        process(comp; kw...)
    end
end

function advance(cn::component_network; kw...)
    for comp in cn.elements
        advance(comp; kw...)
    end
end

#####################################PRINTING AND CHECK
# CHECK FUNCTIONS
function node_enthalpies(node::component_nodes)
    props = calculate_fluid_properties(node)
    @printf "\t %s h = %.2f \n" node.node_num node.h/1e3
end

function node_enthalpies(cm::component_network)
    for n in cm.nodes
        node_enthalpies(n)
    end
end

function show_node_simple(n::node) where node <: component_nodes
    @printf "\t# %i, \t T = %.2f, \t P = %.2f, \t ṁ = %.2f, \t fluid = %s \n"    n.node_num n.T n.P n.mass_flow n.working_fluid
end

function show_node_simple(cn::component_network)
    @printf "NETWORK NODES - %s \n Working fluid: %s\n Tmin = %.2f\n Pmin: %.2f \n" cn.name cn.params.working_fluid cn.params.Tmin cn.params.Pmin
    for n in cn.nodes
        show_node_simple(n)
    end
end

function show_component(c::comp) where comp<:fluid_component
    str = String(nameof(typeof(c)))
    @printf "%s => ΔT = %.2f ,\t ΔP = %.2f ,\t ΔH = %.2f ,\t ΔQ = %.2f,\t ṁ = %.2f  \n" str c.ΔT c.ΔP c.ΔH c.Q c.mass_flow
end

function show_comp_details(cn::component_network)
    for c in cn.elements
        if typeof(c) != mass_flow_source
            show_component(c)
        end
    end
end

function check_energy_balance(cn::component_network; sum_only = false, show_math = false)
    Hsum = 0
    @printf "%s \n" cn.name
    for c in cn.elements
        if typeof(c) != mass_flow_source
            if !sum_only
                if show_math
                    @printf "\t %s  \t ΔH = H_total + H_element = %.2f + %.2f = %.2f \n" c.name Hsum/10^6 c.ΔH/10^6 (Hsum + c.ΔH)/10^6
                else
                    @printf "\t %s\t ΔH = %.2f\n" c.name c.ΔH
                end
            end
            Hsum += c.ΔH
        end
    end
    @printf "\t TOTAL: %s SUM OF ENTHALPY = %.2f\n" cn.name Hsum
end

function check_entropy_balance(cn::component_network; sum_only = false)
    Ssum = 0
    # @printf "%s \n" cn.name
    for c in cn.elements
        if typeof(c) != mass_flow_source
            if !sum_only
                @printf "\t %s\t ΔS = %.2f\n" c.name c.ΔS
            end
            Ssum += c.ΔS
        end
    end
    @printf "\t TOTAL: %s SUM OF ENTROPY = %.2f\n" cn.name Ssum
end

function show_node(n::node) where node <: component_nodes
    ptyp,ttyp = get_node_fix_type(n)
    @printf "Node %i \t T [K] = %.2f %s , P [Bar] = %.2f %s \n" n.node_num (n.T) ttyp (n.P) (ptyp)
end

function show_node(cvec::Vector{<:component_nodes})
    for nod in cvec
        show_node(nod)
    end
end

function show_comp_nodes(c::comp) where comp<:fluid_component
    inP, inT, oP,oT = get_node_fix_type(c.nodes)
    f1str = (inP,inT)
    f2str = (oP,oT)
    str = c.name
    @printf "%s \t \tNode = %s %i => %i %s  \n" str f1str c.nodes.A.node_num c.nodes.B.node_num f2str
end

function show_comp_nodes(cvec::Vector{comp}) where comp<: fluid_component
    for c in cvec
        show_comp_nodes(c)
    end
end

function show_comp_nodes(cn::component_network)
    for c in cn.elements
        show_comp_nodes(c)
    end
end

function component_dp(c::comp) where comp<: fluid_component
    if typeof(c) != gas_regenerator                  
        return (c.nodes.B.P - c.nodes.A.P)
    else
        #evaluating regenerator
        dp1 = (c.hot_nodes.B.P - c.hot_nodes.A.P)
        dp2 =(c.cold_nodes.B.P - c.cold_nodes.A.P)
        return dp1 + dp2
    end
end

function show_node(cn::component_network)
    println(cn.name)
    for nod in cn.nodes
        show_node(nod)
    end
end

function nodal_vs_element_H(cn::component_network)
    Hsum = 0
    @printf "%s \n" cn.name
    for c in cn.elements
        if typeof(c) != mass_flow_source
            ΔT,ΔP,Δu,Δh,Δs = nodal_changes(c.nodes)
            @printf "\t%s\t (Element) ΔH = %.2f,\t(Nodal) ΔH = %.2f\n" c.name c.ΔH/10^6 Δh*c.mass_flow/10^6
            # dh_simple = dh_nodes(c.nodes)
            # @printf "\t%s\t (Element) ΔH = %.2f,\t(Nodal) ΔH = %.2f,\t (Nodal 2) ΔH = %.2f\n" c.name c.ΔH/10^6 Δh*c.mass_flow/10^6 dh_simple*c.mass_flow/10^6
        end
    end
end

function element_io_temp(cn::component_network)
    @printf "%s \n" cn.name
    for c in cn.elements
        @printf "\t %s, \t Tin = %.2f, Tout = %.2f \t Pin = %.2f, Pout = %.2f\n" c.name get_inlet_temp(c.nodes) get_outlet_temp(c.nodes) get_inlet_press(c.nodes) get_outlet_press(c.nodes)
    end
end  

function element_io_any(cn::component_network,io_var::Symbol)
    @printf "%s \n" cn.name
    for c in cn.elements
        inp,outp = calculate_fluid_properties(c.nodes);
        @printf "\t %s, \t %s in = %.2f, %s out = %.2f,\t Difference = %.2f \n" c.name String(io_var) getproperty(inp,io_var) String(io_var)  getproperty(outp,io_var) getproperty(outp,io_var)-getproperty(inp,io_var)
    end
end

function element_QW_simple(cn::component_network; show_math = false)
    for c in cn.elements
        @printf "%s \t\t OUTPUT\t Qin = %.2f \t  Win = %.2f \n" c.name c.Q/10^6 c.W/10^6
    end
end

function element_QW(cn::component_network; show_math = false)
    @printf "%s \n" cn.name
    sumH = 0
    sumW = 0
    sumQ = 0
    Qin = 0; Qout =
    Wout = 0; Win = 0;
    for c in cn.elements
        if show_math == true
            @printf "\t %s \t\t Q_net + W_net = %.3f + %.3f = %.3f \t ΔH = %.3f \n " c.name c.Q/10^6 c.W/10^6 (c.Q/10^6 + c.W/10^6) c.ΔH/10^6
        else
            @printf "\t %s \t\t Q_net + W_net = %e \t ΔH = %e \n " c.name (c.Q + c.W) c.ΔH
        end
        sumQ += c.Q
        sumW += c.W
        sumH += c.ΔH
        c.Q > 0 ? Qin += abs(c.Q) : Qout += abs(c.Q)
        c.W < 0 ? Wout += abs(c.W) : Win += abs(c.W)
    end
    @printf "\t TOTAL \tQ = %.2f \t  W = %.2f \t ΔH = %.2f \n " sumQ/10^6 sumW/10^6 sumH/10^6
    @printf "\t INPUT \tQin = %.2f \t  Win = %.2f \n" Qin/10^6 Win/10^6
    @printf "\t OUTPUT\tQout = %.2f \t  Wout = %.2f \n" Qout/10^6 Wout/10^6
    @printf "\t\t THERMAL EFFECIENCY η = %.2f " sumW/Qin 
end

function calculate_circuit_net(cn::component_network;verbose = false)
    Qin = 0
    Win = 0
    Wout = 0
    Qout = 0
    # Q_hx = ["", 0]
    Q_hx = []
    for c in cn.elements
        if typeof(c) == heat_source         #   heat source
            mode = c.params.source_mode     #   external heating
            if mode == :env 
                c.Q > 0 ? Qin += abs(c.Q) : Qout += abs(c.Q)
            else    # heat exchanger
                # Q_hx =vcat(Q_hx, [c.name,c.Q])
                push!(Q_hx,c.Q)
            end 
        else
            c.W < 0 ? Wout += abs(c.W) : Win += abs(c.W)
        end
    end

    if verbose
        println(cn.name)
        QQ = (Qin/10^6 - Qout/10^6 + sum(Q_hx)/10^6 + Win/10^6 - Wout/10^6)
        @printf "\t Qin = %.2f MW \t Qout = %.2f MW \t Qhx = %.2f \t Win = %.2f MW \t Wout = %.2f MW \n " Qin/10^6 Qout/10^6 sum(Q_hx)/10^6 Win/10^6 Wout/10^6
        @printf "\t\t TOTAL = %.2f \n" QQ
    end
    return Qin, Qout, Win, Wout, Q_hx
end

