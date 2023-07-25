using GeometryBasics

function do_lines_intersect(x1, y1, x2, y2, x3, y3, x4, y4; returnpos = false)
    p1 = Point2(x1, y1);
    p2 = Point2(x2, y2);
    p3 = Point2(x3, y3);
    p4 = Point2(x4, y4);
    l1 = Line(p1,p2)
    l2 = Line(p3,p4)
    if returnpos 
        return GeometryBasics.intersects(l1,l2)
    end
    return GeometryBasics.intersects(l1,l2)[1]
    # d1 = direction(x3, y3, x4, y4, x1, y1)
    # d2 = direction(x3, y3, x4, y4, x2, y2)
    # d3 = direction(x1, y1, x2, y2, x3, y3)
    # d4 = direction(x1, y1, x2, y2, x4, y4)

    # if d1 != d2 && d3 != d4
    #     return true
    # end

    # if d1 == 0 && on_segment(x3, y3, x4, y4, x1, y1)
    #     return true
    # end

    # if d2 == 0 && on_segment(x3, y3, x4, y4, x2, y2)
    #     return true
    # end

    # if d3 == 0 && on_segment(x1, y1, x2, y2, x3, y3)
    #     return true
    # end

    # if d4 == 0 && on_segment(x1, y1, x2, y2, x4, y4)
    #     return true
    # end

    # return false
end

function direction(x1, y1, x2, y2, x3, y3)
    return (x3 - x1) * (y2 - y1) - (x2 - x1) * (y3 - y1)
end

function on_segment(x1, y1, x2, y2, x3, y3)
    return (min(x1, x2) <= x3 <= max(x1, x2) && min(y1, y2) <= y3 <= max(y1, y2))
end

function edge_vertices(edge,layout)
    x1, y1 = layout[src(edge)]
    x2, y2 = layout[dst(edge)]
    return x1,y1,x2,y2
end

function edge_cross_color(graph,layout)
    e_colors = [:black for i=1:ne(graph)]

    for e1 in edges(graph)
        for e2 in edges(graph)
        end
    end


end
#
function edge_crossings(graph, layout; returnpos = false,return_edges = false)
    crossings = 0
    poscros = Point2f[]
    skipped = []
    for edge1 in edges(graph)
        # x1, y1 = layout[src(edge1)]
        # x2, y2 = layout[dst(edge1)]

        x1,y1,x2,y2 = edge_vertices(edge1,layout)

        for edge2 in edges(graph)
            if edge2 == edge1 || edge1.src == edge2.src || edge1.dst == edge2.dst || edge1.dst == edge2.src || edge1.src == edge2.dst
                continue
            end
            x3,y3,x4,y4 = edge_vertices(edge2,layout)
            if returnpos
                retobj = do_lines_intersect(x1, y1, x2, y2, x3, y3, x4, y4; returnpos = returnpos)
                if retobj[1] == true
                    crossings += 1
                    push!(poscros,retobj[2])
                    
                end
            else
                if do_lines_intersect(x1, y1, x2, y2, x3, y3, x4, y4; returnpos = returnpos)
                    crossings += 1
                end
            end
        end
    end
    if returnpos == false
        return crossings ÷ 2  # Each crossing is counted twice
    else
        return crossings ÷ 2 , poscros
    end
end

function minimize_edge_crossings(graph, start_layout; NITER = 100, K = 1.5,cols = :auto, verbose = false, fixednodes = [])
    # Create a graph from the adjacency matrix

    # start_layout = squaregrid(adjacency_matrix(graph); Ptype = Float64, cols = cols)
    # start_layout = Shell(;Ptype =Float32)(graph)
    # Optimize the layout to minimize edge crossings
    crossings = edge_crossings(graph, start_layout)
    best_crossings = crossings
    best_layout = deepcopy(start_layout)

    current_layout  = deepcopy(start_layout)  
    previous_layout = deepcopy(start_layout)

    verbose ? println("Initial crossings = $(best_crossings)") : nothing
    
    for i in 1:NITER
        # Move each vertex randomly
        for edge1 in edges(graph)
            x1,y1,x2,y2 = edge_vertices(edge1,current_layout)


            for edge2 in edges(graph)
                if edge2 == edge1
                    continue
                end
                x1,y1,x2,y2 = edge_vertices(edge1,current_layout)
                x3,y3,x4,y4 = edge_vertices(edge2,layout)
                bool = do_lines_intersect(x1, y1, x2, y2, x3, y3, x4, y4)

                if bool
                    edge1_node = edge1.src;
                    edge2_node = edge2.dst;

                    randBool1 = rand(1)[1]-0.5
                    randBool2 = rand(1)[1]-0.5

                    if randBool1>0
                        edge1_node=edge1.dst
                    end

                    if randBool2>0
                        edge2_node=edge2.src
                    end

                    if edge1_node ∈ fixednodes
                        edge1_node = setdiff([edge1.src,edge1.dst], [edge1_node])[1]
                    end

                    
                    if edge2_node ∈ fixednodes
                        edge2_node = setdiff([edge2.src,edge2.dst], [edge2_node])[1]
                    end


                    stored_var = current_layout[edge1_node]
                    current_layout[edge1_node] = current_layout[edge2_node] + K .* (rand(2) .- 0.5)
                    current_layout[edge2_node] = stored_var + K .* (rand(2) .- 0.5)
                    break
                end
            end
        end

        # Calculate the number of edge crossings
        crossings = edge_crossings(graph, current_layout)
        chag = sum(current_layout-start_layout)
        verbose ? println("     Iter: $(i) , crossings = $(crossings), Δ = $(chag)") : nothing


        # Keep the layout with the fewest edge crossings
        if crossings < best_crossings
            best_crossings = crossings
            best_layout = deepcopy(current_layout)
            previous_layout = deepcopy(current_layout)
        end
        current_layout = deepcopy(previous_layout)
        crossings = edge_crossings(graph, current_layout)
        chag = sum(current_layout-start_layout)
        verbose ? println("     Iter: $(i) , crossings = $(crossings), Δ = $(chag)") : nothing
    end


        # for i in 1:NITER
        #     # Move each vertex randomly
        #     for i = 1:nv(graph)
        #         current_layout[i] += K .* (rand(2) .- 0.5)
        #     end
    
        #     # Calculate the number of edge crossings
        #     crossings = edge_crossings(graph, current_layout)

        #     verbose ? println("     Iter: $(i) , crossings = $(crossings)") : nothing

        #     chag = newlayout-start_layout
        #     println(sum(chag))
    
        #     # Keep the layout with the fewest edge crossings
        #     if crossings < best_crossings
        #         best_crossings = crossings
        #         best_layout = copy(current_layout)
        #     end
        # end
    crossings = edge_crossings(graph, current_layout)
    println("     Iter: $(NITER) , crossings = $(crossings)")
    return best_layout;
end

function minimize_edge_crossings_shell(graph,inipos,nlist; scaling = 1.0, nlist_shuffles = 5, NITER = 50, max_n_perms = 5, max_rotation = pi/4, verbose = false, fixedlevel=[], fixednodes = [])
    # Optimize the layout to minimize edge crossings
    crossings = edge_crossings(graph, inipos)
    best_crossings = crossings
    best_layout = deepcopy(inipos)
    best_nlist = deepcopy(nlist)
    current_layout = best_layout


    for i = 1:NITER
        for lev in nlist
            if lev ∈ fixedlevel
                continue
            end
            lev_nodes = setdiff(lev,fixednodes)
            num_nodes = length(lev_nodes)
            num_perms = min(max_n_perms,num_nodes)

            cpos = best_layout[lev_nodes]
            # @show typeof(cpos)
            rand_pull = Float32.(rand(num_nodes)) 
            rand_angles = rand_pull .* max_rotation .- max_rotation/2
            rand_scale = rand(1)[1]
            # @show typeof(rand_pull)
            rotated_positions = rotate(cpos,rand_angles; scale_factor = rand_scale + 1.0 )
            # @show typeof(rotated_positions)
            rpos = rotated_positions
            templayout = deepcopy(best_layout)
            for nprm in num_perms
                rpos = shuffle(rpos)
                templayout[lev_nodes] = rpos
                # @show del = sum(rpos - cpos)
                crossings = edge_crossings(G,templayout)
                if crossings < best_crossings
                    verbose ? println("ITER $(i) reduced $(best_crossings) to $(crossings)") : nothing

                    best_layout = deepcopy(templayout)
                    best_nlist = cnlist;
                    best_crossings = crossings
                end
            end
        end
    end

    return best_layout
end

function minimize_edge_crossings5(graph, start_layout; NITER = 100, K = 1.5,verbose = true, fixednodes = [])
    # Optimize the layout to minimize edge crossings
    crossings = edge_crossings(graph, start_layout)
    best_crossings = crossings
    best_layout = copy(start_layout)

    for _ in 1:NITER
        # Move each vertex to its average position with connected neighbors
        for v in 1:nv(graph)
            neighbors = all_neighbors(graph, v)
            n = length(neighbors)
            if n == 0 || v ∈ fixednodes
                continue
            end
            verbose ? println(neighbors) : nothing
            verbose ? println(layout[neighbors]) : nothing
            npos = layout[neighbors]
            verbose ? println(npos) : nothing
            
            
            sums = sum(layout[neighbors])
            verbose ? println(sums ./ n) : nothing
            sum_x = sums[1]#sum(layout[neighbors][1])
            sum_y = sums[2]#sum(layout[neighbors][2])
            layout[v] = Point2(sum_x / n, sum_y / n)
        end

        # Calculate the number of edge crossings
        crossings = edge_crossings(graph, layout)

        # Keep the layout with the fewest edge crossings
        if crossings < best_crossings
            best_crossings = crossings
            best_layout = copy(layout)
        end
    end

    return best_layout
end

#   random_swap_nodes(graph,start_layout; NITER = 100, fixednodes = [])
function random_swap_nodes(graph, start_layout; full_check = false, NITER = 100, fixednodes = [], verbose = false)
    # only nmoves nodes not in fixednodes, swaps saved_positions
    #   full check means in checks intersections of entire graph vs just the movable nodes
    adjmat = Array(adjacency_matrix(graph))
    idx_tracker = [1:size(adjmat)[1]...]
    idx_tracker = setdiff(idx_tracker,fixednodes)

    best_full_layout    =  deepcopy(start_layout)
    new_full_layout     = best_full_layout
    best_full_crossings =  edge_crossings(graph, start_layout)
    
    component_adjacencies   = adjmat[idx_tracker,idx_tracker]
    reduced_graph           = DiGraph(component_adjacencies)
    component_pos           = start_layout[idx_tracker]
    best_sys_crossings      = edge_crossings(reduced_graph, component_pos)
    best_sys_crossings      = crossings
    best_sys_layout         = component_pos

    found_improvement = false

    for i = 1:NITER
        new_sys_pos         = shuffle(best_sys_layout)
        sys_crossings       = edge_crossings(reduced_graph, new_sys_pos)

        if full_check == true
            new_full_layout[idx_tracker] .= new_sys_pos
            f_crossings = edge_crossings(graph, new_full_layout)
            
            if f_crossings  < best_full_crossings && (sys_crossings < best_sys_crossings)
                println("ITER $(i) selected reduced $(best_sys_crossings) to $(sys_crossings)")
                println("      full level reduction $(best_full_crossings) to $(f_crossings)")
                best_sys_crossings  = sys_crossings
                best_sys_layout     = new_sys_pos
                best_full_crossings = f_crossings
                best_full_layout    = deepcopy(new_full_layout)
                found_improvement = true
            end
        elseif sys_crossings < best_sys_crossings
            println("ITER $(i) reduced $(best_sys_crossings) to $(sys_crossings)")
            best_sys_crossings  = sys_crossings
            best_sys_layout     = new_sys_pos
            best_full_layout[idx_tracker] .= best_sys_layout
            found_improvement = true
        end
    end

    # if edge_crossings(graph, best_full_layout) 
    if found_improvement
        return best_full_layout
    else
        return start_layout
    end
        
end

function random_swap_nodes2(graph, start_layout, node_list, spot_list; niter1 = 100, niter2 = 100, fixednodes = [], verbose = false)
    # only nmoves nodes not in fixednodes, swaps saved_positions
    #   full check means in checks intersections of entire graph vs just the movable nodes
    #   node list is a list of nodes that can be swapped
    #   node_list = [breeder_idx, wall_idx ...]
    #   spot_list = list of available spots corresponding to each node_list
    best_full_layout    =  deepcopy(start_layout)
    new_full_layout = deepcopy(start_layout)
    best_full_crossings =  edge_crossings(graph, start_layout)
    adjmat = Array(adjacency_matrix(graph))
    idx_tracker = [1:size(adjmat)[1]...]
    idx_tracker = setdiff(idx_tracker,  fixednodes)
    found_improvement = false

    nl2sl = Dict(node_list[i] => spot_list[i] for i =1:length(node_list))
    for nn in niter1
        for nl in node_list
            new_full_layout             = deepcopy(best_full_layout)
            available_spots             = deepcopy(nl2sl[nl])
            for i =1:niter2
                new_sys_pos          = shuffle(available_spots)[1:length(nl)]
                new_full_layout[nl] .= new_sys_pos
                new_sys_crossings    = edge_crossings(G, new_full_layout)

                if new_sys_crossings < best_full_crossings
                    println("ITER $(i) System Intersections Reduced $(best_full_crossings) to $(new_sys_crossings)")
                    best_full_crossings = new_sys_crossings
                    best_full_layout    = deepcopy(new_full_layout)
                    found_improvement = true
                end
            end

        end
    end

    # if edge_crossings(graph, best_full_layout) 
    if found_improvement
        return best_full_layout
    else
        return start_layout
    end
        
end

get_rad(x::Point{2, Float32}) = sqrt(x[1]^2 + x[2]^2)
get_rad(vec::Vector{Point{2, Float32}}) = [get_rad(x) for x in vec]

function rotate(c::T,ang::Real; scale_factor = 1.0, translation = Point2(0.0,0.0)) where T <: Point
    # ONLY USE THIS IF THE COORDS ARE CENTERED AT THE origin
    x = c[1] * scale_factor
    y = c[2] * scale_factor
    sx = x .* cos.(ang) - y .* sin.(ang)
    sy = x .* sin.(ang) + y .* cos.(ang)
    nx = sx .+ translation[1];
    ny = sy .+ translation[2];
    return Point2f(Float32.(nx),Float32.(ny))
end

function rotate(c::Vector{Point2{Float32}},ang::Vector{Float64}; scale_factor = 1.0, translation = Point2(0.0,0.0))
    return [rotate(c[i],ang[i]; scale_factor = scale_factor, translation = translation) for i = 1:length(c)]
end

function edge_rad(edge,layout)
    x1,y1,x2,y2= edge_vertices(edge,layout)
    return (sqrt(x1^2 + y1^2)+sqrt(x2^2 + y2^2))/2
end


function align_nodes(graph, positions)
    n = nv(graph)
    aligned_positions = deepcopy(positions)
    
    for edge in edges(graph)
        src, dest = edge.src, edge.dst
        dx = aligned_positions[dest][1] - aligned_positions[src][1]
        dy = aligned_positions[dest][2] - aligned_positions[src][2]
        
        if dx != 0 && dy != 0
            if abs(dx) > abs(dy)
                # Align horizontally
                aligned_positions[dest] = Point2f(aligned_positions[dest][1], aligned_positions[src][2])
                # aligned_positions[dest][2] = positions[src][2]
            else
                # # Align vertically
                # aligned_positions[dest][1] = positions[src][1]
                aligned_positions[dest] = Point2f(aligned_positions[src][1], aligned_positions[dest][2])
            end
        end
    end
    
    return aligned_positions
end

function straighten_graph(graph,inipos,to_remove; keep_layout = false)
    finalpos = deepcopy(inipos)
    gg = SimpleGraph()
    
    add_vertices!(gg,nv(graph))
    idx_tracker = [1:nv(graph)...]
    idx_tracker = setdiff(idx_tracker,to_remove)
    
    
    for e in edges(graph)
        add_edge!(gg,e.src,e.dst)
    end

    for uc in reverse(sort(to_remove))
        rem_vertex!(gg,uc)
    end
    
    slayout  = finalpos[idx_tracker]
    if keep_layout == false
        slayout = stress(gg)
    end
    
    nlayout = align_nodes(gg,slayout)
    finalpos[idx_tracker] .= nlayout
    return finalpos
end