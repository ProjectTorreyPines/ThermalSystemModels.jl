using .ThermalSystem_Models


function init_g()
TSM = ThermalSystem_Models
sys = TSM.rankine_with_loops()
run_sys(sys)
h=TSM.network2graph2(sys)
gdawg, Edict, Cdict, name_list = network2graph2(sys)
return gdawg, Edict, name_list
end


# using PlotlyJS
    # plotlyjs()
    # pp=PlotlyJS.plot(sankey(
    #
    #     node = attr(
    #       pad = 10,
    #       thickness = 15,
    #       line = attr(color = "black", width = 0.5),
    #       standoff = 25,
    #       label = labs
    #     ),
    #     link = attr(
    #       source = src,     # indices correspond to labels, eg A1, A2, A1, B1, ...
    #       target =  dst,
    #       value  = val,
    #       standoff = 25,
    #       label = ["$(vdict[src[i]]) => $(vdict[dst[i]])" for i =1:ne(G)],
    #       showarrow = true,
    #       arrowlen = attr(arrowlen = 5.0))
    #   ),
    #   Layout(title_text="Basic Sankey Diagram", 
    #             showarrow=true,
    #             width=800, height=500,font_size=15))

    # begin
    #     src = [c.src for c in collect(edges(G))]
    #     dst = [c.dst for c in collect(edges(G))]
    #     val = [(sol[epack[(src[i],dst[i])]][end])/10^6 for i = 1:ne(G)]
    
    #     @show mval = minimum(val)
    #     for (i,v) in enumerate(val)
    #         s = src[i]
    #         d = dst[i]
    
    #         str_end = TSMD.split_variable_sym(epack[(s,d)])
    #         if str_end[end] == "Φ"
    #             val[i] = val[i] + abs(mval)
    #         end
    #     end
    #     # val[findall(x-> x==0,val)] .= 0.1
    #     labs = [String(vdict[i]) for i =1:nv(G)]
    #     labs = vcat("",labs) 
    # end
    