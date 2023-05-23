module DeadTime
    using Revise, Roots, DataFrames, XSteam, Printf
    using Graphs, GraphPlot , GraphRecipes
    using Plots, LinearAlgebra
    using PlotlyJS

    include("01-fluid_data.jl")
    include("02-nodes.jl")
    include("03-utilities.jl")
    include("04-components.jl")
    include("05-networks.jl")

    
    function init_test()
        sys = rankine_with_loops()
        run_sys(sys; NITER = 5)
        gdawg, Edict, Cdict, name_list = TSM.network2graph2(sys; do_plot = true)
        df = eval_sys(sys)
        Wnet_cycle = df.Wout[5]-df.Win[5]
        @show η_cycle = Wnet_cycle/df.Qtx[5]
        @show η_total = (df.Wout[6]-df.Win[6])/df.Qh[6]
        component_info(sys.cycles[4])
        show_node_simple(sys.cycles[4])
        println("\t η_cyle = $(η_cycle)" )
        println(" η_tot = $(η_total)")
        return sys,gdawg
    end
end