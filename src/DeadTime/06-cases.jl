include("04-components.jl")


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
# t,dt,it = run_track_temp(sys.cycles[5],50)
# p=plot(dt)
