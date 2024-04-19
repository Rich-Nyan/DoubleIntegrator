using DataFrames
using Plots

function plotter(str::AbstractString, index,optimal, radius)
    T = 100
    T += 1
    iterations = ["Initial","NO CONSTRAINTS"]
    
    constrains = [2]
    numofconstraints = length(constrains)
    #iterations = ["Initial", "0.001", "0.01", "0.05", "0.1", "0.25"]
    # iterations = ["Initial", "0.05", "0.01", "0.015", "0.02", "0.025", "0.03"]
    numo = 2
    iter = 2
    reader = readlines(str)
    lines = length(reader)
    s = Vector{Float64}(undef,length(reader)*4)
    for i in 1:lines
        input = map(x -> parse(Float64, x), rsplit(reader[i],","))
        for j in 1:4
            s[4*(i-1)+j]=input[j]
        end
    end
    states = zeros(iter,T,4)
    for i in 1:iter
        for j in 1:T
            for k in 1:4
                states[i,j,k] = s[(i - 1) * (4 * T) + 4 * (j - 1) + k]
            end
        end
    end

    colors = palette(:default)[1:iter]

    # x_domain = extrema(states[:,:,1]) .+ (-0.5,0.5)
    # y_domain = extrema(states[:,:,2]) .+ (-0.5,0.5)
    x1_domain = extrema([states[1,1,1],states[1,T,1]]) .+ (-1,1)
    y1_domain = extrema([states[1,1,2],states[1,T,2]]) .+ (-1,1)
    x2_domain = extrema([states[2,1,1],states[2,T,1]]) .+ (-1,1)
    y2_domain = extrema([states[2,1,2],states[2,T,2]]) .+ (-1,1)
    min = minimum([minimum([x1_domain[1], y1_domain[1]]), minimum([x2_domain[1], y2_domain[1]])])
    max = maximum([maximum([x1_domain[2], y1_domain[2]]), maximum([x2_domain[2], y2_domain[2]])])
    domain = [min,max];
    
    gify = @animate for i=1:1:T*numo
        player = Int(ceil(i/T))
        timeframe = Int(i - (player - 1) * T)

        if (optimal)
            # Setup
            Plots.plot(
            linewidth = 4,
            label = false,
            xlabel = 'x',
            ylabel = 'y',
            title = "Optimal Obstacle Avoidance",
            aspectratio = 1,
            ylimits = domain,
            xlimits = domain,
            )
        else
            Plots.plot(
            linewidth = 4,
            label = false,
            xlabel = 'x',
            ylabel = 'y',
            title = "Infeasible Obstacle Avoidance",
            aspectratio = 1,
            ylimits = domain,
            xlimits = domain,
            )
        end


        # # Trajectory
        # for j in 1:2
        #     Plots.plot!(
        #             [states[j,k,1] for k in 1:1:timeframe],
        #             [states[j,k,2] for k in 1:1:timeframe],
        #             j,
        #             linewidth = 1,
        #             linecolor = colors[j],
        #             label = false,
        #             legend=:outertopright
        #     )
        # end
  


        # Plots.scatter!(
        #     [states[1,1,1]],
        #     [states[1,1,2]],
        #     markersize = 5,
        #     color = :black,
        #     label = false
        # )
        # Plots.scatter!(
        #     [states[1,T,1]],
        #     [states[1,T,2]],
        #     markersize = 5,
        #     color = :black,
        #     label = false
        # )
        # Plots.scatter!(
        #     [states[2,1,1]],
        #     [states[2,1,2]],
        #     markersize = 5,
        #     color = :red,
        #     label = false
        # )
        # Plots.scatter!(
        #     [states[2,T,1]],
        #     [states[2,T,2]],
        #     markersize = 5,
        #     color = :red,
        #     label = false
        # )
        # # Points
        # for j in 1:2
        #     Plots.scatter!(
        #         [states[j,timeframe,1]],
        #         [states[j,timeframe,2]],
        #         j,
        #         markersize = 5,
        #         color = colors[j],
        #         label = false
        #     )
        #     # Plots.annotate!(states[j,timeframe,1], states[j,timeframe,2], Plots.text(string(iterations[j]), 8, :black))
        # end
        
        # for j in 1:numofconstraints
        #     constrain = radius
        #     center = [states[j,timeframe,1],states[j,timeframe,2]]
        #     c = x -> (center[1] .+ constrain .* cos.(x), center[2] .+ constrain .* sin.(x)) 
            
        #     Plots.plot!(
        #         [c(t) for t = range(0, stop=2pi, length=100)],
        #         linewidth = 1,
        #         linecolor = :red,
        #         seriesalpha = 1,
        #         fillcolor = :red,
        #         fill = true,
        #         fillalpha = 0.1,
        #         label = false
        #     )
        # end


        # Trajectory
        for j in 1:numo
            Plots.plot!(
                    [states[2*(player-1)+j,k,1] for k in 1:1:timeframe],
                    [states[2*(player-1)+j,k,2] for k in 1:1:timeframe],
                    j,
                    linewidth = 1,
                    linecolor = colors[j],
                    label = false,
                    legend=:outertopright
            )
        end
  


        Plots.scatter!(
            [states[1,1,1]],
            [states[1,1,2]],
            markersize = 5,
            color = :black,
            label = false
        )
        Plots.scatter!(
            [states[1,T,1]],
            [states[1,T,2]],
            markersize = 5,
            color = :black,
            label = false
        )
        Plots.scatter!(
            [states[2,1,1]],
            [states[2,1,2]],
            markersize = 5,
            color = :red,
            label = false
        )
        Plots.scatter!(
            [states[2,T,1]],
            [states[2,T,2]],
            markersize = 5,
            color = :red,
            label = false
        )
        # Points
        for j in 1:numo
            Plots.scatter!(
                [states[2*(player-1)+j,timeframe,1]],
                [states[2*(player-1)+j,timeframe,2]],
                j,
                markersize = 5,
                color = colors[j],
                label = false
            )
            # Plots.annotate!(states[j,timeframe,1], states[j,timeframe,2], Plots.text(string(iterations[j]), 8, :black))
        end
        
        for j in 1:numo
            constrain = radius
            center = [states[2*(player-1)+j,timeframe,1],states[2*(player-1)+j,timeframe,2]]
            c = x -> (center[1] .+ constrain .* cos.(x), center[2] .+ constrain .* sin.(x)) 
            
            Plots.plot!(
                [c(t) for t = range(0, stop=2pi, length=100)],
                linewidth = 1,
                linecolor = :red,
                seriesalpha = 1,
                fillcolor = :red,
                fill = true,
                fillalpha = 0.1,
                label = false
            )
        end

    end

    #filename = split(str, "/")[end]
    # filename = split(str, ".")[1]
    
    # labels=["$i" for i in iterations]
    # Plots.plot(labels, legend=:outertopright)
    if (optimal)
        gif(gify,"MonteCarlo/trial4gifs/Optimal/test"*string(index)*".gif", fps = 15)
    else
        gif(gify,"MonteCarlo/trial4gifs/Infeasible/test"*string(index)*".gif", fps = 15)
    end
end

obstacle_radius = readlines("MonteCarlo/trial4/radius.txt")
obstacle_lines = length(obstacle_radius)

obstacle_vals = Vector{Float64}(undef,obstacle_lines)
for i in 1:obstacle_lines
    obstacle_vals[i] = parse(Float64,obstacle_radius[i])
end

optimal_content = readlines("MonteCarlo/trial4/optimalvalues.txt")
parsed_optimalvalues = map(x -> parse(Int64, x), rsplit(optimal_content[1], ","))
for i in 1:length(parsed_optimalvalues)
    index = parsed_optimalvalues[i]
    plotter("MonteCarlo/trial4/test"*string(index)*".txt",index,true,obstacle_vals[index])
end

infeasible_content = readlines("MonteCarlo/trial4/localinfeasvalues.txt")
parsed_infeasiblevalues = map(x -> parse(Int64, x), rsplit(infeasible_content[1], ","))
for i in 1:length(parsed_infeasiblevalues)
    index = parsed_infeasiblevalues[i]
    plotter("MonteCarlo/trial4/test"*string(index)*".txt",index,false,obstacle_vals[index])
end


