using DataFrames
using Plots

function plotter(str::AbstractString)
    T = 100
    T += 1
    iterations = ["Initial", "0.0001","0.001", "0.01", "0.05", "0.1", "0.25","0.5"]
    #iterations = ["Initial", "0.001", "0.01", "0.05", "0.1", "0.25"]
    # iterations = ["Initial", "0.05", "0.01", "0.015", "0.02", "0.025", "0.03"]
    iter = length(iterations)
    initial_pose = [0,0,-1,-1]
    final_pose = [4,0,-1,-1]
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
    x_domain = extrema([initial_pose[1],final_pose[1]]) .+ (-1,1)
    y_domain = extrema([initial_pose[2],final_pose[2]]) .+ (-1,1)
    domain  = [minimum([x_domain[1],y_domain[1]]),maximum([x_domain[2],y_domain[2]])]
    
    gify = @animate for i=1:1:T
        player = Int(ceil(i/T))
        timeframe = Int(i - (player - 1) * T)


        # Setup
        Plots.plot(
        linewidth = 4,
        label = false,
        xlabel = 'x',
        ylabel = 'y',
        title = "Obstacle Avoidance",
        aspectratio = 1,
        ylimits = domain,
        xlimits = domain,
        )

        # Circle
        constrain = 1.9
        center = (2,0)
        c = x -> (center[1] .+ constrain .* cos.(x), center[2] .+ constrain .* sin.(x)) 
        
        Plots.plot!(
            [c(t) for t = range(0, stop=2pi, length=100)],
            linewidth = 1,
            linecolor = :red,
            seriesalpha = 1,
            # fill = true,
            # fillalpha = 0.1,
            label = false
        )
        # Current Player
        
        # Trajectory
        Plots.plot!(
                [states[1,k,1] for k in 1:1:timeframe],
                [states[1,k,2] for k in 1:1:timeframe],
                linewidth = 1,
                linecolor = colors[1],
                label = string(iterations[1]),
                legend=:outertopright
            )
        for j in 2:length(iterations)
            Plots.plot!(
                [states[j,k,1] for k in 1:1:timeframe],
                [states[j,k,2] for k in 1:1:timeframe],
                linewidth = 1,
                linecolor = colors[j],
                label = "α = " * string(iterations[j]),
                legend=:outertopright
            )
        end

        Plots.scatter!(
            [initial_pose[1]],
            [initial_pose[2]],
            markersize = 5,
            color = :black,
            label = false
        )
        Plots.scatter!(
            [final_pose[1]],
            [final_pose[2]],
            markersize = 5,
            color = :black,
            label = false
        )
        # Points
        for j in 1:length(iterations)
            Plots.scatter!(
                [states[j,timeframe,1]],
                [states[j,timeframe,2]],
                markersize = 5,
                color = colors[j],
                label = false
            )
            # Plots.annotate!(states[j,timeframe,1], states[j,timeframe,2], Plots.text(string(iterations[j]), 8, :black))
        end
        


    end

    filename = split(str, "/")[end]
    filename = split(filename, ".")[1]
    
    # labels=["$i" for i in iterations]
    # Plots.plot(labels, legend=:outertopright)

    gif(gify, "test5.5/" * filename * ".gif", fps = 15)
end

plotter("test5.5/trajectories.txt")

