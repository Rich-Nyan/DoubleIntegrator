using DataFrames
using Plots

function plotter(str::AbstractString)
    T = 100
    T += 1
    iterations = ["Initial","NO CONSTRAINTS"]
    
    constrains = [2]
    numofconstraints = length(constrains)
    #iterations = ["Initial", "0.001", "0.01", "0.05", "0.1", "0.25"]
    # iterations = ["Initial", "0.05", "0.01", "0.015", "0.02", "0.025", "0.03"]
    numo = 2
    iter = 2*numo
    initial_pose1 = [0,0,1,0]
    final_pose1 = [4,4,1,0]
    initial_pose2 = [4,0,1,0]
    final_pose2 = [0,4,1,0]
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
    x1_domain = extrema([initial_pose1[1],final_pose1[1]]) .+ (-1,1)
    y1_domain = extrema([initial_pose1[2],final_pose1[2]]) .+ (-1,1)
    x2_domain = extrema([initial_pose2[1],final_pose2[1]]) .+ (-1,1)
    y2_domain = extrema([initial_pose2[2],final_pose2[2]]) .+ (-1,1)
    min = minimum([minimum([x1_domain[1], y1_domain[1]]), minimum([x2_domain[1], y2_domain[1]])])
    max = maximum([maximum([x1_domain[2], y1_domain[2]]), maximum([x2_domain[2], y2_domain[2]])])
    domain = [min,max];
    
    gify = @animate for i=1:1:T*numo
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
            [initial_pose1[1]],
            [initial_pose1[2]],
            markersize = 5,
            color = :black,
            label = false
        )
        Plots.scatter!(
            [final_pose1[1]],
            [final_pose1[2]],
            markersize = 5,
            color = :black,
            label = false
        )
        Plots.scatter!(
            [initial_pose2[1]],
            [initial_pose2[2]],
            markersize = 5,
            color = :red,
            label = false
        )
        Plots.scatter!(
            [final_pose2[1]],
            [final_pose2[2]],
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
        
        for j in 1:numofconstraints
            constrain = constrains[j]
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
    filename = split(str, ".")[1]
    
    # labels=["$i" for i in iterations]
    # Plots.plot(labels, legend=:outertopright)

    gif(gify, filename * ".gif", fps = 15)
end

plotter("2playertrajectory_1.txt")

