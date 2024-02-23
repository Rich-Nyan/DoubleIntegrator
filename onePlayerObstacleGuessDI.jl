using CairoMakie
using GLMakie
using GLMakie.GLFW
using GLMakie: to_native


function plonk(x,y, obs, t, file)
    
    f = Figure()

    x_domain = extrema(x[1:2]) .+ (-0.5,0.5)
    y_domain = extrema(y[1:2]) .+ (-0.5,0.5)

    ax = Axis(f[1, 1], xlabel = "x label", ylabel = "y label", title = "Initial Guess Plot", 
            aspect = 1,
            xrectzoom = false, yrectzoom = false,
            xzoomlock = true, yzoomlock = true,
            xpanlock = true, ypanlock = true)

    
    min = minimum([x_domain[1], y_domain[1]])
    max = maximum([x_domain[2], y_domain[2]])

    CairoMakie.xlims!(ax, (min,max))
    CairoMakie.ylims!(ax, (min,max))

     # Circle
    constrain = obs[3]
    c = x -> (obs[1] .+ constrain .* cos.(x), obs[2] .+ constrain .* sin.(x)) 
    lines!(
        ax,  
        c.(range(0, stop = 2pi, length = 100)),  # Passing function values directly
        linewidth = 1,
        alpha = 1
    )

    CairoMakie.scatter!(
        ax,
        [x[1]],
        [x[2]],
        markersize = 25,
    )
    text!(
        ax, 
        "I",
        position = (x[1],x[2]),
        align = (:center,:center),
        fontsize = 8
    )

    CairoMakie.scatter!(
        ax,
        [y[1]],
        [y[2]],
        markersize = 25,
    )
    
    text!(
        ax, 
        "F",
        position = (y[1],y[2]),
        align = (:center,:center),
        fontsize = 8
    )

    # points = Observable(Point2f0[])
    points = Vector{Vector{Float64}}()
    lastx = x[1]
    lasty = x[2]

    initial_trajectory = zeros(10*t+4)
    add = true
    register_interaction!(ax, :my_interaction) do event::MouseEvent, axis
        if add
            if event.type === MouseEventTypes.leftclick
                xval = event.data[1]
                yval = event.data[2]
                push!(points,[xval,yval])
                l = x -> (lastx .+ (event.data[1] - lastx) .* x, lasty .+ (event.data[2] - lasty) .* x) 
                lines!(
                    ax,
                    l.(range(0, stop = 1, length = 100)),
                    markersize = 25,
                    color = (:red)
                )
                lastx = event.data[1]
                lasty = event.data[2]

                if ((y[1]-event.data[1])^2 + (y[2]-event.data[2])^2 < 1)
                    push!(points,y[1:2])
                    println(points)
                    l = x -> (event.data[1] .+ (y[1] - event.data[1]) .* x, event.data[2] .+ (y[2] - event.data[2]) .* x) 
                    lines!(
                        ax,
                        l.(range(0, stop = 1, length = 100)),
                        markersize = 25,
                        color = (:red)
                    )
                    v = Int64(floor(t / length(points)))
                    for i in 1:length(points)
                        if (i == 1)
                            xx = b -> x[1] .+ (points[1][1] - x[1]) .* b ./ v
                            yy = b -> x[2] .+ (points[1][2] - x[2]) .* b ./ v
                            for j in 1:v
                                initial_trajectory[4*j-3]=xx(j)
                                initial_trajectory[4*j-2]=yy(j)
                                initial_trajectory[4*j-1]=(points[1][1] - x[1])/v
                                initial_trajectory[4*j]=(points[1][2] - x[2])/v
                            end
                        end
                        if (i > 1)
                            if (i == length(points))
                                finalv = t - (i-1)*v
                                xx = b -> points[i-1][1] .+ (points[i][1] - points[i-1][1]) .* b ./ finalv
                                yy = b -> points[i-1][2] .+ (points[i][2] - points[i-1][2]) .* b ./ finalv
                                dy = points[i][2] - points[i-1][2]
                                dx = points[i][1] - points[i-1][1]
                                for j in 1:finalv
                                    initial_trajectory[4*v*(i-1) + 4*j-3]=xx(j)
                                    initial_trajectory[4*v*(i-1) + 4*j-2]=yy(j)
                                    initial_trajectory[4*v*(i-1) + 4*j-1]=dx/finalv
                                    initial_trajectory[4*v*(i-1) + 4*j]=dy/finalv
                                end
                            end
                            if (i < length(points))
                                xx = b -> points[i-1][1] .+ (points[i][1] - points[i-1][1]) .* b ./ v
                                yy = b -> points[i-1][2] .+ (points[i][2] - points[i-1][2]) .* b ./ v
                                dy = points[i][2] - points[i-1][2]
                                dx = points[i][1] - points[i-1][1]
                                for j in 1:v
                                    initial_trajectory[4*v*(i-1) + 4*j-3]=xx(j)
                                    initial_trajectory[4*v*(i-1) + 4*j-2]=yy(j)
                                    initial_trajectory[4*v*(i-1) + 4*j-1]=dx/v
                                    initial_trajectory[4*v*(i-1) + 4*j]=dy/v
                                end
                            end
                        end
                    end
                    
                    for i in 1:length(points)
                        if (i==1)
                            initial_trajectory[4*t + 1] = initial_trajectory[3]-x[3]
                            initial_trajectory[4*t + 2] = initial_trajectory[4]-x[4]
                        end
                        if (i > 1)
                            initial_trajectory[4*t+2*v*(i-1) + 1] = initial_trajectory[4*v*(i-1)+3] - initial_trajectory[4*v*(i-1)-1]
                            initial_trajectory[4*t+2*v*(i-1) + 2] = initial_trajectory[4*v*(i-1)+4] - initial_trajectory[4*v*(i-1)]
                        end
                    end
                
                    open(file, "w") do io
                        for i in 1:t
                            println(io, initial_trajectory[4*i-3], ",", initial_trajectory[4*i-2], ",", initial_trajectory[4*i-1], ",", initial_trajectory[4*i])
                        end
                    end

                    open("DI_z-guess.txt", "w") do io
                        for i in 1:10*t+4
                            print(io, initial_trajectory[i], " ")
                        end
                    end
                    text!(
                        ax, 
                        "PRESS Q TO QUIT",
                        position = ((max+min)/2,(max+min/2)),
                        align = (:center,:center),
                        fontsize = 30
                    )
                    add = false

                   
                end
            end
        end
    end

    glfw_window = to_native(display(ax.scene))
    register_interaction!(ax, :kill_interaction) do event::KeysEvent, axis
        if (add == false)
            if ispressed(ax.scene, Keyboard.q)
                GLFW.SetWindowShouldClose(glfw_window, true)
            end
        end
    end

f
end

# Optimizer
plonk([0,0,-1,-1],[4,4,-1,-1], [2,2,1.9], 100, "test3/trajectories.txt")
