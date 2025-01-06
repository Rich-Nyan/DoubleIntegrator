using CairoMakie
using GLMakie
using GLMakie.GLFW
using GLMakie: to_native


function plonk(x1,y1, x2,y2, x3,y3, t, file, margin)
    
    f = Figure()
    loop = 1
    x1_domain = extrema(x1[1:2]) .+ (-1*margin,margin)
    y1_domain = extrema(y1[1:2]) .+ (-1*margin,margin)
    x2_domain = extrema(x2[1:2]) .+ (-1*margin,margin)
    y2_domain = extrema(y2[1:2]) .+ (-1*margin,margin)
    x3_domain = extrema(x3[1:2]) .+ (-1*margin,margin)
    y3_domain = extrema(y3[1:2]) .+ (-1*margin,margin)

    ax = Axis(f[1, 1], xlabel = "x label", ylabel = "y label", title = "Initial Guess Plot", 
            aspect = 1,
            xrectzoom = false, yrectzoom = false,
            xzoomlock = true, yzoomlock = true,
            xpanlock = true, ypanlock = true)

    
    min = minimum([minimum([minimum([x1_domain[1], y1_domain[1]]), minimum([x2_domain[1], y2_domain[1]])]),minimum([x3_domain[1], y3_domain[1]])])
    max = maximum([maximum([maximum([x1_domain[2], y1_domain[2]]), maximum([x2_domain[2], y2_domain[2]])]),maximum([x3_domain[2], y3_domain[2]])])

    CairoMakie.xlims!(ax, (min,max))
    CairoMakie.ylims!(ax, (min,max))

    # Player 1
    CairoMakie.scatter!(
        ax,
        [x1[1]],
        [x1[2]],
        markersize = 25,
    )
    text!(
        ax, 
        "I1",
        position = (x1[1],x1[2]),
        align = (:center,:center),
        fontsize = 8
    )

    CairoMakie.scatter!(
        ax,
        [y1[1]],
        [y1[2]],
        markersize = 25,
    )
    
    text!(
        ax, 
        "F1",
        position = (y1[1],y1[2]),
        align = (:center,:center),
        fontsize = 8
    )

    # Player 2
    CairoMakie.scatter!(
        ax,
        [x2[1]],
        [x2[2]],
        markersize = 25,
    )
   
    text!(
        ax, 
        "I2",
        position = (x2[1],x2[2]),
        align = (:center,:center),
        fontsize = 8
    )

    CairoMakie.scatter!(
        ax,
        [y2[1]],
        [y2[2]],
        markersize = 25,
    )
    
    text!(
        ax, 
        "F2",
        position = (y2[1],y2[2]),
        align = (:center,:center),
        fontsize = 8
    )

    # Player 3
    CairoMakie.scatter!(
        ax,
        [x3[1]],
        [x3[2]],
        markersize = 25,
    )
   
    text!(
        ax, 
        "I3",
        position = (x3[1],x3[2]),
        align = (:center,:center),
        fontsize = 8
    )

    CairoMakie.scatter!(
        ax,
        [y3[1]],
        [y3[2]],
        markersize = 25,
    )
    
    text!(
        ax, 
        "F3",
        position = (y3[1],y3[2]),
        align = (:center,:center),
        fontsize = 8
    )
    # points = Observable(Point2f0[])
    points = Vector{Vector{Float64}}()
    
    if loop == 1
        lastx = x1[1]
        lasty = x1[2]
        jj = x1
        k = y1
    elseif loop == 2
        lastx = x2[1]
        lasty = x2[2]
        jj = x2
        k = y2
    else
        lastx = x3[1]
        lasty = x3[2]
        jj = x3
        k = y3
    end
    add = true

    register_interaction!(ax, :my_interaction) do event::MouseEvent, axis
        initial_trajectory = zeros(10*t+4)
        
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

                if ((k[1]-event.data[1])^2 + (k[2]-event.data[2])^2 < .1)
                    push!(points,k[1:2])
                    println(points)
                    l = x -> (event.data[1] .+ (k[1] - event.data[1]) .* x, event.data[2] .+ (k[2] - event.data[2]) .* x) 
                    lines!(
                        ax,
                        l.(range(0, stop = 1, length = 100)),
                        markersize = 25,
                        color = (:red)
                    )
                    v = Int64(floor(t / length(points)))
                    for i in 1:length(points)
                        if (i == 1)
                            xx = b -> jj[1] .+ (points[1][1] - jj[1]) .* b ./ v
                            yy = b -> jj[2] .+ (points[1][2] - jj[2]) .* b ./ v
                            for j in 1:v
                                initial_trajectory[4*j-3]=xx(j)
                                initial_trajectory[4*j-2]=yy(j)
                                initial_trajectory[4*j-1]=(points[1][1] - jj[1])/v
                                initial_trajectory[4*j]=(points[1][2] - jj[2])/v
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
                            initial_trajectory[4*t + 1] = initial_trajectory[3]-jj[3]
                            initial_trajectory[4*t + 2] = initial_trajectory[4]-jj[4]
                        end
                        if (i > 1)
                            initial_trajectory[4*t+2*v*(i-1) + 1] = initial_trajectory[4*v*(i-1)+3] - initial_trajectory[4*v*(i-1)-1]
                            initial_trajectory[4*t+2*v*(i-1) + 2] = initial_trajectory[4*v*(i-1)+4] - initial_trajectory[4*v*(i-1)]
                        end
                    end
                    points = Vector{Vector{Float64}}()
                    if loop == 1
                        open(file, "w") do io
                            println(io, jj[1], ",", jj[2], ",", jj[3], ",", jj[4])
                            for i in 1:t
                                println(io, initial_trajectory[4*i-3], ",", initial_trajectory[4*i-2], ",", initial_trajectory[4*i-1], ",", initial_trajectory[4*i])
                            end
                        end
                    else
                        open(file, "a") do io
                            println(io, jj[1], ",", jj[2], ",", jj[3], ",", jj[4])
                            for i in 1:t
                                println(io, initial_trajectory[4*i-3], ",", initial_trajectory[4*i-2], ",", initial_trajectory[4*i-1], ",", initial_trajectory[4*i])
                            end
                        end
                    end
                    if loop == 1
                        open("DI_z-guess.txt", "w") do io
                            for i in 1:10*t+4
                                print(io, initial_trajectory[i], " ")
                            end
                        end
                    else
                        open("DI_z-guess.txt", "a") do io
                            for i in 1:10*t+4
                                print(io, initial_trajectory[i], " ")
                            end
                        end
                    end
                    loop += 1
                    if loop == 1
                        lastx = x1[1]
                        lasty = x1[2]
                        jj = x1
                        k = y1
                    elseif loop == 2
                        lastx = x2[1]
                        lasty = x2[2]
                        jj = x2
                        k = y2
                    else
                        lastx = x3[1]
                        lasty = x3[2]
                        jj = x3
                        k = y3
                    end
                    if (loop == 4)
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
plonk([-1,0,0,0],[1,0,0,0],[0,-1,0,0],[0,1,0,0],[1,0,0,0],[-1,0,0,0], 100, "3playertrajectory_3.txt",0.5)
