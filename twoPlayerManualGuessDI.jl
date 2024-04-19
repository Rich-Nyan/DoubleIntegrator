
function fedInFunction(initial_pose1,final_pose1,initial_pose2,final_pose2, t)
    points1 = Vector{Vector{Float64}}()
    points2 = Vector{Vector{Float64}}()
    if (initial_pose1[2] > initial_pose2[2])
        push!(points1,[initial_pose1[1],initial_pose1[2]+3])
        push!(points2,[initial_pose2[1],initial_pose2[2]-3])
        if (final_pose1[1] > final_pose2[1])
            push!(points1, [final_pose1[1] + 3, initial_pose1[2] + 3])
            push!(points1, [final_pose1[1] + 3, final_pose1[2]])
            push!(points1, [final_pose1[1], final_pose1[2]])

            push!(points2, [final_pose2[1] - 3, initial_pose2[2] - 3])
            push!(points2, [final_pose2[1] - 3, final_pose2[2]])
            push!(points2, [final_pose2[1], final_pose2[2]])
        else
            push!(points1, [final_pose1[1] - 3, initial_pose1[2] + 3])
            push!(points1, [final_pose1[1] - 3, final_pose1[2]])
            push!(points1, [final_pose1[1], final_pose1[2]])

            push!(points2, [final_pose2[1] + 3, initial_pose2[2] - 3])
            push!(points2, [final_pose2[1] + 3, final_pose2[2]])
            push!(points2, [final_pose2[1], final_pose2[2]])
        end
    else
        push!(points1,[initial_pose1[1],initial_pose1[2]-3])
        push!(points2,[initial_pose2[1],initial_pose2[2]+3])
        if (final_pose1[1] > final_pose2[1])
            push!(points1, [final_pose1[1] + 3, initial_pose1[2] - 3])
            push!(points1, [final_pose1[1] + 3, final_pose1[2]])
            push!(points1, [final_pose1[1], final_pose1[2]])

            push!(points2, [final_pose2[1] - 3, initial_pose2[2] + 3])
            push!(points2, [final_pose2[1] - 3, final_pose2[2]])
            push!(points2, [final_pose2[1], final_pose2[2]])
        else
            push!(points1, [final_pose1[1] - 3, initial_pose1[2] - 3])
            push!(points1, [final_pose1[1] - 3, final_pose1[2]])
            push!(points1, [final_pose1[1], final_pose1[2]])

            push!(points2, [final_pose2[1] + 3, initial_pose2[2] + 3])
            push!(points2, [final_pose2[1] + 3, final_pose2[2]])
            push!(points2, [final_pose2[1], final_pose2[2]])
        end
    end

    initial_trajectory = zeros(2*(10*t+4))
    v = Int64(floor(t / length(points1)))
    # Player 1
    for i in 1:length(points1)
        if (i == 1)
            xx = b -> initial_pose1[1] .+ (points1[1][1] - initial_pose1[1]) .* b ./ v
            yy = b -> initial_pose1[2] .+ (points1[1][2] - initial_pose1[2]) .* b ./ v
            for j in 1:v
                initial_trajectory[4*j-3]=xx(j)
                initial_trajectory[4*j-2]=yy(j)
                initial_trajectory[4*j-1]=(points1[1][1] - initial_pose1[1])/v
                initial_trajectory[4*j]=(points1[1][2] - initial_pose1[2])/v
            end
        end
        if (i > 1)
            if (i == length(points1))
                finalv = t - (i-1)*v
                xx = b -> points1[i-1][1] .+ (points1[i][1] - points1[i-1][1]) .* b ./ finalv
                yy = b -> points1[i-1][2] .+ (points1[i][2] - points1[i-1][2]) .* b ./ finalv
                dy = points1[i][2] - points1[i-1][2]
                dx = points1[i][1] - points1[i-1][1]
                for j in 1:finalv
                    initial_trajectory[4*v*(i-1) + 4*j-3]=xx(j)
                    initial_trajectory[4*v*(i-1) + 4*j-2]=yy(j)
                    initial_trajectory[4*v*(i-1) + 4*j-1]=dx/finalv
                    initial_trajectory[4*v*(i-1) + 4*j]=dy/finalv
                end
            end
            if (i < length(points1))
                xx = b -> points1[i-1][1] .+ (points1[i][1] - points1[i-1][1]) .* b ./ v
                yy = b -> points1[i-1][2] .+ (points1[i][2] - points1[i-1][2]) .* b ./ v
                dy = points1[i][2] - points1[i-1][2]
                dx = points1[i][1] - points1[i-1][1]
                for j in 1:v
                    initial_trajectory[4*v*(i-1) + 4*j-3]=xx(j)
                    initial_trajectory[4*v*(i-1) + 4*j-2]=yy(j)
                    initial_trajectory[4*v*(i-1) + 4*j-1]=dx/v
                    initial_trajectory[4*v*(i-1) + 4*j]=dy/v
                end
            end
        end
    end

    v = Int64(floor(t / length(points2)))
    # Player 2
    for i in 1:length(points2)
        if (i == 1)
            xx = b -> initial_pose2[1] .+ (points2[1][1] - initial_pose2[1]) .* b ./ v
            yy = b -> initial_pose2[2] .+ (points2[1][2] - initial_pose2[2]) .* b ./ v
            for j in 1:v
                initial_trajectory[(10*t+4) + 4*j-3]=xx(j)
                initial_trajectory[(10*t+4) + 4*j-2]=yy(j)
                initial_trajectory[(10*t+4) + 4*j-1]=(points2[1][1] - initial_pose2[1])/v
                initial_trajectory[(10*t+4) + 4*j]=(points2[1][2] - initial_pose2[2])/v
            end
        end
        if (i > 1)
            if (i == length(points2))
                finalv = t - (i-1)*v
                xx = b -> points2[i-1][1] .+ (points2[i][1] - points2[i-1][1]) .* b ./ finalv
                yy = b -> points2[i-1][2] .+ (points2[i][2] - points2[i-1][2]) .* b ./ finalv
                dy = points2[i][2] - points2[i-1][2]
                dx = points2[i][1] - points2[i-1][1]
                for j in 1:finalv
                    initial_trajectory[(10*t+4) + 4*v*(i-1) + 4*j-3]=xx(j)
                    initial_trajectory[(10*t+4) + 4*v*(i-1) + 4*j-2]=yy(j)
                    initial_trajectory[(10*t+4) + 4*v*(i-1) + 4*j-1]=dx/finalv
                    initial_trajectory[(10*t+4) + 4*v*(i-1) + 4*j]=dy/finalv
                end
            end
            if (i < length(points2))
                xx = b -> points2[i-1][1] .+ (points2[i][1] - points2[i-1][1]) .* b ./ v
                yy = b -> points2[i-1][2] .+ (points2[i][2] - points2[i-1][2]) .* b ./ v
                dy = points2[i][2] - points2[i-1][2]
                dx = points2[i][1] - points2[i-1][1]
                for j in 1:v
                    initial_trajectory[(10*t+4) + 4*v*(i-1) + 4*j-3]=xx(j)
                    initial_trajectory[(10*t+4) + 4*v*(i-1) + 4*j-2]=yy(j)
                    initial_trajectory[(10*t+4) + 4*v*(i-1) + 4*j-1]=dx/v
                    initial_trajectory[(10*t+4) + 4*v*(i-1) + 4*j]=dy/v
                end
            end
        end
    end

    v = Int64(floor(t / length(points1)))
    # Player 1
    for i in 1:length(points1)
        if (i==1)
            initial_trajectory[4*t + 1] = initial_trajectory[3]-initial_pose1[3]
            initial_trajectory[4*t + 2] = initial_trajectory[4]-initial_pose1[4]
        end
        if (i > 1)
            initial_trajectory[4*t+2*v*(i-1) + 1] = initial_trajectory[4*v*(i-1)+3] - initial_trajectory[4*v*(i-1)-1]
            initial_trajectory[4*t+2*v*(i-1) + 2] = initial_trajectory[4*v*(i-1)+4] - initial_trajectory[4*v*(i-1)]
        end
    end

    v = Int64(floor(t / length(points1))) 
    # Player 2
    for i in 1:length(points2)
        if (i==1)
            initial_trajectory[(10*t+4) + 4*t + 1] = initial_trajectory[(10*t+4) + 3]-initial_pose2[3]
            initial_trajectory[(10*t+4) + 4*t + 2] = initial_trajectory[(10*t+4) + 4]-initial_pose2[4]
        end
        if (i > 1)
            initial_trajectory[(10*t+4) + 4*t+2*v*(i-1) + 1] = initial_trajectory[(10*t+4) + 4*v*(i-1)+3] - initial_trajectory[(10*t+4) + 4*v*(i-1)-1]
            initial_trajectory[(10*t+4) + 4*t+2*v*(i-1) + 2] = initial_trajectory[(10*t+4) + 4*v*(i-1)+4] - initial_trajectory[(10*t+4) + 4*v*(i-1)]
        end
    end

    open("DI_test.txt", "w") do io
        println(io, initial_pose1[1], ",", initial_pose1[2], ",", initial_pose1[3], ",", initial_pose1[4])
        for i in 1:t
            println(io, initial_trajectory[4*i-3], ",", initial_trajectory[4*i-2], ",", initial_trajectory[4*i-1], ",", initial_trajectory[4*i])
        end

        println(io, initial_pose2[1], ",", initial_pose2[2], ",", initial_pose2[3], ",", initial_pose2[4])
        for i in 1:t
            println(io, initial_trajectory[10*t+4+4*i-3], ",", initial_trajectory[10*t+4+4*i-2], ",", initial_trajectory[10*t+4+4*i-1], ",", initial_trajectory[10*t+4+4*i])
        end
    end
    return initial_trajectory
end

fedInFunction([0,0,0,0],[4,0,0,0],[2,0,0,0],[0,3,0,0], 100)