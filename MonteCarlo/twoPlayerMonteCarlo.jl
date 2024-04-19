using Ipopt
using Optim
using LinearAlgebra
using ForwardDiff
using Plots
using JuMP
using Random

# Starts at time stamp 0

# x
# [1,5,9,...,4T-3]: x 
# [2,6,10,...,4T-2]: y 
# [3,7,11,...,4T-1]: deltax
# [4,8,12,...,4T]: deltay

# u
# [1,3,5,...,2T-1]: xaccel
# [2,4,6,...,2T]: yaccel

# Dynamics
# jacob : gets the derivative wrt x axis, y axis, deltax, and deltay
# x : holds x,y,deltax,deltay at that specific timestamp (state vector)
# u : xaccel and y accel

# Values
T = 100
dt = 0.1
α = 0.05
# initial_pose1 = [0,0,1,0]
# final_pose1 = [4,4,1,0]
# initial_pose2 = [4,0,1,0]
# final_pose2 = [0,4,1,0]

# obstacle = 2
z_guess = zeros(20 * T + 8)

iterations = [1]
iter = length(iterations)

function objective_function(z)
    s = sum(x -> x^2, z[4*T+1:6*T]) + sum(x -> x^2, z[4*T+1 + 6*T:6*T + 6*T])
    # s = sum(x -> x^2, z[4*T+1:6*T]) + sum(t -> penalty(z, t), 0:T-1)
    return s
end

function distance(x,y)
    return α * dot(x .- y, x .- y)
end

function penalty(z,t)
    x1_t = [z[4*t+1], z[4*t+2]]
    x2_t = [z[10*T+4+4*t+1], z[10*T+4+4*t+2]]
    return -log(distance(x1_t,x2_t) - obstacle^2)
end

# Constraints h(z)
# h(z)
# T vectors of length 4 (x)
# T vectors of length 2 (u)
function dynamic_feasibility(z, model, T,initial_pose1,initial_pose2,final_pose1,final_pose2)

    h_z = zeros(QuadExpr, 8*T+8)

    for a in 1:2
        x = z[1+(6*T)*(a-1):4*T+(6*T)*(a-1)]
        u = z[4*T+1+(6*T)*(a-1):6*T*(a)]

        

        # Initial Pose 
        if a == 1
            h_z[1+(4*T+4)*(a-1)] = (x[1] - (initial_pose1[1] + dt * initial_pose1[3]))
            h_z[2+(4*T+4)*(a-1)] = (x[2] - (initial_pose1[2] + dt * initial_pose1[4]))
            h_z[3+(4*T+4)*(a-1)] = (x[3] - (initial_pose1[3] + dt * u[1]))
            h_z[4+(4*T+4)*(a-1)] = (x[4] - (initial_pose1[4] + dt * u[2]))
        else
            h_z[1+(4*T+4)*(a-1)] = (x[1] - (initial_pose2[1] + dt * initial_pose2[3]))
            h_z[2+(4*T+4)*(a-1)] = (x[2] - (initial_pose2[2] + dt * initial_pose2[4]))
            h_z[3+(4*T+4)*(a-1)] = (x[3] - (initial_pose2[3] + dt * u[1]))
            h_z[4+(4*T+4)*(a-1)] = (x[4] - (initial_pose2[4] + dt * u[2]))
        end

        # Intermediate Poses
        for i in 1:T-1
            h_z[4*i+1+(4*T+4)*(a-1)] = (x[4*i+1] - (x[4*i-3] + dt * x[4*i-1]))
            h_z[4*i+2+(4*T+4)*(a-1)] = (x[4*i+2] - (x[4*i-2] + dt * x[4*i]))
            h_z[4*i+3+(4*T+4)*(a-1)] = (x[4*i+3] - (x[4*i-1] + dt * u[2*i+1]))
            h_z[4*i+4+(4*T+4)*(a-1)] = (x[4*i+4] - (x[4*i] + dt * u[2*i+2]))
        end

        # Final Pose
        if a == 1
            h_z[4*T+1+(4*T+4)*(a-1)] = x[4*T-3] - final_pose1[1]
            h_z[4*T+2+(4*T+4)*(a-1)] = x[4*T-2] - final_pose1[2]
            h_z[4*T+3+(4*T+4)*(a-1)] = x[4*T-1] - final_pose1[3]
            h_z[4*T+4+(4*T+4)*(a-1)] = x[4*T] - final_pose1[4]
        else
            h_z[4*T+1+(4*T+4)*(a-1)] = x[4*T-3] - final_pose2[1]
            h_z[4*T+2+(4*T+4)*(a-1)] = x[4*T-2] - final_pose2[2]
            h_z[4*T+3+(4*T+4)*(a-1)] = x[4*T-1] - final_pose2[3]
            h_z[4*T+4+(4*T+4)*(a-1)] = x[4*T] - final_pose2[4]
        end
    end

    return h_z
end

function lagrangian(z,λ, h_z)
    f_z = objective_function(z)
    return f_z + dot(λ, h_z)
end

function gradient_f(z)
    grad_f = zeros(QuadExpr, 2*6*T)
    for i in 4*T+1:6*T
        grad_f[i] = 2 * z[i]
    end
    for i in 6*T+4*T+1:2*6*T
        grad_f[i] = 2 * z[i]
    end
    return grad_f
end

function gradientObs_f(z,model, α, obstacle)
    grad_f = zeros(QuadExpr, 2*6*T)
    for i in 4*T+1:6*T
        grad_f[i] = 2 * z[i]
    end
    for i in 6*T+4*T+1:2*6*T
        grad_f[i] = 2 * z[i]
    end
    denom = @variable(model,[1:T])
    for i in 1:T        
        @NLconstraint(model, denom[i] == 1/((z[4*i-3]-z[6*T+4*i-3])^2+(z[4*i-2]-z[6*T+4*i-2])^2-obstacle^2))
        grad_f[4*i-3] = α * denom[i] * 2 * (z[4*i-3] - z[6*T+4*i-3])
        grad_f[4*i-2] = α * denom[i] * 2 * (z[4*i-2] - z[6*T+4*i-2])
    end
    return grad_f
end


function gradient_h(z,model)
    grad_h_z = zeros(QuadExpr,2*6*T,2*(4*T+4))
    for i in 1:2*6*T
        for j in 1:2*(4*T+4)
            grad_h_z[i,j] = 0
        end
    end
    for a in 1:2
        for i in 1:4*T
            grad_h_z[i+(a-1)*6*T,i+(a-1)*(4*T+4)] = 1
        end
        for i in 5:4*T
            grad_h_z[i-4+(a-1)*6*T,i+(a-1)*(4*T+4)] = -1
        end
        for i in 4T+1:4T+4
            grad_h_z[i-4+(a-1)*6*T,i+(a-1)*(4*T+4)] = 1
        end

        for i in 1:4*T
            if (i % 4 == 3)
                grad_h_z[4*T+2*div(i-3,4)+1+(a-1)*6*T,i+(a-1)*(4*T+4)] = -dt
                grad_h_z[4*T+2*div(i-3,4)+2+(a-1)*6*T,i+1+(a-1)*(4*T+4)] = -dt
            end
        end

        for i in 5:4*T
            if (i % 4 == 1)
                grad_h_z[i-2+(a-1)*6*T,i+(a-1)*(4*T+4)] = -dt
            end
            if (i % 4 == 2)
                grad_h_z[i-2+(a-1)*6*T,i+(a-1)*(4*T+4)] = -dt
            end
        end
    end

    return grad_h_z
end

#Optimize Function

function optimize(initial_pose1,initial_pose2,final_pose1,final_pose2)
    model = Model(optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 5))

    @variable(model, z[i = 1:12 * T])

    for i in 1:length(z)
        if i <= 6 * T
            set_start_value(z[i], z_guess[i])
        else
            set_start_value(z[i], z_guess[10 * T + 4 + i - 6*T])
        end
    end

    h_z = dynamic_feasibility(z, model, T, initial_pose1,initial_pose2,final_pose1,final_pose2)
    # Objective
    # @objective(model, Min, objective_function(z))

    # Lagrange gradient
    grad_f = gradient_f(z)
    grad_h = gradient_h(z, model)

    @variable(model, λ[i = 1:8*T+8])

    for i in 1:length(λ)
        if i <= 4 * T + 4
            set_start_value(z[i], z_guess[6*T+i])
        else
            set_start_value(z[i], z_guess[16 * T + 4 + i - 4*T-4])
        end
    end
    for k in 1:2*6*T
        λ_quad = @NLexpression(model, sum(λ[i] * grad_h[k, i] for i in 1:8*T+8))
        grad_f_quad = @NLexpression(model, grad_f[k])
        @NLconstraint(model, grad_f_quad - λ_quad == 0)
    end

    for i in 1:8*T+8
        @constraint(model,h_z[i] == 0)
    end
    vals = zeros(20*T+8)
    for i in 1:20*T+8
        vals[i] = values(i)
    end

    # Optimize the model
    optimize!(model)

    return vals

end

# Manual Feed
function fedInFunction(initial_pose1,final_pose1,initial_pose2,final_pose2, t, num)
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
    open("MonteCarlo/trial5/test" * string(num) * ".txt", "w") do io
        println(io, initial_pose1[1], ",", initial_pose1[2], ",", initial_pose1[3], ",", initial_pose1[4])
        for i in 1:T
                println(io, initial_trajectory[4*(i-1)+1], ",", initial_trajectory[4*(i-1)+2], ",", initial_trajectory[4*(i-1)+3], ",", initial_trajectory[4*(i-1)+4])
        end
        println(io, initial_pose2[1], ",", initial_pose2[2], ",", initial_pose2[3], ",", initial_pose2[4])
        for i in 1:T
                println(io, initial_trajectory[(10*T+4)+ 4*(i-1)+1], ",", initial_trajectory[(10*T+4)+ 4*(i-1)+2], ",", initial_trajectory[(10*T+4)+ 4*(i-1)+3], ",", initial_trajectory[(10*T+4)+ 4*(i-1)+4])
        end
    end
    return initial_trajectory
end
function optimizerObstacle(num, initial_pose1,initial_pose2,final_pose1,final_pose2, obstacle, optimalcount, optimalvals, localinfeasvals)
    model = Model(optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 5))

    # arr = optimize(initial_pose1,initial_pose2,final_pose1,final_pose2)
    arr = fedInFunction(initial_pose1,initial_pose2,final_pose1,final_pose2,T, num)
    @variable(model, z[i = 1:12 * T])

    for i in 1:length(z)
        if i <= 6 * T
            set_start_value(z[i], arr[i])
        else
            set_start_value(z[i], arr[10 * T + 4 + i - 6*T])
        end
    end

    
    h_z = dynamic_feasibility(z, model, T, initial_pose1,initial_pose2,final_pose1,final_pose2)
    # Objective
    # @objective(model, Min, objective_function(z))

    # Lagrange gradient
    grad_f = gradientObs_f(z, model, α, obstacle)
    grad_h = gradient_h(z, model)

    @variable(model, λ[i = 1:8*T+8])

    for i in 1:length(λ)
        if i <= 4 * T + 4
            set_start_value(z[i], arr[6*T+i])
        else
            set_start_value(z[i], arr[16 * T + 4 + i - 4*T-4])
        end
    end
    for k in 1:2*6*T
        λ_quad = @NLexpression(model, sum(λ[i] * grad_h[k, i] for i in 1:8*T+8))
        grad_f_quad = @NLexpression(model, grad_f[k])
        @NLconstraint(model, grad_f_quad - λ_quad == 0)
    end

    for i in 1:8*T+8
        @constraint(model,h_z[i] == 0)
    end
    states = zeros(2, T+1, 4)
    for i in 1:iter 
        # Optimize the model
        optimize!(model)

        # Calculate h_z
        # z_guess = value.(z)

        # Player 1
        states[1, 1, 1] = initial_pose1[1]
        states[1, 1, 2] = initial_pose1[2]
        states[1, 1, 3] = initial_pose1[3]
        states[1, 1, 4] = initial_pose1[4]
        for j in 2:T+1
            states[1, j, 1] = value(z[4*(j-1)-3])
            states[1, j, 2] = value(z[4*(j-1)-2])
            states[1, j, 3] = value(z[4*(j-1)-1])
            states[1, j, 4] = value(z[4*(j-1)])
        end
         # Player 2
         states[2, 1, 1] = initial_pose2[1]
         states[2, 1, 2] = initial_pose2[2]
         states[2, 1, 3] = initial_pose2[3]
         states[2, 1, 4] = initial_pose2[4]
         for j in 2:T+1
             states[2, j, 1] = value(z[4*(j-1)-3+6*T])
             states[2, j, 2] = value(z[4*(j-1)-2+6*T])
             states[2, j, 3] = value(z[4*(j-1)-1+6*T])
             states[2, j, 4] = value(z[4*(j-1)+6*T])
         end
         
        
        println("Iteration $i:")
        println("Objective: ", objective_value(model))
        z_value = value.(λ[1:end])
        println("Constraints: ", norm(z_value))
    end
    

    if is_solved_and_feasible(model)
        println("Optimization converged to the optimal solution.")
        optimalcount += 1;
        optimalvals = push!(optimalvals,num)
    else
        println("Optimization stopped due to local infeasibility.")
        localinfeasvals = push!(localinfeasvals,num)
    end
    # Write the trajectory to a file
    open("MonteCarlo/trial5/test" * string(num) * ".txt", "a") do io
        for i in 1:2
            for j in 1:T+1
                println(io, states[i, j, 1], ",", states[i, j, 2], ",", states[i, j, 3], ",", states[i, j, 4])
            end
        end
    end
end

function rand_float(a, b)
    return a + (b - a) * rand()
end

Random.seed!(200)
optimalvalues = Set()
localinfeasvalues = Set()
optimalcount = 0;
for i in 1:100
    a = -5
    b = 5
    x1 = rand_float(a,b)
    x2 = rand_float(a,b)
    x3 = rand_float(a,b)
    x4 = rand_float(a,b)
    y1 = rand_float(a,b)
    y2 = rand_float(a,b)
    y3 = rand_float(a,b)
    y4 = rand_float(a,b)
    initial_pose1 = [x1,y1,1,0]
    initial_pose2 = [x2,y2,1,0]
    final_pose1 = [x3,y3,1,0]
    final_pose2 = [x4,y4,1,0]
    distance1 = (x2-x1)^2+(y2-y1)^2
    distance2 = (x4-x3)^2+(y4-y3)^2
    distance = min(distance1,distance2)
    obstacle = sqrt(distance)*rand()
    if (i == 1)
        open("MonteCarlo/trial5/radius.txt", "w") do io
            println(io,obstacle)
        end
    else
        open("MonteCarlo/trial5/radius.txt", "a") do io
            println(io,obstacle)
        end
    end
    optimizerObstacle(i, initial_pose1,initial_pose2,final_pose1,final_pose2, obstacle,optimalcount,optimalvalues,localinfeasvalues)
end
print(optimalcount/1000.0)
print(optimalvalues)
print(localinfeasvalues)