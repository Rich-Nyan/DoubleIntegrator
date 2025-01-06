using Ipopt
using Optim
using LinearAlgebra
using ForwardDiff
using Plots
using JuMP

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

initial_pose1 = [0,0,1,0]
final_pose1 = [4,4,1,0]
initial_pose2 = [4,0,1,0]
final_pose2 = [0,4,1,0]

obstacle = 2
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
function dynamic_feasibility(z, model, T)

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

function gradient_f(z,model, α)
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

function file_to_array(str::AbstractString)
    # Initialize an empty array
    arr = []
    inde = 1
    
    open(str, "r") do file
        reader = readlines(str)
        lines = length(reader)
        for i in 1:lines
            input = split(reader[i])
            for j in 1:length(input)
                z_guess[inde] = parse(Float64,input[j])
                inde += 1
            end
        end
    end

    return arr
end
#Optimize Function

function optimizer(file, α)
    model = Model(optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 5))

    file_to_array("DI_z-guess.txt")
    
    @variable(model, z[i = 1:12 * T])

    for i in 1:length(z)
        if i <= 6 * T
            set_start_value(z[i], z_guess[i])
        else
            set_start_value(z[i], z_guess[10 * T + 4 + i - 6*T])
        end
    end

    
    h_z = dynamic_feasibility(z, model, T)
    # Objective
    # @objective(model, Min, objective_function(z))

    # Lagrange gradient
    grad_f = gradient_f(z, model, α)
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
    
    # Write the trajectory to a file
    open(file, "a") do io
        for i in 1:2
            for j in 1:T+1
                println(io, states[i, j, 1], ",", states[i, j, 2], ",", states[i, j, 3], ",", states[i, j, 4])
            end
        end
    end
end

α = [0.01]
for i in α
    optimizer("2playertrajectory_1.txt",i)
end