# simulate the stationary distribution 
# using Monte Carlo simulations
################################################################################

using Printf
using Statistics
using LinearAlgebra
using DelimitedFiles
using TickTock
using IterativeSolvers
using SparseArrays

rng = MersenneTwister(1234)

# generate a two-community network
function two_community(N::Int64)
    N1 = Int64(round(0.5*N))
    edge_seq = zeros(Int64,N1*(N1+1),3)
    id = 1
    for i = 1:N1-1
        for j = i+1:N1
            edge_seq[id,1] = i
            edge_seq[id,2] = j
            edge_seq[id,3] = 1
            id = id + 1
            edge_seq[id,1] = j
            edge_seq[id,2] = i
            edge_seq[id,3] = 1
            id = id + 1
        end
    end

    for i = 2:N1
        edge_seq[id, :] = [N1+1 N1+i 1]
        edge_seq[id+1, :] = [N1+i N1+1 1]
        id = id+2
    end
    edge_seq[id, :] = [N1+1 N1 1]
    edge_seq[id+1, :] = [N1 N1+1 1]
    return edge_seq
end

# function to get reproductive value of each individual in each layer
# input: edge_seq is the replacement graph
# any directed structure
function pi_directed_DB(edge_seq::Array{Int64,2})
    N = maximum(edge_seq[:,1:2])
    M = sparse(edge_seq[:,1],edge_seq[:,2],edge_seq[:,3],N,N)
    M = M./sum(M,dims=1)

    MatA = spzeros(Float64, N, N)
    MatA = M - sparse(I,N,N)
    MatB_reduced = -MatA[1:N-1,N]
    MatA = MatA - MatA[:,N]*ones(Float64,1,N)
    MatA_reduced = MatA[1:N-1,1:N-1]

    pi_solution_reduced = idrs(MatA_reduced,MatB_reduced)

    pi_solution = zeros(Float64,N)
    pi_solution[1:N-1] = pi_solution_reduced
    pi_solution[N] = 1.0-sum(pi_solution)

    return pi_solution
end

# solve a set of equations to obtain beta_i1i2_j1j2
# for rhoA > rhoB under a uniform distribution of initial mutant
# input: edge_seq is based on replacement graph
function beta_directed_DB_uniform(edge_seq::Array{Int64,2})
    N = maximum(edge_seq[:,1:2])
    M = sparse(edge_seq[:,1],edge_seq[:,2],edge_seq[:,3],N,N)
    M = M./sum(M,dims=1)

    MT = transpose(M)
    inds = findall(!iszero, MT)
    a = getindex.(inds, 1)
    b = getindex.(inds, 2)
    vec = [i for i = 1:N]
    vec = reshape(vec,1,N)
    X1 = N*(a.-1)*ones(Int64,1,N)+ones(Int64,length(a))*vec
    Y1 = N*(b.-1)*ones(Int64,1,N)+ones(Int64,length(b))*vec
    W1 = MT[inds]*ones(Int64,1,N)

    vec = [(i-1)*N for i = 1:N]
    vec = reshape(vec,1,N)
    X2 = a*ones(Int64,1,N)+ones(Int64,size(a,1))*vec
    Y2 = b*ones(Int64,1,N)+ones(Int64,size(b,1))*vec
    W2 = MT[inds]*ones(Int64,1,N)

    X11 = reshape(X1,:,1)
    Y11 = reshape(Y1,:,1)
    W11 = reshape(W1,:,1)
    X22 = reshape(X2,:,1)
    Y22 = reshape(Y2,:,1)
    W22 = reshape(W2,:,1)
    matA = sparse(X11[:,1],Y11[:,1],W11[:,1],N^2,N^2)+sparse(X22[:,1],Y22[:,1],W22[:,1],N^2,N^2)

    matA = matA/2 - sparse(I, N^2, N^2)
    vec = setdiff(1:N^2,[(i-1)*N+i for i = 1:N])
    matA_reduced = matA[vec, vec]
    matB_reduced = -ones(Float64,N^2-N)

    beta_solution_reduced = idrs(matA_reduced,matB_reduced)
    beta_solution = spzeros(Float64, N^2)
    beta_solution[vec] = beta_solution_reduced
    beta_solution = reshape(beta_solution, N, N)

    return transpose(beta_solution)
end

# get (b/c)^* under death-birth updating
# donation and dispersal the same direction
# input: edge_seq is based on interaction graph
function bc_directed_DB_uniform_dirsame(edge_seq::Array{Int64,2})
    pi = pi_directed_DB(edge_seq)
    beta = beta_directed_DB_uniform(edge_seq)
    N = maximum(edge_seq[:,1:2])
    M = sparse(edge_seq[:,1],edge_seq[:,2],edge_seq[:,3],N,N)
    Mp = M./sum(M,dims=2)

    n1 = sum(pi.*beta.*Mp)
    n3 = sum(pi.*beta.*Mp^3)
    n2 = sum(pi.*beta.*Mp^2)

    return n2/(n3-n1), n1, n2, n3
end

# get (b/c)^* under death-birth updating
# donation and dispersal the same direction
# input: edge_seq is based on interaction graph
function eta_structural_coefficient(edge_seq::Array{Int64,2})
    pi = pi_directed_DB(edge_seq)
    beta = beta_directed_DB_uniform(edge_seq)
    N = maximum(edge_seq[:,1:2])
    M = sparse(edge_seq[:,1],edge_seq[:,2],edge_seq[:,3],N,N)
    Mp = M./sum(M,dims=2)
    pjj2 = Mp^2
    MM = diag(pjj2)

    n1 = sum(pi.*(beta.*transpose(MM)).*Mp)
    n3 = sum(pi.*(beta.*transpose(MM)).*Mp^3)
    n2 = sum(pi.*(beta.*transpose(MM)).*Mp^2)

    return n1, n2, n3
end

# input: a network
# output: identify the category of the network input (see Extended Data Figure 3)
function network_category(edge_seq::Array{Int64,2})
    N = maximum(edge_seq[:,1])
    G = sparse(edge_seq[:,1], edge_seq[:,2], edge_seq[:,3], N, N)
    deg = sum(G, dims = 2)
    P = spzeros(Float64, N, N)
    P[:,:] = G./deg
    P2 = P^2
    p2 = sum(diag(P2))/N

    bc, n01, n02, n03 = bc_directed_DB_uniform_dirsame(edge_seq)
    n11, n12, n13 = eta_structural_coefficient(edge_seq)

    error = 10^(-10)
    n01 = abs(n01)>error ? n01 : 0
    n02 = abs(n02)>error ? n02 : 0
    n03 = abs(n03)>error ? n03 : 0
    n11 = abs(n11)>error ? n11 : 0
    n12 = abs(n12)>error ? n12 : 0
    n13 = abs(n13)>error ? n13 : 0
    b0 = 0 
    b1 = 0
    if abs(n03-n01) > error
        b0 = n02/(n03-n01)
    end
    if abs(n13-n11) > error
        b1 = n12/(n13-n11)
    end
    
    # case 1
    sign = 0 
    b_critical = 0
    if n03-n01 > error && n13-n11 > error
        if b0-b1 > error
            sign = 11
            p = -b1^2/3 + (n02*b1)/(n12*b0*p2)
            q = -2*b1^3/27 + (n02*b1^2)/(3*n12*b0*p2) - (n02*b1)/(n12*p2)
            D = p^3/27 + q^2/4
            b_critical = b1/3 + cbrt(-q/2+D^(1/2)) + cbrt(-q/2-D^(1/2))
        elseif b0-b1 < -error 
            sign = 12
        else 
            sign = 13 
        end
    end
    # case 2
    if n03-n01 > error && n13-n11 < -error
        sign = 2
    end
    # case 3 
    if n03-n01 > error && abs(n13-n11) < error
        sign = 3
    end
    # case 4
    if n03-n01 < -error && n13-n11 > error
        p = -b1^2/3 + (n02*b1)/(n12*b0*p2)
        q = -2*b1^3/27 + (n02*b1^2)/(3*n12*b0*p2) - (n02*b1)/(n12*p2)
        D = p^3/27 + q^2/4
        if D > 0
            sign = 41
            b_critical = b1/3 + cbrt(-q/2+D^(1/2)) + cbrt(-q/2-D^(1/2))
        else 
            sign = 42
            b_critical = b1/3 + 2*(-p/3)^(1/2)*cos(acos((3*q)/(2*p)*(-3/p)^(1/2))/3)
        end
    end
    # case 5
    if n03-n01 < -error && n13-n11 < -error
        if b0-b1 > error
            sign = 51
        elseif b0-b1 < -error 
            sign = 52
        else 
            sign = 53 
        end 
    end
    # case 6
    if n03-n01 < -error && abs(n13-n11) < error
        sign = 6
    end
    # case 7
    if abs(n03-n01) < error && n13-n11 > error
        sign = 7
        p = -b1^2/3
        q = -2*b1^3/27 - (n02*b1)/(n12*p2)
        D = p^3/27 + q^2/4
        b_critical = b1/3 + cbrt(-q/2+D^(1/2)) + cbrt(-q/2-D^(1/2))
    end
    # case 8
    if abs(n03-n01) < error && n13-n11 < -error
        sign = 8
    end
    # case 9
    if abs(n03-n01) < error && abs(n13-n11) < error
        sign = 9
    end

    return sign, b_critical, b0, b1, n01, n02, n03, n11, n12, n13
end

# get fixation probability by Monte Carlo Simulations
# get fixation probability of A-individuals by Monte Carlo Simulations
function rhoA_simulation_DB(edge_seq::Array{Int64,2}, psi1::Float64, psi2::Float64, b::Float64, delta::Float64, delta_rep::Float64, generation::Int64, sample_times::Int64)
    N = maximum(edge_seq[:,1:2])
    G = sparse(edge_seq[:,1], edge_seq[:,2], edge_seq[:,3], N, N)
    deg = sum(G, dims = 2)
    P = spzeros(Float64, N, N)
    P[:,:] = G./deg
    P2 = P^2

    fixation_times = 0
    # simulation
    for t = 1:sample_times 
        # initial SVO designation
        SVO_list = zeros(Int64,N)
        SVO_list[rand(1:N)] = 1

        for g = 1:generation
            # selecting the node to be updated
            node_updated = rand(1:N)
            node_competitor = findall(!iszero, P[node_updated,:])

            # check if the population enters into an aborbing state
            if g % 1000 == 0
                if sum(SVO_list) == N
                    fixation_times = fixation_times + 1
                    # gg = open("fixation_time","a")
                    # writedlm(gg,[1 g])
                    # close(gg)
                    break
                end
                if sum(SVO_list) == 0
                    # gg = open("fixation_time","a")
                    # writedlm(gg,[0 g])
                    # close(gg)
                    break
                end
            end
            if sum(SVO_list[node_competitor]) == length(node_competitor)*SVO_list[node_updated]
                continue
            end

            # calculate the payoff
            psi_seq = SVO_list*psi1 + (-SVO_list .+ 1)*psi2
            x = -cos.(psi_seq) + b*sin.(psi_seq).*diag(P2)
            payoff_list = (-x + P*x*b)*delta/4 .+ (b-1)/2 
            payoff_list = exp.(delta_rep*payoff_list)

            payoff_competitor = P[node_updated,:].*payoff_list
            threshold = rand(1)*sum(payoff_competitor)
            sum1 = 0
            for j = 1:length(node_competitor)
                sum1 = sum1 + payoff_competitor[node_competitor[j]]
                if sum1 >= threshold[1]
                    SVO_list[node_updated] = SVO_list[node_competitor[j]]
                    break
                end
            end
        end
    end

    return fixation_times/sample_times
end

# get fixation probability by Monte Carlo Simulations
# get fixation probability of B-individuals by Monte Carlo Simulations
function rhoB_simulation_DB(edge_seq::Array{Int64,2}, psi1::Float64, psi2::Float64, b::Float64, delta::Float64, delta_rep::Float64, generation::Int64, sample_times::Int64)
    N = maximum(edge_seq[:,1:2])
    G = sparse(edge_seq[:,1], edge_seq[:,2], edge_seq[:,3], N, N)
    deg = sum(G, dims = 2)
    P = spzeros(Float64, N, N)
    P[:,:] = G./deg
    P2 = P^2

    fixation_times = 0
    # simulation
    for t = 1:sample_times 
        # initial SVO designation
        SVO_list = ones(Int64,N)
        SVO_list[rand(1:N)] = 0

        for g = 1:generation
            # selecting the node to be updated
            node_updated = rand(1:N)
            node_competitor = findall(!iszero, P[node_updated,:])

            # check if the population enters into an aborbing state
            if g % 1000 == 0
                if sum(SVO_list) == N
                    break
                end
                if sum(SVO_list) == 0
                    fixation_times = fixation_times + 1
                    break
                end
            end
            if sum(SVO_list[node_competitor]) == length(node_competitor)*SVO_list[node_updated]
                continue
            end

            # calculate the payoff
            psi_seq = SVO_list*psi1 + (-SVO_list .+ 1)*psi2
            x = -cos.(psi_seq) + b*sin.(psi_seq).*diag(P2)
            payoff_list = (-x + P*x*b)*delta/4 .+ (b-1)/2 
            payoff_list = exp.(delta_rep*payoff_list)

            payoff_competitor = P[node_updated,:].*payoff_list
            threshold = rand(1)*sum(payoff_competitor)
            sum1 = 0
            for j = 1:length(node_competitor)
                sum1 = sum1 + payoff_competitor[node_competitor[j]]
                if sum1 >= threshold[1]
                    SVO_list[node_updated] = SVO_list[node_competitor[j]]
                    break
                end
            end
        end
    end

    return fixation_times/sample_times
end


tick()
N = 12
beta = 0.01
delta = 0.1
generation = 1000000
sample_times = 100000000

# generate a two-community network as shown in Figure 4a
edge_seq = two_community(N)
sign1, b_critical, b0, b1, n01, n02, n03, n11, n12, n13 = network_category(edge_seq)

b_example = 10.0
# obtain the fittest SVO
tan_phi = -b_example*(-n12+b_example*(n13-n11))/(-n02+b_example*(n03-n01));
psi_critical = atan(tan_phi)

# run the simulation to get the fixation probability of psi_instance over psi_critical and psi_critical over psi_instance
psi_instance = -pi/2
rho_critical = rhoA_simulation_DB(edge_seq, psi_critical, psi_instance, b_example, beta, delta, generation, sample_times)
rho_instance = rhoB_simulation_DB(edge_seq, psi_critical, psi_instance, b_example, beta, delta, generation, sample_times)
tock()

println("rho_critical = ", rho_critical)
println("rho_instance = ", rho_instance)