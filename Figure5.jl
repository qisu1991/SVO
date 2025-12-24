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

# get the trajectory by Monte Carlo Simulations
# output the cooperation rate throughtout the direct evolutionary process
function xC_direct_trajectory(edge_seq::Array{Int64,2}, b::Float64, delta::Float64, generation::Int64, sample_interval::Int64, mutation::Float64)
    N = maximum(edge_seq[:,1:2])
    G = sparse(edge_seq[:,1], edge_seq[:,2], edge_seq[:,3], N, N)
    deg = sum(G, dims = 2)
    P = spzeros(Float64, N, N)
    P[:,:] = G./deg

    xC_overall = 0.0 
    str_seq = ones(Float64,N)*rand([0 1])

    for g = 1:generation
        if g%sample_interval == 0
            xC = sum(str_seq)/N
            xC_overall = xC_overall + xC
        end

        node_updated = rand(1:N)
        node_competitor = findall(!iszero, P[node_updated,:])
        ismutation = 0 
        if rand() < mutation
            ismutation = 1 
        end

        if length(findall(isequal(str_seq[node_updated]),str_seq[node_competitor])) == length(node_competitor) && ismutation == 0
            continue
        end 
        
        payoff_list = -str_seq + b*P*str_seq
        payoff_list = exp.(delta*payoff_list)
        payoff_competitor = P[node_updated,:].*payoff_list
        threshold = rand(1)*sum(payoff_competitor)
        sum1 = 0
        for j = 1:length(node_competitor)
            sum1 = sum1 + payoff_competitor[node_competitor[j]]
            if sum1 >= threshold[1]
                if ismutation == 0
                    str_seq[node_updated] = str_seq[node_competitor[j]]
                else
                    str_seq[node_updated] = rand([0 1])
                end
                break
            end
        end
    end

    return (xC_overall*sample_interval)/generation
end

# get the trajectory by Monte Carlo Simulations
# output one's psi throughtout the indirect evolutionary process
function xC_indirect_trajectory(edge_seq::Array{Int64,2}, b::Float64, delta::Float64, delta_rep::Float64, generation::Int64, sample_interval::Int64, mutation::Float64)
    N = maximum(edge_seq[:,1:2])
    G = sparse(edge_seq[:,1], edge_seq[:,2], edge_seq[:,3], N, N)
    deg = sum(G, dims = 2)
    P = spzeros(Float64, N, N)
    P[:,:] = G./deg
    P2 = P^2

    xC_overall = 0.0 
    psi_overall = 0.0
    psi_all = -pi*9/10:pi/10:pi 
    psi_seq = ones(Float64,N)*psi_all[rand(1:20)]

    for g = 1:generation
        if g%sample_interval == 0
            xC = sum(-cos.(psi_seq) + b*sin.(psi_seq).*diag(P2))/N
            xC = xC*delta/4 + 1/2
            xC_overall = xC_overall + xC
        
            psi_mean = sum(psi_seq)/N
            psi_overall = psi_overall + psi_mean
        end

        node_updated = rand(1:N)
        node_competitor = findall(!iszero, P[node_updated,:])
        ismutation = 0 
        if rand() < mutation
            ismutation = 1 
        end

        if length(findall(isequal(psi_seq[node_updated]),psi_seq[node_competitor])) == length(node_competitor) && ismutation == 0
            continue
        end 
        x = -cos.(psi_seq) + b*sin.(psi_seq).*diag(P2)
        payoff_list = (-x + P*x*b)*delta/4 .+ (b-1)/2 
        payoff_list = exp.(delta_rep*payoff_list)
        payoff_competitor = P[node_updated,:].*payoff_list
        threshold = rand(1)*sum(payoff_competitor)
        sum1 = 0
        for j = 1:length(node_competitor)
            sum1 = sum1 + payoff_competitor[node_competitor[j]]
            if sum1 >= threshold[1]
                if ismutation == 0
                    psi_seq[node_updated] = psi_seq[node_competitor[j]]
                else
                    psi_seq[node_updated] = psi_all[rand(1:20)]
                end
                break
            end
        end
    end

    return (xC_overall*sample_interval)/generation, (psi_overall*sample_interval)/generation
end

tick()
N = 12
delta = 0.1
beta = 0.1
mutation = 0.001
generation = 100000000
sample_interval = 1000

# generate a two-community network as shown in Figure 4a
edge_seq = two_community(N)

# obtain the cooperation frequency in the direct approach
b = 10.0
xC_direct_average = xC_direct_trajectory(edge_seq, b, delta, generation, sample_interval, mutation)

xC_indirect_average, psi_average = xC_indirect_trajectory(edge_seq, b, beta, delta, generation, sample_interval, mutation)

println("xC_direct_average = ", xC_direct_average)
println("xC_indirect_average = ", xC_indirect_average)
tock()