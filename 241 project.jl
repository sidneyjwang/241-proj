using LinearAlgebra

# the probability that the user actually clicks on a link from the current page
# instead of picking one at random
const DAMPING = 0.99
const EPSILON = 10^-5

# function to compare equality of two floats
function FloatEqual(x1::ComplexF64, x2::ComplexF64)
    return abs(x1 - x2) < EPSILON
end

# creates a map of each site to its index in the input
function MkIndexMap(sites::Vector{String}) 
    result = Dict{String, Int}()
    for i in 1:size(sites, 1)
        result[sites[i]] = i
    end
    return result
end

# this function takes in a list of websites, 
# and a list containing a list of outgoing links for each website (adjacency list),
# returning an unweighted Adjacency Matrix representation of the graph they form
# where 1 = directed edge from col index -> row index and 0 = no edge
# links is assumed to contain [] for a specific vertex's adjacency list if it has
# no out-neighbors
function MkGraph(sites::Vector{String}, links::Vector{Vector{Any}})
    # initialize an empty matrix of zeros
    result = zeros(size(sites, 1), size(sites, 1))
    indexMap = MkIndexMap(sites)
    for i in 1:size(links, 1)
        for j in 1:length(links[i])
            result[indexMap[links[i][j]], i] = 1
        end
    end
    return result
end    

# this function takes in an adjacency matrix and returns the weighted transition
# matrix; normalizes all col sums to be 1
# takes into account random jumps made by the user from page A to an unlinked page
# invariant: every row in the returned matrix has sum 1
function MkTransition(graph::Matrix{Float64})
    len = size(graph, 1)
    # keep track of how many vertices are not immediate neighbors of v
    zeroes = zeros(len)
    # total probability of leaving a page is 1, so sum of each row should be 1
    for i in 1:len # column
        curr = sum(graph[:, i])
        for j in 1:len # row
            if graph[j, i] != 0 # prevent zero divisions if a vertex is isolated
                # normalize all values and account for randomness damping factor
                graph[j, i] /= curr
                graph[j, i] *= DAMPING
            else
                zeroes[i] += 1
            end
        end
    end
    # now add in the possibility to jump to a random page, for anything that's not linked
    # any "ghost edge" should have a weight of (1-damping) / (# of ghost edges out of v)
    # unless it's a lone vertex-- then visit any other page with equal probability
    for i in 1:len # column
        for j in 1:len # row
            if graph[j, i] == 0 && zeroes[i] != len
                graph[j, i] = (1-DAMPING) / zeroes[i]
            elseif graph[j, i] == 0
                # disconnected vertex, jump to a random page
                graph[j, i] = 1/len
            end
        end
    end
    return graph
end

function FindXIdx(list::Vector{ComplexF64}, x::ComplexF64)
    for i in 1:length(list)
        if FloatEqual(x, list[i])
            return i
        end
    end
end

# PageRank takes in a series of indexed websites, as well as a vector containing
# vectors of pages that each indexed site links to
function PageRank(sites::Vector{String}, links::Vector{Vector{Any}})
    graph = MkGraph(sites, links)
    transition = MkTransition(graph)
    val, vec = eigen(transition)
    # identify the eigenvector associated w/ eigenvalue 1-- this is the solution
    # process the data into a well formed ranking of the pages
    idx = FindXIdx(convert(Vector{ComplexF64}, val), convert(ComplexF64, 1.0))
    soln = convert(Vector{Float64}, vec[:, idx])
    # soln contains the pageranks for each page
    # now, convert them to a ranked order
    enumerate = [(sites[i], soln[i]) for i in 1:length(sites)]
    final = sort(enumerate, by = x -> x[2], rev=true)
    return [final[i][1] for i in 1:length(final)]
end













# decide whether to give more impact to pages that are hubs vs authorities
const HUB_WEIGHT = 0.2
const AUTHORITY_WEIGHT = 1

# builds a hash table of vertex -> all its in-neighbors
function MkInGraph(sites::Vector{String}, links::Vector{Vector{Any}})
    res = Dict{String, Set{String}}()
    for i in 1:length(sites)
        curr = links[i]
        for j in 1:length(curr)
            if !haskey(res, curr[j])
                res[curr[j]] = Set()
            end
            push!(res[curr[j]], sites[i])
        end
    end
    return res
end

# builds a hash table of vertex -> all its out-neighbors
function MkOutGraph(sites::Vector{String}, links::Vector{Vector{Any}})
    res = Dict{String, Set{String}}()
    for i in 1:length(sites)
        res[sites[i]] = Set(links[i])
    end
    return res
end

# sum of hub scores of x's in-neighbors
function GetAuthority(x::String, sites::Vector{String}, inGraph::Dict{String, Set{String}}, hub::Vector{Float64})
    total = 0
    for i in 1:length(sites)
        if sites[i] in inGraph[x]
            total += hub[i]
        end
    end
    return total
end

# sum of authority scores of out-neighbors
function GetHub(x::String, sites::Vector{String}, outGraph::Dict{String, Set{String}}, authority::Vector{Float64})
    total = 0
    for i in 1:length(sites)
        if sites[i] in outGraph[x]
            total += authority[i]
        end
    end
    return total
end

function HITS(sites::Vector{String}, links::Vector{Vector{Any}})
    iNbors = MkInGraph(sites, links)
    oNbors = MkOutGraph(sites, links)
    # initialize all hub and authority values to be 1
    authority = normalize([1.0 for i in 1:length(sites)])
    hub = normalize([1.0 for i in 1:length(sites)])
    while true
        newAuthority = normalize([GetAuthority(sites[i], sites, iNbors, hub) for i in 1:length(authority)])
        newHub = normalize([GetHub(sites[i], sites, oNbors, newAuthority) for i in 1:length(hub)])
        if abs(sum(newAuthority) - sum(authority)) + abs(sum(newHub) - sum(hub)) < EPSILON
            authority = newAuthority
            hub = newHub
            break
        end 
        authority = newAuthority
        hub = newHub
    end
    # account for a vertex's "importantness" by taking the sum of its two components
    println(authority)
    println(hub)
    final = [(AUTHORITY_WEIGHT * authority[i] + HUB_WEIGHT * hub[i], sites[i]) for i in 1:length(sites)]
    sorted = sort(final, by = x -> x[1], rev=true)
    return [sorted[i][2] for i in 1:length(sites)]
end

sites = ["A", "B", "C", "D", "E", "F"]
links = [[], ["C"], ["B"], ["A", "B"], ["D", "B", "F"], ["E", "B"]]

# println(MkInGraph(sites, links))
# println(MkOutGraph(sites, links))

println(HITS(sites, links))
println(PageRank(sites, links))