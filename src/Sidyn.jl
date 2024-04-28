module Sidyn

using Graphs, SimpleWeightedGraphs
using LinearAlgebra

function hierarchy(n=5, teamsize=5)
    # First make a random tree on n vertices.
    g = uniform_tree(n) |> DiGraph
    # First highest degree vertex becomes the root.
    root = g |> degree |> argmax
    # Identify the leaves of the tree and their number.
    leaves = findall(==(2), degree(g))
    nleaves = leaves |> length
    # Add cliques to each leaf.
    cliquev = n + 1 # Clique-vertex to connect to.
    for leaf in leaves
        g = blockdiag(g, complete_digraph(teamsize)) # Add clique as component.
        add_edge!(g, leaf, cliquev) # Connect the two components one way.
        add_edge!(g, cliquev, leaf) # Connect the two components the other way.
        cliquev += teamsize # The next clique-vertex to connect to.
    end
    # Get the adjacency matrix.
    W = g |> adjacency_matrix |> collect |> float
    # Distances from root.
    d = gdistances(g, root)
    # Weight the edges.
    W[root, root] = 1.0 # Root listens to root.
    for edge in edges(g)
        # Pointing _down_ the hierarchy.
        if d[edge.dst] - d[edge.src] > 0
            W[edge.src, edge.dst] = 0.5
        end
        # Pointing _up_ the hierarchy.
        if d[edge.dst] - d[edge.src] < 0
            W[edge.src, edge.dst] = 5.0
        end
        # Pointing _level_ in the hierarchy.
        if d[edge.dst] - d[edge.src] == 0
            W[edge.src, edge.dst] = 1.0
        end
    end
    # Normalise the rows to make the matrix stochastic.
    for i ∈ 1:size(W)[1]
        W[i, :] /= sum(W[i, :])
    end
    # Insert the external identities.
    g = WDiGraph(W)
    for team in 1:nleaves
        v_team = n + team*teamsize # The last vertex in the team.
        add_vertex!(g) # Add the external identity vertex.
        v_extid = size(g)[1] # Get the external identity vertex.
        add_edge!(g, v_extid, v_extid, 1.0) # Listens only to itself.
        add_edge!(g, v_team, v_extid, 1.0) # Is listened to.
        # Re-normalise --- v_extid gets listened to by .5.
        for v in outneighbors(g, v_team)
            w_new = get_weight(g, v_team, v) / 2.0 # Two, as weight of 1 added.
            add_edge!(g, v_team, v, w_new)
        end
    end
    # Deliver.
    return g
end


function strongly_ER(n=5)
    g = DiGraph()
    while !is_strongly_connected(g)
        # Find a random number of edges between n and n(n-1).
        m = rand(n:(n/2)*(n-1))
        # Generate an Erdos-Renyi graph accordingly and make it a digraph.
        g = erdos_renyi(n, m) |> DiGraph
    end
    # Deliver.
    return g
end

angle(a, b) = acos(a ⋅ b / (norm(a)*norm(b)))
ϕ = angle

polar(r, phi) = r * [ cos(phi), sin(phi) ]

function ru(ndims, i=1, bias=nothing) # .3)
    v = randn(ndims)/10
    if !isnothing(bias)
        v[i] = (randn(1)[1])/5 +  bias
    end
    v /= norm(v)
    return v
end

function collapse(ψ::Vector, dims)
    # Prepare the new vector with zeroes everywhere.
    z = zeros(length(ψ))
    # Take over the entries of the selected dimensions.
    for i ∈ dims
        z[i] = ψ[i]
    end
    # Re-normalise the result and deliver.
    return z / norm(z)
end

"""
    project(ψ::Vector, A::Matrix; yield=:coordinates)

Project the identity vector `ψ` onto the column space of `A`. By default the
coordinates against the columns of `A` are returned. To instead have the actual
projection vector returned, pass `yield=:vector.
"""
function project(ψ::Vector, A::Matrix; yield=:coordinates)
    # Determine whether coordinates against A or the projection of ψ are needed.
    check = 0
    yield == :coordinates ?  P = ((A'*A)^-1) * A'     : check += 1
    yield == :vector      ?  P = A * ((A'*A)^-1) * A' : check += 1
    if check > 1
        error("Wrong value for argument `yield`: $yield")
    end
    # Compute coordinates of v against A's columns or the vector they represent.
    z = P*ψ
    # Deliver.
    return z
end

"""
    collapse(ψ::Vector, A::Matrix; yield=:coordinates)

Collapse the identity vector `ψ` onto the column space of `A`. Collapsing is
orthogonal projection followed by rescaling to unit euclidean length. By
default the coordinates against the columns of `A` are returned. To instead
have the actual projection vector returned, pass `yield=:vector.
"""
function collapse(ψ::Vector, A::Matrix; yield=:coordinates)
    # Project, re-normalise, deliver.
    return project(ψ, A, yield=yield) / norm(project(ψ, A, yield=:vector))
end


# Export everything this module has to offer --- temporary convenience.
for n in names(@__MODULE__; all=true)
    if Base.isidentifier(n) && n ∉ (Symbol(@__MODULE__), :eval, :include)
        @eval export $n
    end
end


end # module Sidyn
