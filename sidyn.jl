using Sidyn
using Graphs, SimpleWeightedGraphs, LinearAlgebra

W = (complete_digraph(10) |> adjacency_matrix) / 9
W∞ = W ^ 10000
Ψ = [ ru(8, 2)'; ru(8, 2)'; ru(8, 2)'; ru(8, 2)'; ru(8, 2)'
    ; ru(8, 2)'; ru(8, 2)'; ru(8, 2)'; ru(8, 2)'; ru(8, 2)' ]
ψˣ = (W∞ * Ψ)[1, :]
