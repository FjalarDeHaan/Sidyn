using StaticArrays
using Random

struct Autophore{N, Float64} <: StaticVector{N, Float64}
    values::NTuple{N, Float64}
end

@inline Base.getindex(ψ::Autophore, i::Int64) = ( @boundscheck checkbounds(ψ,i)
                                                ; ψ.values[Base.to_index(i)] )

function (::Type{Autophore{N}})(d::NTuple{N,Float64}) where {N, Float64}
    return Autophore{N,Float64}(d)
end

StaticArrays.similar_type( ::Type{A}
                         , ::Type{Float64}
                         , ::S ) where {A<:Autophore, S<:(Size)} = Autophore

function Autophore(values::AbstractVector; normalised=true)
    if normalised
        return Autophore{length(values)}(values/norm(values))
    else
        return Autophore{length(values)}(values)
    end
end

function Autophore(rng::AbstractRNG, n::Int)
    return Autophore(rand(rng, n) .- .5)
end

function Autophore(n::Int)
    return Autophore(rand(n) .- .5)
end

function Base.show(io::IO, ψ::Autophore)
    if norm(ψ) > 1
    println( io
           , length(ψ.values)
           , "-dimensional identity vector "
           , "(norm ≈ ", (@sprintf "%0.2f" norm(ψ)), "):" )
        for i in 1:length(ψ.values)
            if ψ.values[i] < 0
                println(io, " ", ψ.values[i])
            else
                println(io, "  ", ψ.values[i])
            end
        end
    else
    println( io
           , length(ψ.values)
           , "-dimensional identity vector (normalised):" )
        for i in 1:length(ψ.values)
            s = string(ψ.values[i])
            if ψ.values[i] == 1.0
                println(io, "  1")
            elseif ψ.values[i] == -1.0
                println(io, " -1")
            elseif ψ.values[i] == 0.0
                println(io, "  0")
            else
                if s[1] == string(-1)[1]
                    println(io, " - ", s[3:end])
                else
                    println(io, "   ", s[2:end])
                end
            end
        end
    end
end

Base.display(ψ::Autophore) = show(ψ)

function normalise(ψ::Autophore)
    return ψ / norm(ψ)
end
