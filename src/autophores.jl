struct Autophore <: AbstractVector{Float64}
    v::Vector{Float64} # Components of the autophore.
    λ::Float64 # Susceptibility, i.e. 1 - stubbornness.
    Autophore(v, λ) = new(v/norm(v), λ)
end

Autophore(v::Vector) = Autophore(v, 1.0)

Base.size(ψ::Autophore) = (length(ψ.v),)
Base.getindex(ψ::Autophore, i::Int) = ψ.v[i]

function Base.show(io::IO, ψ::Autophore)
    println( io
           , length(ψ.v), "-dimensional identity vector "
           , "(λ = ", ψ.λ, "):" )
    for i in 1:length(ψ.v)
        if ψ.v[i] == 1.0
            println(io, " 1") #, @sprintf "%0.6f" ψ.v[i])
        elseif ψ.v[i] == 0.0
            println(io, " 0")
        else
            println(io, "  ", (@sprintf "%0.6f" ψ.v[i])[2:end])
        end
    end
end

Base.display(ψ::Autophore) = show(ψ)
