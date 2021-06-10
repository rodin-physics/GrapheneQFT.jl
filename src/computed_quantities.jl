include("defects.jl")
"""
    δG_R(z::ComplexF64,
     pairs::Vector{Tuple{GrapheneCoord,GrapheneCoord}},
     s::GrapheneSystem)

The correction to the real-space graphene Green's function in the presence of
defects as a function of complex energy `z`.

The function returns a vector of `ComplexF64` for each [`GrapheneCoord`](@ref)
    in `pairs`.

# Arguments
* `z`: complex energy
* `pairs`: pairs of [`GrapheneCoord`](@ref)'s for which `δG_R` is calculated
* `s`: [`GrapheneSystem`](@ref) for which `δG_R` is calculated
"""
function δG_R(
    z::ComplexF64,
    pairs::Vector{Tuple{GrapheneCoord,GrapheneCoord}},
    s::GrapheneSystem,
)
    # Calculate the scattering matrix D which is the same
    # for every pair of coordinates
    prop_mat = propagator_matrix(z, s.scattering_atoms)
    if length(s.imps) == 0
        D =
            s.Δ *
            inv(Diagonal(ones(length(s.scattering_atoms))) .- prop_mat * s.Δ)
    else
        Γ0 = 1 ./ (z .- s.imps) |> Diagonal
        D =
            (s.Δ .+ s.V * Γ0 * adjoint(s.V)) * inv(
                Diagonal(ones(length(s.scattering_atoms))) .-
                prop_mat * (s.Δ .+ s.V * Γ0 * adjoint(s.V)),
            )
    end

    PropVectorRs = [
        [graphene_propagator(x, p[2], z) for x in s.scattering_atoms] for
        p in pairs
    ]

    PropVectorLs = [
        pairs[idx][1] == pairs[idx][2] ? permutedims(PropVectorRs[idx]) :
        [
            graphene_propagator(pairs[idx][1], x, z) for
            x in s.scattering_atoms
        ] |> permutedims for idx = 1:length(pairs)
    ]

    res = [(PropVectorLs[ii]*D*PropVectorRs[ii])[1] for ii = 1:length(pairs)]
    return res
end

"""
    G_R(z::ComplexF64,
    pairs::Vector{Tuple{GrapheneCoord,GrapheneCoord}},
    s::GrapheneSystem)

The full real-space graphene Green's function in the presence of
defects as a function of complex energy `z`.

# Arguments
* `z`: complex energy
* `pairs`: pairs of [`GrapheneCoord`](@ref)'s for which `G_R` is calculated
* `s`: [`GrapheneSystem`](@ref) for which `G_R` is calculated
"""
function G_R(
    z::ComplexF64,
    pairs::Vector{Tuple{GrapheneCoord,GrapheneCoord}},
    s::GrapheneSystem,
)
    res =
        [graphene_propagator(p[1], p[2], z) for p in pairs] + δG_R(z, pairs, s)
    return res
end

"""
    δΓ(z::ComplexF64, s::GrapheneSystem)

The correction to the impurity Green's function due to the impurities'
interaction with graphene.

# Arguments
* `z`: complex energy
* `s`: [`GrapheneSystem`](@ref) for which `δΓ` is calculated
"""
function δΓ(z::ComplexF64, s::GrapheneSystem)
    if isempty(s.imps)
        error("No impurity states in the system")
    else
        Γ0 = 1 ./ (z .- s.imps) |> Diagonal
        prop_mat = propagator_matrix(z, s.scattering_atoms)
        Λ =
            prop_mat +
            prop_mat *
            s.Δ *
            inv(Diagonal(ones(length(s.scattering_atoms))) - prop_mat * s.Δ) *
            prop_mat
        res =
            Γ0 *
            adjoint(s.V) *
            Λ *
            inv(
                Diagonal(ones(length(s.scattering_atoms))) -
                s.V * Γ0 * adjoint(s.V) * Λ,
            ) *
            s.V *
            Γ0
        return res
    end
end

"""
    Γ(z::ComplexF64, s::GrapheneSystem)

The full impurity Green's function with the correction due to the impurities'
interaction with graphene.

# Arguments
* `z`: complex energy
* `s`: [`GrapheneSystem`](@ref) for which `Γ` is calculated
"""
function Γ(z::ComplexF64, s::GrapheneSystem)
    Γ0 = 1 ./ (z .- s.imps) |> Diagonal
    res = Γ0 + δΓ(z, s)
    return res
end
