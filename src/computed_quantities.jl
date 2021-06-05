include("defects.jl")
"""
    δG_R(z::ComplexF64, a1::GrapheneCoord, a2::GrapheneCoord, s::GrapheneSystem)

The correction to the real-space graphene Green's function in the presence of
defects as a function of complex energy `z`.

# Arguments
* `z`: complex energy
* `a1`: [`GrapheneCoord`](@ref) of the first atom
* `a2`: [`GrapheneCoord`](@ref) of the second atom
* `s`: [`GrapheneSystem`](@ref) for which `δG_R` is calculated
"""
function δG_R(z::ComplexF64, a1::GrapheneCoord, a2::GrapheneCoord, s::GrapheneSystem)
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
    PropVectorR = map(x -> graphene_propagator(x, a2, z), s.scattering_atoms)
    if a1 == a2
        PropVectorL = permutedims(PropVectorR)
    else
        PropVectorL =
            map(x -> graphene_propagator(a1, x, z), s.scattering_atoms) |>
            permutedims

    end
    res = (PropVectorL*D*PropVectorR)[1]
    return res
end

"""
    G_R(z::ComplexF64, a1::GrapheneCoord, a2::GrapheneCoord, s::GrapheneSystem)

The full real-space graphene Green's function in the presence of
defects as a function of complex energy `z`.

# Arguments
* `z`: complex energy
* `a1`: [`GrapheneCoord`](@ref) of the first atom
* `a2`: [`GrapheneCoord`](@ref) of the second atom
* `s`: [`GrapheneSystem`](@ref) for which `G_R` is calculated
"""
function G_R(z::ComplexF64, a1::GrapheneCoord, a2::GrapheneCoord, s::GrapheneSystem)
    res =
        (a1 == a2) *
        graphene_propagator(graphene_A(0, 0), graphene_A(0, 0), z) +
        δG_R(z, a1, a2, s)
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
