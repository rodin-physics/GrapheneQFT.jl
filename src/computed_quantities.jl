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
    imps_len = length(s.imps)
    scatter_len = length(s.scattering_atoms)
    # Calculate the scattering matrix D which is the same
    # for every pair of coordinates
    prop_mat = propagator_matrix(z, s.scattering_atoms)
    if imps_len == 0
        D = s.Δ * inv(Diagonal(ones(2*scatter_len)) .- prop_mat * s.Δ)
    else
        Γ0_init = 1 ./ (z .- s.imps) |> Diagonal
        Γ0 = spin_expand(Γ0_init)

        D =
            (s.Δ .+ s.J .+ s.V * Γ0 * adjoint(s.V)) * inv(
                Diagonal(ones(2*scatter_len)) .-
                prop_mat * (s.Δ .+ s.J .+ s.V * Γ0 * adjoint(s.V)),
            )
    end

    PropVectorRs =
        [[graphene_propagator(x, p[2], z) for x in s.scattering_atoms] for p in pairs]

    PVR = map(spin_expand, PropVectorRs)

    PropVectorLs = [
        pairs[idx][1] == pairs[idx][2] ? permutedims(PropVectorRs[idx]) :
            [graphene_propagator(pairs[idx][1], x, z) for x in s.scattering_atoms] |> permutedims
        for idx = 1:length(pairs)
    ]

    PVL = map(spin_expand, PropVectorLs)

    res = [(PVL[ii]*D*PVR[ii]) for ii = 1:length(pairs)]
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
    res = [spin_expand(graphene_propagator(p[1], p[2], z)) for p in pairs] + δG_R(z, pairs, s)
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
        Γ0 = 1 ./ (z .- s.imps) |> Diagonal |> spin_expand
        prop_mat = propagator_matrix(z, s.scattering_atoms)
        Λ =
            prop_mat +
            prop_mat *
            (s.Δ .+ s.J) *
            inv(Diagonal(ones(2*length(s.scattering_atoms))) - prop_mat * (s.Δ .+ s.J)) *
            prop_mat
        res =
            Γ0 *
            adjoint(s.V) *
            Λ *
            inv(Diagonal(ones(2*length(s.scattering_atoms))) - s.V * Γ0 * adjoint(s.V) * Λ) *
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
    Γ0 = 1 ./ (z .- s.imps) |> Diagonal |> spin_expand
    res = Γ0 + δΓ(z, s)
    return res
end

"""
    δρ_R_graphene(loc::GrapheneCoord, spin::Int64, s::GrapheneSystem)

The correction to charge density in graphene induced by impurities at a given `GrapheneCoord`. A `spin` value of 0 returns the total charge density correction, while 1 and -1 return the spin-up and spin-down charge density correction, respectively.

# Arguments
* `loc`: [`GrapheneCoord`](@ref) for which `δρ_R_graphene` is calculated
* `spin`: Value determining whether the spin-up, spin-down or total charge density correction is calculated
* `s`: [`GrapheneSystem`](@ref) for which `δρ_R_graphene` is calculated

"""

function δρ_R_graphene(loc::GrapheneCoord, spin::Int64, s::GrapheneSystem)
    spin ∈ Int64[1, 0, -1] || error("Invalid spin value")

    if spin == 0
        a = b = true
    else
        a = (spin == 1)
        b = !a
    end

    if s.T == 0
        res = quadgk(
            x -> real(G_R(s.μ + 1im * x, [(loc, loc)], s)[1]),
            0,
            Inf,
            rtol = 1e-2,
        )

        return ((res[1][1,1] * a + res[1][2,2] * b) / π)::Float64
    else
        error("Finite T given")
    end
end
