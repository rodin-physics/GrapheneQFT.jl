include("../defects.jl")

# Integrand used to calculate the local density correction in graphene
function δρ_Graphene_Integrand(
    loc::GrapheneCoord,
    z,
    imps::Vector{ImpurityState},
    mod_atoms::Vector{PerturbedAtom},
)
    # Scattering atoms
    atoms = scattering_atoms(imps, mod_atoms)
    D = scattering_matrix(z, imps, atoms)

    # Coordinates of the neighbors. Needed due to the orbital overlap
    neighbors_ = neighbors(loc)
    # Self-propagator
    PropVectorR = map(x -> propagator(x.coord, loc, z), atoms)
    res = (transpose(PropVectorR)*D*PropVectorR)[1]

    return (res)

end

# Local density function in graphene. For T = 0, we integrate along the
# imaginary axis. Otherwise, we perform a contour integration to enclose the
# Matsubara frequencies.
function δρ_Graphene(loc, s::GrapheneSystem)
    μ = s.μ
    T = s.T
    imps = s.imps
    mod_atoms = s.mod_atoms
    if T == 0
        res = quadgk(
            x -> real(δρ_Graphene_Integrand(loc, μ + 1im * x, imps, mod_atoms)),
            0,
            Inf,
            rtol = ν,
        )
        return (res[1] / π)::Float64
    else
        res = quadgk(
            x ->
                -imag(
                    δρ_Graphene_Integrand(
                        loc,
                        x + μ + 1im * η,
                        imps,
                        mod_atoms,
                    ),
                ) * nF(x, T),
            -Inf,
            Inf,
            rtol = ν,
        )
        return (res[1] / π)
    end
end

# Integrand used to calculate the graphene-induced density in impurity states
function δρ_Impurity_Integrand(
    n,
    z,
    imps::Vector{ImpurityState},
    mod_atoms::Vector{PerturbedAtom},
)
    atoms = scattering_atoms(imps, mod_atoms)
    n_atoms = length(atoms)
    Δ_ = Δ(atoms)
    V_ = V(atoms, imps)

    prop_mat = propagator_matrix(z, map(x -> x.coord, atoms))

    Λ =
        prop_mat .+
        prop_mat *
        Δ_ *
        inv(Diagonal(ones(n_atoms, n_atoms)) .- prop_mat * Δ_) *
        prop_mat

    Γ0 = map(x -> z - x.ϵ, imps) |> Diagonal |> inv
    res =
        Γ0 *
        transpose(V_) *
        Λ *
        inv(Diagonal(ones(n_atoms, n_atoms)) .- V_ * Γ0 * transpose(V_) * Λ) *
        V_ *
        Γ0

    return (res[n]::ComplexF64)
end

# Graphene-induced density in impurity states. For T = 0, we integrate along the
# imaginary axis. Otherwise, we perform a contour integration to enclose the
# Matsubara frequencies.
function δρ_Impurity(n, s::GrapheneSystem)
    μ = s.μ
    T = s.T
    imps = s.imps
    mod_atoms = s.mod_atoms
    if T == 0
        res = quadgk(
            x -> real(δρ_Impurity_Integrand(n, μ + 1im * x, imps, mod_atoms)),
            0,
            Inf,
            rtol = ν,
        )
        return (res[1] / π)
    else
        res = quadgk(
            x ->
                -imag(
                    δρ_Impurity_Integrand(n, x + μ + 1im * η, imps, mod_atoms),
                ) * nF(x, T),
            -Inf,
            Inf,
            rtol = ν,
        )
        return (res[1] / π)
    end
end

# Total impurity density is the sum of the isolated impurity filling factor
# plus a correction δρ
function ρ_Impurity(n, s::GrapheneSystem)
    en = map(x -> x.ϵ, s.imps)
    return δρ_Impurity(n, s) + nF(en[n] - s.μ, s.T)
end

# Total impurity density for all the impurities
function ρ_All(s::GrapheneSystem)
    res = map(n -> ρ_Impurity(n, s), 1:length(s.imps))
end

function spectral_graphene(ω, loc, s::GrapheneSystem)
    pristine_spectral = -(propagator(loc, loc, ω + 1im * η)) / π |> imag
    correction_spectral =
        -δρ_Graphene_Integrand(loc, ω + 1im * η, s.imps, s.mod_atoms) / π |>
        imag
    return (pristine_spectral + correction_spectral)
end

function spectral_impurity(ω, n, s::GrapheneSystem)
    Γ0_inv = map(x -> ω + 1im * η - x.ϵ, s.imps)
    pristine_spectral = -1 / (Γ0_inv[n]) / π |> imag
    correction_spectral =
        -δρ_Impurity_Integrand(n, ω + 1im * η, s.imps, s.mod_atoms) / π |> imag
    return (pristine_spectral + correction_spectral)
end
