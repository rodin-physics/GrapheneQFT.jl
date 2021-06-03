include("pristine_graphene.jl")

"""
    struct Coupling
        V::ComplexF64           # Coupling to a graphene atom
        coord::GrapheneCoord    # Location of the graphene atom
    end

A structure describing the coupling `V` (in eV) between an impurity state and a
graphene atom at `coord`.
"""
struct Coupling
    V::ComplexF64           # Coupling to a graphene atom
    coord::GrapheneCoord    # Location of the graphene atom
end

"""
    struct ImpurityState
        ϵ::Float64                  # Impurity state energy
        coupling::Vector{Coupling}  # Coupling array for the impurity
    end

A structure describing an impurity state with energy `ϵ` (in eV) and containing
a list of all its couplings to graphene atoms.
"""
struct ImpurityState
    ϵ::Float64                  # Impurity state energy
    coupling::Vector{Coupling}  # Coupling array for the impurity
end

"""
    struct Perturbation
        pert::Dict{Tuple{GrapheneCoord,GrapheneCoord},ComplexF64}
    end

A structure wrapping a dictionary describing the coupling between graphene atoms
and local potential. See [`new_perturbation()`](@ref) and
[`add_perturbation`](@ref) for more details.
"""
struct Perturbation
    pert::Dict{Tuple{GrapheneCoord,GrapheneCoord},ComplexF64}  # Direct perturbation
end

"""
    new_perturbation()

A helper initialization function , giving an empty dictionary wrapped in the
    [`Perturbation`](@ref) structure. The keys for the dictionary are tuples of two
    [`GrapheneCoord`](@ref), while the values are `ComplexF64`.

# Arguments:
* None

# Output:
* [`Perturbation`](@ref) wrapping an empty `Dict{Tuple{GrapheneCoord,GrapheneCoord},ComplexF64}`
"""
function new_perturbation()
    return Perturbation(Dict{Tuple{GrapheneCoord,GrapheneCoord},ComplexF64}())
end

"""
    function add_perturbation(
        p::Perturbation,
        a::GrapheneCoord,
        b::GrapheneCoord,
        Δ::ComplexF64,
    )

A function that adds coupling between atoms `a` and `b` or, if `a == b`, an
    on-site potential to an existing [`Perturbation`](@ref) structure. When
    adding a coupling between two atoms, two entries are added to the dictionary:
    one `(a,b)=>Δ` and another one `(b,a)=>conj(Δ)`. Passing the same pair more
    than once overwrites the previous entries.

# Arguments
* `p`: [`Perturbation`](@ref) to which an additional element is to be added
* `a`: [`GrapheneCoord`](@ref) of the first atom
* `b`: [`GrapheneCoord`](@ref) of the second atom
* `Δ`: coupling between the atoms

# Output
* [`Perturbation`](@ref) with the new coupling added
"""
function add_perturbation(
    p::Perturbation,
    a::GrapheneCoord,
    b::GrapheneCoord,
    Δ::ComplexF64,
)
    c = p.pert
    c[(a, b)] = Δ
    c[(b, a)] = conj(Δ)
    return Perturbation(c)
end

"""
    struct GrapheneSystem
        μ::Float64                              # Chemical potential
        T::Float64                              # Temperature
        Δ::Array{ComplexF64,2}                  # Δ matrix
        V::Array{ComplexF64,2}                  # V Matrix
        scattering_atoms::Vector{GrapheneCoord} # List of all perturbed atoms
        imps::Vector{Float64}                   # Impurity energies
    end

A structure describing the perturbed graphene system. See
[`mk_GrapheneSystem`](@ref) for details.
"""
struct GrapheneSystem
    μ::Float64                              # Chemical potential
    T::Float64                              # Temperature
    Δ::Array{ComplexF64,2}                  # Δ matrix
    V::Array{ComplexF64,2}                  # V Matrix
    scattering_atoms::Vector{GrapheneCoord} # List of all perturbed atoms
    imps::Vector{Float64}                   # Impurity energies
end

"""
    function mk_GrapheneSystem(
        μ::Float64,
        T::Float64,
        imps::Vector{ImpurityState},
        p::Perturbation,
    )

A function for constructing [`GrapheneSystem`](@ref).

# Arguments
* `μ`: chemical potential
* `T`: Temperature
* `imps`: a list of [`ImpurityState`]'s
* `p`: [`Perturbation`](@ref) struct containing the carbon-carbon coupling

# Output
* [`GrapheneSystem`](@ref) with the Δ and V matrix computed from `imps` and `p`.
In addition, a list of all [`GrapheneCoord`](@ref) that are perturbed and a
lise of impurity energies are included.
"""
function mk_GrapheneSystem(
    μ::Float64,
    T::Float64,
    imps::Vector{ImpurityState},
    p::Perturbation,
)
    # Get the coordinates of all the directly-perturbed atoms
    perturbed_atoms = map(x -> x[1], p.pert |> keys |> collect) |> unique
    # Get the coordinates of atoms coupled to impurities
    coupled_atoms = map(
        y -> y.coord,
        reduce(vcat, map(x -> x.coupling, imps), init = Coupling[]),
    )
    # Combine the two sets of coordinates and keep only unique entries
    all_atoms = vcat(perturbed_atoms, coupled_atoms) |> unique |> sort
    # Assemble Δ
    all_atoms_M = repeat(all_atoms, 1, length(all_atoms))
    all_atoms_M_T = permutedims(all_atoms_M)
    Δ = map((x, y) -> get!(p.pert, (x, y), 0.0), all_atoms_M, all_atoms_M_T)
    # Assemble V
    V_array = map(
        imp -> map(
            atom -> sum(
                map(c -> ((atom == c.coord) * c.V), imp.coupling),
            )::ComplexF64,
            all_atoms,
        ),
        imps,
    )
    if !isempty(V_array)
        V = reduce(hcat, V_array)
    else
        V = Array{GrapheneCoord}(undef, 0, 0)
    end
    imp_energies = map(x -> x.ϵ, imps)
    return GrapheneSystem(μ, T, Δ, V, all_atoms, imp_energies)
end

## Green's functions

"""
    δG_R(z, a1::GrapheneCoord, a2::GrapheneCoord, s::GrapheneSystem)

The correction to the real-space graphene Green's function in the presence of
defects as a function of complex energy `z`.

# Arguments
* `z`: complex energy
* `a1`: [`GrapheneCoord`](@ref) of the first atom
* `a2`: [`GrapheneCoord`](@ref) of the second atom
* `s`: [`GrapheneSystem`](@ref) for which `δG_R` is calculated

# Output
* `ComplexF64`
"""
function δG_R(z, a1::GrapheneCoord, a2::GrapheneCoord, s::GrapheneSystem)
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
    G_R(z, a1::GrapheneCoord, a2::GrapheneCoord, s::GrapheneSystem)

The full real-space graphene Green's function in the presence of
defects as a function of complex energy `z`.

# Arguments
* `z`: complex energy
* `a1`: [`GrapheneCoord`](@ref) of the first atom
* `a2`: [`GrapheneCoord`](@ref) of the second atom
* `s`: [`GrapheneSystem`](@ref) for which `G_R` is calculated

# Output
* `ComplexF64`
"""
function G_R(z, a1::GrapheneCoord, a2::GrapheneCoord, s::GrapheneSystem)
    res =
        (a1 == a2) *
        graphene_propagator(graphene_A(0, 0), graphene_A(0, 0), z) +
        δG_R(z, a1, a2, s)
end

"""
    δΓ(z, s::GrapheneSystem)

The correction to the impurity Green's function due to the impurities' interaction
with graphene.

# Arguments
* `z`: complex energy
* `s`: [`GrapheneSystem`](@ref) for which `δΓ_R` is calculated

# Output
* `Matrix{ComplexF64}`
"""
function δΓ(z, s::GrapheneSystem)
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
    Γ(z, s::GrapheneSystem)

The full impurity Green's function with the correction due to the impurities'
interaction with graphene.

# Arguments
* `z`: complex energy
* `s`: [`GrapheneSystem`](@ref) for which `Γ_R` is calculated

# Output
* `Matrix{ComplexF64}`
"""
function Γ(z, s::GrapheneSystem)
    Γ0 = 1 ./ (z .- s.imps) |> Diagonal
    res = Γ0 + δΓ(z, s)
    return res
end
