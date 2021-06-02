include("pristine_graphene.jl")

"""
    struct Coupling
        V::Float64              # Coupling to a graphene atom
        coord::GrapheneCoord    # Location of the graphene atom
    end

A structure describing the coupling `V` (in eV) between an impurity state and a
graphene atom at `coord`.
"""
struct Coupling
    V::Float64              # Coupling to a graphene atom
    coord::GrapheneCoord    # Location of the graphene atom
end

"""
    mutable struct ImpurityState
        ϵ::Float64                  # Impurity state energy
        coupling::Vector{Coupling}  # Coupling array for the impurity
    end

A structure describing an impurity state with energy `ϵ` (in eV) and containing
a list of all its couplings to graphene atoms.
"""
mutable struct ImpurityState
    ϵ::Float64                  # Impurity state energy
    coupling::Vector{Coupling}  # Coupling array for the impurity
end

"""
    new_impurity(ϵ::Float64)

A function to create a new impurity.

# Arguments
* `ϵ`: energy of the impurity state

# Output
* [`ImpurityState`](@ref) with an empty coupling list
"""
function new_impurity(ϵ)
    return ImpurityState(ϵ, Array{Coupling}(undef, 0))
end

"""
    add_coupling!(imp::ImpurityState, V::Float64, coord::GrapheneCoord)

A function to add a coupling to an existing [`ImpurityState`](@ref).

# Arguments
* `imp`: an existing [`ImpurityState`](@ref)
* `V`: coupling energy (in eV)
* `coord`: [`GrapheneCoord`](@ref) to which one wishes to couple the impurity

# Output
* [`ImpurityState`](@ref) with an added element in the coupling list
"""
function add_coupling!(imp::ImpurityState, V, coord)
    if (coord in map(x -> x.coord, imp.coupling))
        error("The impurity is already coupled to this atom")
    else
        imp.coupling = push!(imp.coupling, Coupling(V, coord))
    end
end

"""
    remove_coupling!(imp::ImpurityState, ind::Int)

A function to remove a coupling from an existing [`ImpurityState`](@ref).

# Arguments
* `imp`: an existing [`ImpurityState`](@ref)
* `ind`: the index of the coupling in the coupling list to be removed

# Output
* [`ImpurityState`](@ref) with the desired coupling removed
"""
function remove_coupling!(imp::ImpurityState, ind::Int)
    imp.coupling = deleteat!(imp.coupling, ind)
end

"""
    mutable struct GrapheneSystem
        μ::Float64                          # Chemical potential
        T::Float64                          # Temperature
        imps::Vector{ImpurityState}         # Impurity states in the system
        pert::Dict{Tuple{GrapheneCoord,GrapheneCoord},Float64}  # Direct perturbation
        Δ::Array{Float64,2}                 # Δ matrix
        V::Array{Float64,2}                 # V Matrix
        scattering_atoms::Vector{GrapheneCoord} # List of all perturbed atoms
    end

A structure describing the perturbed graphene system. Whenever an impurity or
a direct perturbation (coupling between graphene atoms) is added or removed,
the `Δ` and `V` matrices are updated, as is the list of the perturbed atoms
`scattering_atoms`.
"""
mutable struct GrapheneSystem
    μ::Float64                          # Chemical potential
    T::Float64                          # Temperature
    imps::Vector{ImpurityState}         # Impurity states in the system
    pert::Dict{Tuple{GrapheneCoord,GrapheneCoord},Float64}  # Δ perturbation
    Δ::Array{Float64,2}
    V::Array{Float64,2}
    scattering_atoms::Vector{GrapheneCoord}
end

"""
    new_graphene_system()

Create a new [`GrapheneSystem`](@ref).

# Arguments
* None

# Output
* [`GrapheneSystem`](@ref) with ``T = 0``, ``μ = 0``, and no defects.
"""
function new_graphene_system()
    return GrapheneSystem(
        0.0,
        0.0,
        Array{ImpurityState}(undef, 0),
        Dict{Tuple{GrapheneCoord,GrapheneCoord},Float64}(),
        Array{Float64}(undef, 0, 0),
        Array{Float64}(undef, 0, 0),
        Array{GrapheneCoord}(undef, 0),
    )
end

"""
    set_T!(s::GrapheneSystem, new_T::Float64)

Modify the temperature of an existing [`GrapheneSystem`](@ref)

# Arguments
* `s`: an existing [`GrapheneSystem`](@ref)
* `new_T`: new temperature

# Output
* [`GrapheneSystem`](@ref) with the new temperature `new_T`
"""
function set_T!(s::GrapheneSystem, new_T::Float64)
    s.T = new_T
end

"""
    set_μ!(s::GrapheneSystem, new_μ::Float64)

Modify the chemical potential of an existing [`GrapheneSystem`](@ref)

# Arguments
* `s`: an existing [`GrapheneSystem`](@ref)
* `new_μ`: new chemical potential

# Output
* [`GrapheneSystem`](@ref) with the new chemical potential `new_μ`
"""
function set_μ!(s::GrapheneSystem, new_μ::Float64)
    s.μ = new_μ
end

"""
    add_perturbation!(
        s::GrapheneSystem,
        a::GrapheneCoord,
        b::GrapheneCoord,
        Δ::Float64,
    )

Introduce a direct coupling between graphene atoms `a` and `b` or, if `a==b`,
a local potential

# Arguments
* `s`: an existing [`GrapheneSystem`](@ref)
* `a`: [`GrapheneCoord`](@ref) of the first atom
* `b`: [`GrapheneCoord`](@ref) of the second atom
* `Δ`: coupling between the atoms

# Output
* [`GrapheneSystem`](@ref) with the newly-added coupling
"""
function add_perturbation!(
    s::GrapheneSystem,
    a::GrapheneCoord,
    b::GrapheneCoord,
    Δ::Float64,
)
    c = s.pert
    c[(a, b)] = Δ
    c[(b, a)] = Δ
    s.pert = c
    scattering!(s)
end

"""
    remove_perturbation!(
        s::GrapheneSystem,
        a::GrapheneCoord,
        b::GrapheneCoord,
    )

Remove a coupling between graphene atoms `a` and `b` or, if `a==b`, a local potential

# Arguments
* `s`: an existing [`GrapheneSystem`](@ref)
* `a`: [`GrapheneCoord`](@ref) of the first atom
* `b`: [`GrapheneCoord`](@ref) of the second atom

# Output
* [`GrapheneSystem`](@ref) with the coupling removed
"""
function remove_perturbation!(
    s::GrapheneSystem,
    a::GrapheneCoord,
    b::GrapheneCoord,
)
    c = s.pert
    delete!(c, (a, b))
    delete!(c, (b, a))
    s.pert = c
    scattering!(s)
end

"""
    add_impurity!(s::GrapheneSystem, imp::ImpurityState)

A function to add an impurity to an existing [`GrapheneSystem`](@ref).

# Arguments
* `s`: an existing [`GrapheneSystem`](@ref)
* `imp`: the [`ImpurityState`](@ref) to be added

# Output
* [`GrapheneSystem`](@ref) with an added element in the impurity list
"""
function add_impurity!(s::GrapheneSystem, imp::ImpurityState)
    s.imps = push!(s.imps, imp)
    scattering!(s)
end

"""
    remove_impurity!(s::GrapheneSystem, ind::Int)

A function to remove an impurity from an existing [`GrapheneSystem`](@ref).

# Arguments
* `s`: an existing [`GrapheneSystem`](@ref)
* `ind`: the index of the impurity in the impurity list to be removed

# Output
* [`GrapheneSystem`](@ref) with the desired impurity removed
"""
function remove_impurity!(s::GrapheneSystem, ind::Int)
    s.imps = deleteat!(s.imps, ind)
    scattering!(s)
end

function scattering!(s::GrapheneSystem)
    # Get the coordinates of all the directly-perturbed atoms
    perturbed_atoms = map(x -> x[1], s.pert |> keys |> collect) |> unique
    # Get the coordinates of atoms coupled to impurities
    coupled_atoms = map(
        y -> y.coord,
        reduce(vcat, map(x -> x.coupling, s.imps), init = Coupling[]),
    )
    # Combine the two sets of coordinates and keep only unique entries
    all_atoms = vcat(perturbed_atoms, coupled_atoms) |> unique
    s.scattering_atoms = all_atoms
    # Assemble Δ
    all_atoms_M = repeat(all_atoms, 1, length(all_atoms))
    all_atoms_M_T = permutedims(all_atoms_M)
    Δ = map((x, y) -> get!(s.pert, (x, y), 0.0), all_atoms_M, all_atoms_M_T)
    # Assemble V
    V_array = map(
        imp -> map(
            atom -> sum(
                map(c -> ((atom == c.coord) * c.V), imp.coupling),
            )::Float64,
            all_atoms,
        ),
        s.imps,
    )
    s.Δ = Δ
    if !isempty(V_array)
        s.V = reduce(hcat, V_array)
    else
        s.V = Array{GrapheneCoord}(undef, 0, 0)
    end
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
    Γ0 = map(x -> z - x.ϵ, s.imps) |> Diagonal |> inv
    prop_mat = propagator_matrix(z, s.scattering_atoms)
    D =
        (s.Δ .+ s.V * Γ0 * transpose(s.V)) * inv(
            Diagonal(ones(length(s.scattering_atoms))) .-
            prop_mat * (s.Δ .+ s.V * Γ0 * transpose(s.V)),
        )
    PropVectorR = map(x -> graphene_propagator(x, a2, z), s.scattering_atoms)
    if a1 == a2
        PropVectorL = permutedims(PropVectorR)
    else
        PropVectorL =
            map(x -> propagator(a1, x, z), s.scattering_atoms) |> permutedims

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
    Γ0 = map(x -> z - x.ϵ, s.imps) |> Diagonal |> inv
    prop_mat = propagator_matrix(z, s.scattering_atoms)
    Λ =
        prop_mat +
        prop_mat *
        s.Δ *
        inv(Diagonal(ones(length(s.scattering_atoms))) - prop_mat * s.Δ) *
        prop_mat
    res =
        Γ0 *
        transpose(s.V) *
        Λ *
        inv(
            Diagonal(ones(length(s.scattering_atoms))) -
            s.V * Γ0 * transpose(s.V) * Λ,
        ) *
        s.V *
        Γ0
    return res
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
    Γ0 = map(x -> z - x.ϵ, s.imps) |> Diagonal |> inv
    res = Γ0 + δΓ(z, s)
    return res
end
