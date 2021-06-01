include("pristine_graphene.jl")

# Coupling for impurity states
struct Coupling
    V::Float64              # Coupling to a graphene atom
    coord::GrapheneCoord    # Location of the graphene atom
end

# Impurity state defined by its on-site energy and its coupling to graphene
mutable struct ImpurityState
    ϵ::Float64                  # Impurity state energy
    coupling::Vector{Coupling}  # Coupling array for the impurity
end

# Initialize an impurity
function new_impurity(ϵ)
    return ImpurityState(ϵ, Array{Coupling}(undef, 0))
end

# Add coupling to an existing impurity
function add_coupling!(imp::ImpurityState, V, coord)
    if (coord in map(x -> x.coord, imp.coupling))
        error("The impurity is already coupled to this atom")
    else
        imp.coupling = push!(imp.coupling, Coupling(V, coord))
    end
end

# Remove coupling from an impurity
function remove_coupling!(imp::ImpurityState, ind::Int)
    imp.coupling = deleteat!(imp.coupling, ind)
end

mutable struct GrapheneSystem
    μ::Float64                          # Chemical potential
    T::Float64                          # Temperature
    imps::Vector{ImpurityState}         # Impurity states in the system
    pert::Dict{Tuple{GrapheneCoord,GrapheneCoord},Float64}  # Δ perturbation
    Δ::Array{Float64,2}
    V::Array{Float64,2}
    scattering_atoms::Vector{GrapheneCoord}
end

# Setting up the graphene system
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

# Changing the temperature of an initialized system
function set_T!(s::GrapheneSystem, new_T::Float64)
    s.T = new_T
end

# Changing the chemical potential of an initialized system
function set_μ!(s::GrapheneSystem, new_μ::Float64)
    s.μ = new_μ
end

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

function add_impurity!(s::GrapheneSystem, imp::ImpurityState)
    s.imps = push!(s.imps, imp)
    scattering!(s)
end

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
        s.V = []
    end
end

## Green's functions

function δG_R(z, a1::GrapheneCoord, a2::GrapheneCoord, s::GrapheneSystem)
    Γ0 = map(x -> z - x.ϵ, s.imps) |> Diagonal |> inv
    prop_mat = propagator_matrix(z, map(x -> x.coord, s.scattering_atoms))
    D =
        (s.Δ .+ s.V * Γ0 * transpose(s.V)) * inv(
            Diagonal(ones(length(s.scattering_atoms))) .-
            prop_mat * (s.Δ .+ s.V * Γ0 * transpose(s.V)),
        )
    PropVectorR = map(x -> propagator(x.coord, a2, z), s.scattering_atoms)
    if a1 == a2
        PropVectorL = permutedims(PropVectorR)
    else
        PropVectorL = map(x -> propagator(a1, x.coord, z), s.scattering_atoms)
    end
    res = (PropVectorL*D*PropVectorR)[1]
    return res
end

function G_R(z, a1::GrapheneCoord, a2::GrapheneCoord, s::GrapheneSystem)
    res =
        (a1 == a2) *
        graphene_propagator(graphene_A(0, 0), graphene_A(0, 0), z) +
        δG_R(z, a1, a2, s)
end

function δΓ(z, s::GrapheneSystem)
    Γ0 = map(x -> z - x.ϵ, s.imps) |> Diagonal |> inv
    prop_mat = propagator_matrix(z, map(x -> x.coord, s.scattering_atoms))
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

function Γ(z, s::GrapheneSystem)
    Γ0 = map(x -> z - x.ϵ, s.imps) |> Diagonal |> inv
    res = Γ0 + δΓ(z, s)
    return res
end
# imp1 = new_impurity(.3)
# imp2 = new_impurity(.7)
# add_coupling!(imp1,2.0,graphene_A(0,0))
# add_perturbation!(zz,graphene_A(0,0), graphene_B(2,1), 3.4)
# #
# add_coupling!(imp2, 1.4, graphene_B(5,5))
# # imp1
# zz = new_graphene_system()
# add_impurity!(zz,imp1)
# add_impurity!(zz,imp2)
# zz
# # Array{Float64}(undef, 0, 0)
# zz
# add_perturbation!(zz, graphene_A(0,0), graphene_B(5,5),7.0)
#
