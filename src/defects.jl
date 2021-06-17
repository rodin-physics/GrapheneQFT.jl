include("pristine_graphene.jl")

"""
    ImpurityState(ϵ::Float64, coupling::Vector{Tuple{Float64,GrapheneCoord}})

An impurity state of energy `ϵ` (in eV) coupled to the graphene system. The
tuples in the `coupling` field contain all the coupling energies (in eV) and the
 corresponding [`GrapheneCoord`](@ref)'s.
"""
struct ImpurityState
    ϵ::Float64
    coupling::Vector{Tuple{Float64,GrapheneCoord}}
end

"""
    const noimps = ImpurityState[]

An empty array to be used in constructing the [`GrapheneSystem`](@ref) if
there are no impurity states.
"""
const noimps = ImpurityState[]

"""
    const nopert = Tuple{GrapheneCoord,GrapheneCoord,Float64}[]

An empty array to be used in constructing the [`GrapheneSystem`](@ref) if
there are no direct perturbation terms.
"""
const nopert = Tuple{GrapheneCoord,GrapheneCoord,Float64}[]

"""
    GrapheneSystem(
        μ::Float64,
        T::Float64,
        Δ::Array{Float64,2},
        V::Array{Float64,2},
        scattering_atoms::Vector{GrapheneCoord},
        imps::Vector{Float64},
    )

A structure describing the perturbed graphene system.

See also [`mkGrapheneSystem`](@ref).
"""
struct GrapheneSystem
    μ::Float64                              # Chemical potential
    T::Float64                              # Temperature
    Δ::Array{Float64,2}                     # Δ matrix
    V::Array{Float64,2}                     # V Matrix
    scattering_atoms::Vector{GrapheneCoord} # List of all perturbed atoms
    imps::Vector{Float64}                   # Impurity energies
end

"""
    mkGrapheneSystem(
        μ::Float64,
        T::Float64,
        imps::Vector{ImpurityState},
        pert::Vector{Tuple{GrapheneCoord,GrapheneCoord,Float64}},
    )

Construct [`GrapheneSystem`](@ref).

# Arguments
* `μ`: chemical potential
* `T`: temperature
* `imps`: a list of [`ImpurityState`]'s
* `pert`: a list of 3-tuples describing the coupling between carbon atoms

When supplying `pert`, ensure that each coordinate pair appears only once as
repeated pairs with different couplings will cause eariler values to be
overwritten. The order of the coordinates does not matter. A tuple with the same
coordinate given twice creates and on-site potential for the coordinate.

The function constructs a [`GrapheneSystem`](@ref) with the `Δ` and `V` matrix
computed from `imps` and `pert`. In addition, a list of all
[`GrapheneCoord`](@ref) that are perturbed (`scattering_atoms` field in
[`GrapheneCoord`](@ref)) and a list of
impurity energies (`imps` field in [`GrapheneCoord`](@ref)) are included.
"""
function mkGrapheneSystem(
    μ::Float64,
    T::Float64,
    imps::Vector{ImpurityState},
    pert::Vector{Tuple{GrapheneCoord,GrapheneCoord,Float64}},
)
    # Create a coupling dictionary
    if isempty(pert)
        pert_dict = Dict{Tuple{GrapheneCoord,GrapheneCoord},Float64}()
    else
        pert_dict =
            reduce(
                vcat,
                map(x -> [((x[1], x[2]), x[3]), ((x[2], x[1]), x[3])], pert),
            ) |> Dict
    end

    # Get the coordinates of all the directly-perturbed atoms
    perturbed_atoms = map(x -> x[1], pert_dict |> keys |> collect) |> unique
    # Get the coordinates of atoms coupled to impurities
    coupled_atoms = map(
        y -> y[2],
        reduce(
            vcat,
            map(x -> x.coupling, imps),
            init = Vector{Tuple{Float64,GrapheneCoord}}[],
        ),
    )
    # Combine the two sets of coordinates and keep only unique entries
    all_atoms = vcat(perturbed_atoms, coupled_atoms) |> unique |> sort
    # Assemble Δ
    all_atoms_M = repeat(all_atoms, 1, length(all_atoms))
    all_atoms_M_T = permutedims(all_atoms_M)
    Δ = map((x, y) -> get!(pert_dict, (x, y), 0.0), all_atoms_M, all_atoms_M_T)
    # Assemble V
    V_array = map(
        imp -> map(
            atom -> sum(
                map(c -> ((atom == c[2]) * c[1]), imp.coupling),
            )::Float64,
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
# mkGrapheneSystem(0.0, 0.0, GrapheneQFT.noimps, GrapheneQFT.nopert)
#
# pert_dict =
#     Dict{Tuple{GrapheneQFT.GrapheneCoord,GrapheneQFT.GrapheneCoord},Float64}()
# pert_dict|>keys
#
# pert_dict_2 = [((graphene_A(0,0), graphene_B(0,0)),2.5)]|>Dict
#
# typeof(pert_dict) == typeof(pert_dict_2)
# typeof(pert_dict_2)
#
# typeof(pert_dict)
#
