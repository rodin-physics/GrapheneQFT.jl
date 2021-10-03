include("pristine_graphene.jl")

"""
    Defect

Abstract type for all defect types.
"""
abstract type Defect end

"""
    LocalSpin(x::Float64, y::Float64, z::Float64, coord::GrapheneCoord) <: Defect

A local spin with components `x`, `y`, and `z` located at `coord`. See
[`GrapheneCoord`](@ref).
"""
struct LocalSpin <: Defect
    x::Float64
    y::Float64
    z::Float64
    coord::GrapheneCoord
end

"""
    ImpurityState(ϵ::Float64, coupling::Vector{Tuple{Float64,GrapheneCoord}}) <: Defect

An impurity state of energy `ϵ` (in eV) coupled to the graphene system. The
tuples in the `coupling` field contain all the coupling energies (in eV) and the
 corresponding [`GrapheneCoord`](@ref)'s.
"""
struct ImpurityState <: Defect
    ϵ::Float64
    coupling::Vector{Tuple{Float64,GrapheneCoord}}
end

"""
    Hopping(c1::GrapheneCoord, c2::GrapheneCoord, Δ::ComplexF64) <: Defect

Hopping modification (in eV) between two [`GrapheneCoord`](@ref)'s. If `c1==c2`,
 this quantity corresponds to a local energy modification.
"""
struct Hopping <: Defect
    c1::GrapheneCoord
    c2::GrapheneCoord
    Δ::ComplexF64
end

"""
    GrapheneSystem(
        μ::Float64,
        T::Float64,
        Δ::Array{ComplexF64,2},
        V::Array{Float64,2},
        scattering_states::Vector{GrapheneState},
        imps::Vector{Float64},
    )

A structure describing the perturbed graphene system.

See also [`mkGrapheneSystem`](@ref).
"""
struct GrapheneSystem
    μ::Float64                                  # Chemical potential
    T::Float64                                  # Temperature
    Δ::Array{ComplexF64,2}                      # Δ matrix
    V::Array{Float64,2}                         # V Matrix
    scattering_states::Vector{GrapheneState}    # List of all perturbed states
    imps::Vector{Float64}                       # Impurity energies
end

"""
    mkGrapheneSystem(
        μ::Float64,
        T::Float64,
        defects::Vector{Defect},
    )

Construct [`GrapheneSystem`](@ref).

# Arguments
* `μ`: chemical potential
* `T`: temperature
* `defects`: a list of [`Defect`]'s

When supplying []`Hopping`](@ref) in `defects`, ensure that each coordinate pair
appears only once as repeated pairs with different couplings will cause eariler
values to be overwritten. The order of the coordinates does not matter.

The function constructs a [`GrapheneSystem`](@ref) with the `Δ` and `V` matrices.
 In addition, a list of all
[`GrapheneState`](@ref) that are perturbed (`scattering_states` field in
[`GrapheneCoord`](@ref)) and a list of
impurity energies (`imps` field in [`GrapheneCoord`](@ref)) are included.
"""
function mkGrapheneSystem(μ::Float64, T::Float64, defects::Vector{Defect})

    hops = filter(x -> typeof(x) == Hopping, defects)
    spins = filter(x -> typeof(x) == LocalSpin, defects)
    imps = filter(x -> typeof(x) == ImpurityState, defects)

    # Create a coupling dictionary
    if isempty(hops)
        hops_dict = Dict{Tuple{GrapheneState,GrapheneState},ComplexF64}()
    else
        hops_dict =
            reduce(
                vcat,
                [
                    [
                        (
                            (
                                GrapheneState(x.c1, SpinUp),
                                GrapheneState(x.c2, SpinUp),
                            ),
                            x.Δ,
                        ),
                        (
                            (
                                GrapheneState(x.c2, SpinUp),
                                GrapheneState(x.c1, SpinUp),
                            ),
                            conj(x.Δ),
                        ),
                        (
                            (
                                GrapheneState(x.c1, SpinDown),
                                GrapheneState(x.c2, SpinDown),
                            ),
                            x.Δ,
                        ),
                        (
                            (
                                GrapheneState(x.c2, SpinDown),
                                GrapheneState(x.c1, SpinDown),
                            ),
                            conj(x.Δ),
                        ),
                    ] for x in hops
                ],
            ) |> Dict
    end

    if isempty(spins)
        spins_dict = Dict{Tuple{GrapheneState,GrapheneState},ComplexF64}()
    else
        spins_dict =
            reduce(
                vcat,
                [
                    [
                        (
                            (
                                GrapheneState(s.coord, SpinUp),
                                GrapheneState(s.coord, SpinUp),
                            ),
                            s.z,
                        ),
                        (
                            (
                                GrapheneState(s.coord, SpinDown),
                                GrapheneState(s.coord, SpinDown),
                            ),
                            -s.z,
                        ),
                        (
                            (
                                GrapheneState(s.coord, SpinUp),
                                GrapheneState(s.coord, SpinDown),
                            ),
                            s.x - 1im * s.y,
                        ),
                        (
                            (
                                GrapheneState(s.coord, SpinDown),
                                GrapheneState(s.coord, SpinUp),
                            ),
                            s.x + 1im * s.y,
                        ),
                    ] for s in spins
                ],
            ) |> Dict
    end

    pert_dict = mergewith(+, hops_dict, spins_dict)

    # Get the coordinates of all the directly-perturbed atoms
    perturbed_atoms =
        [x[1].coord for x in pert_dict |> keys |> collect] |> unique
    # Get the coordinates of atoms coupled to impurities
    coupled_atoms = [
        y[2] for y in reduce(
            vcat,
            [x.coupling for x in imps],
            init = Vector{Tuple{Float64,GrapheneCoord}}[],
        )
    ]
    # Combine the two sets of coordinates and keep only unique entries
    all_atoms = vcat(perturbed_atoms, coupled_atoms) |> unique |> sort
    all_states = vcat(
        [GrapheneState(x, SpinUp) for x in all_atoms],
        [GrapheneState(x, SpinDown) for x in all_atoms],
    )

    # Assemble Δ
    Δ = [
        get!(pert_dict, (x, y), 0.0 + 0.0im) for x in all_states,
        y in all_states
    ]

    # Assemble V
    V_per_spin = [
        [
            (
                sum([
                    ((atom == c[2]) * c[1]) for c in imp.coupling
                ])::Float64
            ) for atom in all_atoms
        ] for imp in imps
    ]

    if !isempty(V_per_spin)
        V_per_spin = reduce(hcat, V_per_spin)
        V = [
            1 0
            0 1
        ] ⊗ V_per_spin |> collect
    else
        V = Array{Float64}(undef, 0, 0)
    end

    imp_energies = [x.ϵ for x in imps]
    return GrapheneSystem(μ, T, Δ, V, all_states, imp_energies)
end


"""
    peierls_phase(vec_pot, a1::GrapheneCoord, a2::GrapheneCoord)

Calculate the phase used in the Peierls substitution to include the effects of
the magnetic field.

To make the units work out better, the magnetic
field ``\\mathbf{B}(x, y, z) = \\Phi_0/\\mathcal{V} \\mathbf{f}(x, y, z)``,
where ``\\mathcal{V}`` is the area of the graphene unit cell in Å², ``\\Phi_0 = h / 2e``
is the magnetic flux quantum, and ``\\mathbf{f}(x, y, z) =
\\nabla \\times \\mathbf{g}(x,y,z)`` is a dimensionless vector function. Note that
``\\mathbf{f} = 1`` produces a field of about 40000 T.

For the vector potential, we have
``\\mathbf{A}(x, y, z) =\\Phi_0 / \\mathcal{V} \\mathbf{g}(x, y, z)``, where
``\\mathbf{g}(x, y, z)`` has the units of Å. Using the definition of the Peierls
phase, one gets

```math
\\phi = -\\frac{\\pi}{ \\Phi_0} \\int \\mathbf{A}\\cdot d\\mathbf{l} =
-\\frac{\\pi}{ \\mathcal{V}} \\int  \\mathbf{g}(x, y, z) \\cdot d\\mathbf{l}
\\rightarrow
-\\frac{\\pi}{ \\mathcal{V}} \\int  \\mathbf{g}_{xy}(x, y) \\cdot d\\mathbf{l}\\,.
```

The last step follows from the fact that the graphene system resides in the ``xy``
plane, so one needs to retain only the ``x`` and ``y`` componends of ``\\mathbf{g}``,
as denoted by the subscript ``xy``.

# Arguments
* `vec_pot(x,y)`: ``\\mathbf{g}_{xy}(x, y)`` in Å with `x` and `y` in Å.
* `a1`: [`GrapheneCoord`](@ref) of the "from" atom.
* `a2`: [`GrapheneCoord`](@ref) of the "to" atom.

`vec_pot(x,y)` needs to return a tuple corresponding to the vector potential in
    ``x`` and ``y`` directions.
"""
function peierls_phase(vec_pot, a1::GrapheneCoord, a2::GrapheneCoord)
    a1_loc = crystal_to_cartesian(a1)
    a2_loc = crystal_to_cartesian(a2)
    disp = a2_loc - a1_loc
    res = quadgk(
        u ->
            vec_pot(a1_loc[1] + u * disp[1], a1_loc[2] + u * disp[2]) .* disp |> sum,
        0,
        1,
    )[1]
    return (-π * res / UC_area)
end
