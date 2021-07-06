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
    const nopert = Tuple{GrapheneCoord,GrapheneCoord,ComplexF64}[]

An empty array to be used in constructing the [`GrapheneSystem`](@ref) if
there are no direct perturbation terms.
"""
const nopert = Tuple{GrapheneCoord,GrapheneCoord,ComplexF64}[]

"""
    GrapheneSystem(
        μ::Float64,
        T::Float64,
        Δ::Array{ComplexF64,2},
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
    Δ::Array{ComplexF64,2}                  # Δ matrix
    V::Array{Float64,2}                     # V Matrix
    scattering_atoms::Vector{GrapheneCoord} # List of all perturbed atoms
    imps::Vector{Float64}                   # Impurity energies
end

"""
    mkGrapheneSystem(
        μ::Float64,
        T::Float64,
        imps::Vector{ImpurityState},
        pert::Vector{Tuple{GrapheneCoord,GrapheneCoord,ComplexF64}},
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
    pert::Vector{Tuple{GrapheneCoord,GrapheneCoord,ComplexF64}},
)
    # Create a coupling dictionary
    if isempty(pert)
        pert_dict = Dict{Tuple{GrapheneCoord,GrapheneCoord},ComplexF64}()
    else
        pert_dict =
            reduce(
                vcat,
                map(x -> [((x[1], x[2]), x[3]), ((x[2], x[1]), conj(x[3]))], pert),
            ) |> Dict
    end

    # Get the coordinates of all the directly-perturbed atoms
    perturbed_atoms = [x[1] for x in pert_dict |> keys |> collect] |> unique
    # Get the coordinates of atoms coupled to impurities
    coupled_atoms = [
        y[2]
        for
        y in reduce(
            vcat,
            [x.coupling for x in imps],
            init = Vector{Tuple{Float64,GrapheneCoord}}[],
        )
    ]
    # Combine the two sets of coordinates and keep only unique entries
    all_atoms = vcat(perturbed_atoms, coupled_atoms) |> unique |> sort
    # Assemble Δ
    Δ = [get!(pert_dict, (x, y), 0.0 + 0.0im) for x in all_atoms, y in all_atoms]
    # Assemble V
    V_array = [
        [
            sum([((atom == c[2]) * c[1]) for c in imp.coupling])::Float64 for atom in all_atoms
        ] for imp in imps
    ]
    if !isempty(V_array)
        V = reduce(hcat, V_array)
    else
        V = Array{GrapheneCoord}(undef, 0, 0)
    end
    imp_energies = [x.ϵ for x in imps]
    return GrapheneSystem(μ, T, Δ, V, all_atoms, imp_energies)
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
        u -> vec_pot(a1_loc[1] + u * disp[1], a1_loc[2] + u * disp[2]) .* disp |> sum,
        0,
        1,
    )[1]
    return (-π * res / UC_area)
end
