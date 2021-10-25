const NN_hopping = 2.8
const graphene_lattice_constant = 2.46
const sublattice_shift = -1 / √(3) * graphene_lattice_constant
const graphene_d1 =
    [graphene_lattice_constant / 2, graphene_lattice_constant * √(3) / 2]
const graphene_d2 =
    [-graphene_lattice_constant / 2, graphene_lattice_constant * √(3) / 2]
const UC_area =
    graphene_d1[1] * graphene_d2[2] - graphene_d1[2] * graphene_d2[1] |> abs

@data Sublattice begin
    A
    B
end

function Base.:!(s::Sublattice)
    @match s begin
        A => B
        B => A
    end
end

function Base.isless(s1::Sublattice, s2::Sublattice)
    @match (s1, s2) begin
        (A, B) => true
        _ => false
    end
    end

@data Spin begin
    SpinUp
    SpinDown
end

function Base.:!(s::Spin)
    @match s begin
        SpinUp => SpinDown
        SpinDown => SpinUp
    end
end

function Base.isless(s1::Spin, s2::Spin)
    @match (s1, s2) begin
        (SpinDown, SpinUp) => true
        _ => false
    end
end

"""
    GrapheneCoord(u::Int, v::Int, sublattice::Sublattice)

Lattice coordinate of a carbon atom.

Each coordinate contains the sublattice index `A` or `B`, as well as
the integer coefficients of the two basis vectors
``d\\times(\\pm 1 \\hat{x} + \\sqrt{3}\\hat{y}) / 2``
(`u` for `+`, `v` for `-`), with ``d = 2.46``Å as the lattice constant.
"""
struct GrapheneCoord
    u::Int
    v::Int
    sublattice::Sublattice
end

function Base.show(io::IO, a::GrapheneCoord)
    print(io, "|", a.u, ", ", a.v, ", ", a.sublattice, "⟩")
end

function Base.isless(a1::GrapheneCoord, a2::GrapheneCoord)
    if a1.sublattice != a2.sublattice
        return isless(a1.sublattice, a2.sublattice)
    else
        if a1.v == a2.v
            return isless(a1.u, a2.u)
        else
            return isless(a1.v, a2.v)
        end
    end
end

"""
    GrapheneState(coord::GrapheneCoord, spin::Spin)

Quantum state of an electron in graphene, denoted by 
``|u, v, L\\rangle\\otimes |\\sigma\\rangle`` in the drivation.

The state is given by the [`GrapheneCoord`](@ref) of the orbital, as well as
the electronic spin, which can take values `SpinUp` and `SpinDown`.
"""
struct GrapheneState
    coord::GrapheneCoord
    spin::Spin
end

function Base.show(io::IO, a::GrapheneState)
    print(io, a.coord, "⊗|", a.spin, "⟩")
end

function Base.isless(a1::GrapheneState, a2::GrapheneState)
    if a1.spin != a2.spin
        return isless(a1.spin, a2.spin)
    else
        return isless(a1.coord, a2.coord)
    end
end

"""
    crystal_to_cartesian(coord::GrapheneCoord)

Convert a [`GrapheneCoord`](@ref) to a cartesian point with lengths in Å, where
`GrapheneCoord(0, 0, A)` is at the origin.
"""
function crystal_to_cartesian(coord::GrapheneCoord)
    u = coord.u
    v = coord.v
    x = graphene_d1[1] * u + graphene_d2[1] * v
    y = graphene_d1[2] * u + graphene_d2[2] * v

    return ((x, y + (coord.sublattice == B) * sublattice_shift))
end

"""
    graphene_neighbors(atom::GrapheneCoord)

Determine the nearest neighbors of an `atom` and return a vector of the
corresponding [`GrapheneCoord`](@ref)'s.
"""
function graphene_neighbors(atom::GrapheneCoord)
    u = atom.u
    v = atom.v
    s = !atom.sublattice
    return [
        GrapheneCoord(u, v, s),
        GrapheneCoord(u + (-1)^(s == A), v, s),
        GrapheneCoord(u, v + (-1)^(s == A), s),
    ]
end

"""
    graphene_multiple_neighbors(atom::GrapheneCoord, n::Int)

Return an array of [`GrapheneCoord`](@ref)'s obtained by iteratively running
[`graphene_neighbors`](@ref) `n` times, applying it only to the
newly-added [`GrapheneCoord`](@ref)'s from the past iteration. The entries in
the result are unique and sorted by their distance from `atom`.

# Arguments
* `atom`: [`GrapheneCoord`](@ref) from which the iteration begins
* `n`: number of iterations of [`graphene_neighbors`](@ref).
"""
function graphene_multiple_neighbors(atom::GrapheneCoord, idx::Int)
    res = [atom]
    atoms = [atom]
    for ii = 1:idx
        atoms =
            filter(x -> x ∉ res, mapreduce(graphene_neighbors, vcat, atoms)) |>
            unique
        res = vcat(res, atoms)
    end
    res = sort(
        res,
        by=a ->
            norm(crystal_to_cartesian(atom) - crystal_to_cartesian(a)),
    )
    return res
end

## Propagator
# When computing Ω, occasionally the integrand becomes small enough to give NaN
# The integrand helper functions help to catch these instances
@inline function Ω_Integrand(z, u, v, x::Float64)
    W = ((z / NN_hopping)^2 - 1.0) / (4.0 * cos(x)) - cos(x) |> Complex
    res =
        exp(1.0im * (u - v) * x) / cos(x) *
        ((W - √(W - 1) * √(W + 1))^abs.(u + v)) / (√(W - 1) * √(W + 1))
    return (isnan(res) ? 0.0 + 0.0im : res)
end

@inline function Ω(z, u, v)
    return ((quadgk(
        x -> 2 * Ω_Integrand(z, u, v, x) / (8.0 * π * NN_hopping^2),
        0.0,
        π / 3,
        2 * π / 3,
        π,
        atol=1e-16,
        rtol=1e-4,
    ))[1])
end

@inline function Ωp_Integrand(z, u, v, x::Float64)
    W = ((z / NN_hopping)^2 - 1.0) / (4.0 * cos(x)) - cos(x) |> Complex
    res =
        2 *
        exp(1.0im * (u - v) * x) *
        ((W - √(W - 1) * √(W + 1))^abs.(u + v + 1)) / (√(W - 1) * √(W + 1))
    return (isnan(res) ? 0.0 + 0.0im : res)
end

@inline function Ωp(z, u, v)
    return ((quadgk(
        x -> 2 * Ωp_Integrand(z, u, v, x) / (8.0 * π * NN_hopping^2),
        0.0,
        π / 3,
        2 * π / 3,
        π,
        atol=1e-16,
        rtol=1e-4,
    ))[1])
end

@inline function Ωn_Integrand(z, u, v, x::Float64)
    W = ((z / NN_hopping)^2 - 1.0) / (4.0 * cos(x)) - cos(x) |> Complex
    res =
        2 *
        exp(1.0im * (u - v) * x) *
        ((W - √(W - 1) * √(W + 1))^abs.(u + v - 1)) / (√(W - 1) * √(W + 1))
    return (isnan(res) ? 0.0 + 0.0im : res)
end

@inline function Ωn(z, u, v)
    return ((quadgk(
        x -> 2 * Ωn_Integrand(z, u, v, x) / (8.0 * π * NN_hopping^2),
        0.0,
        π / 3,
        2 * π / 3,
        π,
        atol=1e-16,
        rtol=1e-4,
    ))[1])
end

# The propagator function picks out the correct element of the Ξ matrix based
# on the sublattices of the graphene coordinates.
function graphene_propagator(a_l::GrapheneState, a_m::GrapheneState, z)
    if a_l.spin != a_m.spin
        return 0.0
    else
        u = (a_l.coord).u - (a_m.coord).u
        v = (a_l.coord).v - (a_m.coord).v
        @match ((a_l.coord).sublattice, (a_m.coord).sublattice) begin
            (A, A) => z * Ω(z, u, v)
            (B, B) => z * Ω(z, u, v)
            (A, B) => -NN_hopping * (Ω(z, u, v) + Ωp(z, u, v))
            (B, A) => -NN_hopping * (Ω(z, u, v) + Ωn(z, u, v))
        end
    end
end

        function propagator_matrix(z, States::Vector{GrapheneState})
    precomputed = Dict{Tuple{Int,Int,GrapheneQFT.Sublattice,GrapheneQFT.Sublattice,GrapheneQFT.Spin,GrapheneQFT.Spin,},ComplexF64,}()

    len_states = length(States)
    out = zeros(ComplexF64, len_states, len_states)
    for ii = 1:len_states
        @inbounds for jj = ii:len_states
            key_ = (
                (States[ii].coord).u - (States[jj].coord).u,
                (States[ii].coord).v - (States[jj].coord).v,
                (States[ii].coord).sublattice,
                (States[jj].coord).sublattice,
                States[ii].spin,
                States[jj].spin,
            )
            c = get(precomputed, key_, NaN)
            if isnan(c)
                res = graphene_propagator(States[ii], States[jj], z)
                out[ii, jj] = res
                out[jj, ii] = res
                precomputed[key_] = res
            else
                out[ii, jj] = c
                out[jj, ii] = c
            end
        end
    end
    return out
end
