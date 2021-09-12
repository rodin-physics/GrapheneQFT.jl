# Nearest-neighbor hopping parameter
const NN_hopping = 2.8
const P = 0.065

# Graphene basis vectors
# Distance between neighboring carbon atoms
const graphene_lattice_constant = 2.46
const sublattice_shift = -1 / √(3) * graphene_lattice_constant

const graphene_d1 = [graphene_lattice_constant / 2, graphene_lattice_constant * √(3) / 2]
const graphene_d2 = [-graphene_lattice_constant / 2, graphene_lattice_constant * √(3) / 2]

const UC_area = graphene_d1[1] * graphene_d2[2] - graphene_d1[2] * graphene_d2[1] |> abs

#  Algebraic data type for graphene sublattices
@data Sublattice begin
    A
    B
end

# Define a not function for the two sublattices
function Base.:!(s::Sublattice)
    @match s begin
        A => B
        B => A
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

# Define a comparison function between the coordinates to help with testing
function Base.isless(a1::GrapheneCoord, a2::GrapheneCoord)
    if a1.sublattice == A && a2.sublattice == B
        return true
    elseif a1.sublattice == B && a2.sublattice == A
        return false
    else
        if a1.v == a2.v
            return isless(a1.u, a2.u)
        else
            return isless(a1.v, a2.v)
        end
    end
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
        atoms = filter(x -> x ∉ res, mapreduce(graphene_neighbors, vcat, atoms)) |> unique
        res = vcat(res, atoms)
    end
    res = sort(res, by = a -> norm(crystal_to_cartesian(atom) - crystal_to_cartesian(a)))
    return res
end


# Helper function to expand dimensions to include spin effects
function spin_expand(array)
    array_init = [array[ii,jj]*I(2) for ii in 1:size(array,1), jj in 1:size(array,2)]
    final = hvcat(size(array_init,1), array_init...)
    return size(final,1) == size(final, 2) ? final : permutedims(final)

end

## Propagator
# Integrals used in computing the propagator

# When computing Ω, occasionally the integrand becomes small enough to give NaN
# The integrand helper functions help to catch these instances
@inline function Ω_Integrand(z, u, v, x::Float64)
    NNwOverlap = NN_hopping + z * P
    W = ((z / NNwOverlap)^2 - 1.0) / (4.0 * cos(x)) - cos(x) |> Complex
    res =
        exp(1.0im * (u - v) * x) / cos(x) * ((W - √(W - 1) * √(W + 1))^abs.(u + v)) /
        (√(W - 1) * √(W + 1))
    return (isnan(res) ? 0.0 + 0.0im : res)
end

@inline function Ω(z, u, v)
    NNwOverlap = NN_hopping + z * P
    return ((quadgk(
        x -> Ω_Integrand(z, u, v, x) / (8.0 * π * NNwOverlap^2),
        0.0,
        π / 3,
        2 * π / 3,
        π,
        atol = 1e-16,
        rtol = 1e-4,
    ))[1])
end

@inline function Ωp_Integrand(z, u, v, x::Float64)
    NNwOverlap = NN_hopping + z * P
    W = ((z / NNwOverlap)^2 - 1.0) / (4.0 * cos(x)) - cos(x) |> Complex
    res =
        2 * exp(1.0im * (u - v) * x) * ((W - √(W - 1) * √(W + 1))^abs.(u + v + 1)) /
        (√(W - 1) * √(W + 1))
    return (isnan(res) ? 0.0 + 0.0im : res)
end

@inline function Ωp(z, u, v)
    NNwOverlap = NN_hopping + z * P
    return ((quadgk(
        x -> Ωp_Integrand(z, u, v, x) / (8.0 * π * NNwOverlap^2),
        0.0,
        π / 3,
        2 * π / 3,
        π,
        atol = 1e-16,
        rtol = 1e-4,
    ))[1])
end

@inline function Ωn_Integrand(z, u, v, x::Float64)
    NNwOverlap = NN_hopping + z * P
    W = ((z / NNwOverlap)^2 - 1.0) / (4.0 * cos(x)) - cos(x) |> Complex
    res =
        2 * exp(1.0im * (u - v) * x) * ((W - √(W - 1) * √(W + 1))^abs.(u + v - 1)) /
        (√(W - 1) * √(W + 1))
    return (isnan(res) ? 0.0 + 0.0im : res)
end

@inline function Ωn(z, u, v)
    NNwOverlap = NN_hopping + z * P
    return ((quadgk(
        x -> Ωn_Integrand(z, u, v, x) / (8.0 * π * NNwOverlap^2),
        0.0,
        π / 3,
        2 * π / 3,
        π,
        atol = 1e-16,
        rtol = 1e-4,
    ))[1])
end

# The propagator function picks out the correct element of the Ξ matrix based
# on the sublattices of the graphene coordinates.
function graphene_propagator(a_l::GrapheneCoord, a_m::GrapheneCoord, z)
    NNwOverlap = NN_hopping + z * P
    u = a_l.u - a_m.u
    v = a_l.v - a_m.v
    @match [a_l.sublattice, a_m.sublattice] begin
        [A, A] => z * Ω(z, u, v)
        [B, B] => z * Ω(z, u, v)
        [A, B] => -NNwOverlap * (Ω(z, u, v) + Ωp(z, u, v))
        [B, A] => -NNwOverlap * (Ω(z, u, v) + Ωn(z, u, v))
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

    return ([x, y + (coord.sublattice == B) * sublattice_shift])
end

# Given a list of [`GrapheneCoord`](@ref), this functiohn returns a Ξ(z) propagator
# matrix. The calculation is sped up using the fact that the matrix is symmetric.
function propagator_matrix(z, Coords::Vector{GrapheneCoord})
    precomputed =
        Dict{Tuple{Int,Int,GrapheneQFT.Sublattice,GrapheneQFT.Sublattice},ComplexF64}()

    len_coords = length(Coords)
    out = zeros(ComplexF64, len_coords, len_coords)
    for ii = 1:len_coords
        @inbounds for jj = ii:len_coords
            key_ = (
                Coords[ii].u - Coords[jj].u,
                Coords[ii].v - Coords[jj].v,
                Coords[ii].sublattice,
                Coords[jj].sublattice,
            )
            c = get(precomputed, key_, 0.0)
            if c == 0.0
                res = graphene_propagator(Coords[ii], Coords[jj], z)
                out[ii, jj] = res
                out[jj, ii] = res
                precomputed[key_] = res
            else
                out[ii, jj] = c
                out[jj, ii] = c
            end
        end
    end

    return spin_expand(out)
end
