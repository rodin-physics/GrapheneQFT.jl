## Graphene parameters
# All lengths are in Å
const graphene_lattice_constant = 2.46
# Nearest-neighbor hopping parameter
const NN_hopping = 2.8

"""
    Location(x::Float64, y::Float64, z::Float64)

A structure describing a point in 3D space, with the lengths in Å
"""
struct Location
    x::Float64
    y::Float64
    z::Float64
end

# Graphene basis vectors
const graphene_d1 = Location(
    graphene_lattice_constant / 2,
    graphene_lattice_constant * √(3) / 2,
    0.0,
)

const graphene_d2 = Location(
    -graphene_lattice_constant / 2,
    graphene_lattice_constant * √(3) / 2,
    0.0,
)

# Distance between neighboring carbon atoms
const sublattice_shift = -1 / √(3) * graphene_lattice_constant

#  Algebraic data type for graphene sublattices
@data Sublattice begin
    A
    B
end

"""
Lattice coordinate of a carbon atom, generated using [`graphene_A`](@ref) or
[`graphene_B`](@ref). Each coordinate contains the sublattice index, as well as
the integer coefficients of the two basis vectors ``d\\times(\\pm 1 \\hat{x} +
\\sqrt{3}\\hat{y}) / 2``, where ``d = 2.46``Å is the graphene lattice constant.
"""
struct GrapheneCoord
    u::Int
    v::Int
    sublattice::Sublattice
end

"""
    crystal_to_cartesian(coord::GrapheneCoord)

Convert a crystal coordinate to a cartesian one.

# Arguments
* `coord`: a [`GrapheneCoord`](@ref) that needs to be converted to a 3D [`Location`](@ref)

# Output
* [`Location`](@ref) of the carbon atom
"""
function crystal_to_cartesian(coord::GrapheneCoord)
    u = coord.u
    v = coord.v
    x = graphene_d1[1] * u + graphene_d2[1] * v
    y = graphene_d1[2] * u + graphene_d2[2] * v

    return (Location(x, y + (coord.sublattice == B) * shft, 0.0))
end

"""
    graphene_A(u::Int, v::Int)

Create a [`GrapheneCoord`](@ref) for an atom belonging to sublattice A at the unit cell (u, v)

# Arguments
* `u`: coefficient of basis vector ``d\\times(1 \\hat{x} + \\sqrt{3}\\hat{y}) / 2``
* `v`: coefficient of basis vector ``d\\times(-1 \\hat{x} + \\sqrt{3}\\hat{y}) / 2``

# Output
* [`GrapheneCoord`](@ref) of the carbon atom
"""
function graphene_A(u::Int, v::Int)
    return GrapheneCoord(u, v, A)
end

"""
    graphene_B(u::Int, v::Int)

Create a [`GrapheneCoord`](@ref) for an atom belonging to sublattice B at the unit cell (u, v)

# Arguments
* `u`: coefficient of basis vector ``d\\times(1 \\hat{x} + \\sqrt{3}\\hat{y}) / 2``
* `v`: coefficient of basis vector ``d\\times(-1 \\hat{x} + \\sqrt{3}\\hat{y}) / 2``

# Output
* [`GrapheneCoord`](@ref) of the carbon atom
"""
function graphene_B(u::Int, v::Int)
    return GrapheneCoord(u, v, B)
end

"""
    graphene_neighbors(atom::GrapheneCoord)

Determine the nearest neighbors of an atom.

# Arguments
* `atom`: [`GrapheneCoord`](@ref) whose nearest neighbors are to be determined

# Output
* A Vector of [`GrapheneCoord`](@ref) containing the neighbors of `atom`
"""
function graphene_neighbors(atom::GrapheneCoord)
    u = atom.u
    v = atom.v
    if atom.sublattice == A
        return [
            graphene_B(u, v)
            graphene_B(u + 1, v)
            graphene_B(u, v + 1)
        ]
    elseif atom.sublattice == B
        return [
            graphene_A(u, v)
            graphene_A(u - 1, v)
            graphene_A(u, v - 1)
        ]
    else
        error("Illegal sublattice parameter")
    end
end

## Propagator
# Integrals used in computing the propagator
@inline function Ω_Integrand(z, u, v, x::Float64)
    t = NN_hopping
    W = ((z / t)^2 - 1.0) / (4.0 * cos(x)) - cos(x)
    return (
        exp(1.0im * (u - v) * x) / cos(x) *
        ((W - √(W - 1) * √(W + 1))^abs.(u + v)) / (√(W - 1) * √(W + 1))
    )
end

@inline function Ω(z, u, v)
    t = NN_hopping
    return ((quadgk(
        x -> Ω_Integrand(z, u, v, x) / (8.0 * π * t^2),
        0.0,
        2.0 * π,
    ))[1])
end

@inline function Ωp_Integrand(z, u, v, x::Float64)
    t = NN_hopping
    W = ((z / t)^2 - 1.0) / (4.0 * cos(x)) - cos(x)
    return (
        2 *
        exp(1.0im * (u - v) * x) *
        ((W - √(W - 1) * √(W + 1))^abs.(u + v + 1)) / (√(W - 1) * √(W + 1))
    )
end

@inline function Ωp(z, u, v)
    t = NN_hopping
    return ((quadgk(
        x -> Ωp_Integrand(z, u, v, x) / (8.0 * π * t^2),
        0.0,
        2.0 * π,
    ))[1])
end

@inline function Ωn_Integrand(z, u, v, x::Float64)
    t = NN_hopping
    W = ((z / t)^2 - 1.0) / (4.0 * cos(x)) - cos(x)
    return (
        2 *
        exp(1.0im * (u - v) * x) *
        ((W - √(W - 1) * √(W + 1))^abs.(u + v - 1)) / (√(W - 1) * √(W + 1))
    )
end

@inline function Ωn(z, u, v)
    t = NN_hopping
    return ((quadgk(
        x -> Ωn_Integrand(z, u, v, x) / (8.0 * π * t^2),
        0.0,
        2.0 * π,
    ))[1])
end


"""
    graphene_propagator(atom1::GrapheneCoord, atom2::GrapheneCoord, z)

The propagator function picks out the correct element of the Ξ matrix based
on the sublattices of the graphene coordinates.
"""
function graphene_propagator(a_l::GrapheneCoord, a_m::GrapheneCoord, z)
    t = NN_hopping
    u = a_l.u - a_m.u
    v = a_l.v - a_m.v
    if a_l.sublattice == a_m.sublattice
        return (z * Ω(z, u, v))
    elseif ([a_l.sublattice, a_m.sublattice] == [A, B])
        return (-t * (Ω(z, u, v) + Ωp(z, u, v)))
    elseif ([a_l.sublattice, a_m.sublattice] == [B, A])
        return (-t * (Ω(z, u, v) + Ωn(z, u, v)))
    else
        error("Illegal sublattice parameter")
    end
end

"""
    graphene_propagator_matrix(z, coords::Vector{GrapheneCoord})

Given a list of [`GrapheneCoord`](@ref), this functiohn returns a Ξ(z) propagator
matrix. The calculation is sped up using the fact that the matrix is symmetric.
"""
function propagator_matrix(z, Coords::Vector{GrapheneCoord})
    len_coords = length(Coords)
    out = zeros(ComplexF64, len_coords, len_coords)
    for ii = 1:len_coords
        @inbounds for jj = ii:len_coords
            out[ii, jj] = propagator(Coords[ii], Coords[jj], z)
            out[jj, ii] = out[ii, jj]
        end
    end
    return out
end