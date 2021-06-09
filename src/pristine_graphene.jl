# Nearest-neighbor hopping parameter
const NN_hopping = 2.8

# Graphene basis vectors
# Distance between neighboring carbon atoms
const graphene_lattice_constant = 2.46
const sublattice_shift = -1 / √(3) * graphene_lattice_constant

const graphene_d1 =
    [graphene_lattice_constant / 2, graphene_lattice_constant * √(3) / 2]

const graphene_d2 =
    [-graphene_lattice_constant / 2, graphene_lattice_constant * √(3) / 2]

#  Algebraic data type for graphene sublattices
@data Sublattice begin
    A
    B
end

"""
    GrapheneCoord(u::Int, v::Int, sublattice::Sublattice)

Lattice coordinate of a carbon atom, generated using [`graphene_A`](@ref) or
[`graphene_B`](@ref).

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
    graphene_A(u::Int, v::Int)

Create a [`GrapheneCoord`](@ref) for an atom belonging to sublattice A at the
unit cell (u, v).
"""
function graphene_A(u::Int, v::Int)
    return GrapheneCoord(u, v, A)
end

"""
    graphene_B(u::Int, v::Int)

Create a [`GrapheneCoord`](@ref) for an atom belonging to sublattice B at the
unit cell (u, v).
"""
function graphene_B(u::Int, v::Int)
    return GrapheneCoord(u, v, B)
end

## Propagator
# Integrals used in computing the propagator

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
        π,
        atol = 1e-16,
        rtol = 1e-4,
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
        π,
        atol = 1e-16,
        rtol = 1e-4,
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
        π,
        atol = 1e-16,
        rtol = 1e-4,
    ))[1])
end



# The propagator function picks out the correct element of the Ξ matrix based
# on the sublattices of the graphene coordinates.
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


# Given a list of [`GrapheneCoord`](@ref), this functiohn returns a Ξ(z) propagator
# matrix. The calculation is sped up using the fact that the matrix is symmetric.
function crystal_to_cartesian(coord::GrapheneCoord)
    u = coord.u
    v = coord.v
    x = graphene_d1[1] * u + graphene_d2[1] * v
    y = graphene_d1[2] * u + graphene_d2[2] * v

    return ([x, y + (coord.sublattice == B) * sublattice_shift, 0.0])
end


function propagator_matrix(z, Coords::Vector{GrapheneCoord})
    precomputed = Dict{Float64,ComplexF64}()
    len_coords = length(Coords)
    out = zeros(ComplexF64, len_coords, len_coords)
    for ii = 1:len_coords
        @inbounds for jj = ii:len_coords
            dist = norm(
                crystal_to_cartesian(Coords[ii]) -
                crystal_to_cartesian(Coords[jj]),
            )
            c = get(precomputed, dist, 0.0)
            if c == 0.0
                res = graphene_propagator(Coords[ii], Coords[jj], z)
                out[ii, jj] = res
                out[jj, ii] = res
                precomputed[dist] = res
            else
                out[ii, jj] = c
                out[jj, ii] = c
            end
        end
    end
    return out
end
