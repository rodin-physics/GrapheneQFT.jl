# """
#     crystal_to_cartesian(coord::GrapheneCoord)
#
# Convert a crystal coordinate to a cartesian one for plotting.
#
# # Arguments
# * `coord`: a [`GrapheneCoord`](@ref) that needs to be converted to a 3D [`Location`](@ref)
#
# # Output
# * [`Location`](@ref) of the carbon atom
# """
# function crystal_to_cartesian(coord::GrapheneCoord)
#     u = coord.u
#     v = coord.v
#     x = graphene_d1[1] * u + graphene_d2[1] * v
#     y = graphene_d1[2] * u + graphene_d2[2] * v
#
#     return (Location(x, y + (coord.sublattice == B) * shft, 0.0))
# end
#


# """
#     struct Location
#         x::Float64
#         y::Float64
#         z::Float64
#     end
#
# A structure describing a point in 3D space, with the lengths in Å
# """
# struct Location
#     x::Float64
#     y::Float64
#     z::Float64
# end
#
# # Graphene basis vectors
# const graphene_d1 = Location(
#     graphene_lattice_constant / 2,
#     graphene_lattice_constant * √(3) / 2,
#     0.0,
# )
#
# const graphene_d2 = Location(
#     -graphene_lattice_constant / 2,
#     graphene_lattice_constant * √(3) / 2,
#     0.0,
# )
#
# # Distance between neighboring carbon atoms
# const sublattice_shift = -1 / √(3) * graphene_lattice_constant
# All lengths are in Å
# const graphene_lattice_constant = 2.46
