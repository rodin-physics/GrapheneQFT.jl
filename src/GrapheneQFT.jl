"""
    GrapheneQFT
``t = 2.8``eV
"""
module GrapheneQFT

using MLStyle
using Cubature
using QuadGK

export graphene_A,
    graphene_B,
    crystal_to_cartesian,
    graphene_neighbors,
    graphene_propagator,
    graphene_propagator_matrix,
    new_graphene_system,
    set_T!,
    set_μ!,
    add_perturbation!,
    remove_perturbation!,
    add_coupling!,
    remove_coupling!,
    add_impurity!,
    remove_impurity!

include("defects.jl")

end # module
