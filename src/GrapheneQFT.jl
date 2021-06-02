"""
This package provides provides a set of functions to facilitate the field-theoretic
treatment of monolayer graphene using the tight-binding model. The Hamiltonian
employed by this package includes only the nearest-neighbor hopping term with
``t = 2.8``eV. The derivation of the formalism is available
[here](https://arxiv.org/pdf/2007.06984.pdf).
"""
module GrapheneQFT

using MLStyle
using Cubature
using QuadGK
using LinearAlgebra

export graphene_A,
    graphene_B,
    new_impurity,
    new_graphene_system,
    new_perturbation,
    add_perturbation,
    mk_GrapheneSystem,
    δG_R,
    G_R,
    δΓ,
    Γ

include("defects.jl")

end # module
