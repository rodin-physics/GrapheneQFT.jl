"""
This package provides provides a set of functions to facilitate the field-theoretic
treatment of monolayer graphene using the tight-binding model. The Hamiltonian
employed by this package includes only the nearest-neighbor hopping term with
``t = 2.8``eV. The derivation of the formalism is available
[here](https://arxiv.org/pdf/2007.06984.pdf).
"""
module GrapheneQFT

using Kronecker
using LinearAlgebra
using MLStyle
using QuadGK

export A,
    B,
    SpinUp,
    SpinDown,
    GrapheneCoord,
    GrapheneState,
    crystal_to_cartesian,
    graphene_neighbors,
    graphene_multiple_neighbors,
    Defect,
    ImpurityState,
    LocalSpin,
    Hopping,
    GrapheneSystem,
    mkGrapheneSystem,
    δG_R,
    G_R,
    δΓ,
    Γ,
    δρ_R_graphene,
    peierls_phase,
    δF

include("computed_quantities.jl")

end # module
