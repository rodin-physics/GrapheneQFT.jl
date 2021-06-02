"""
This package provides provides a set of functions to facilitate the field-theoretic
treatment of monolayer graphene using the tight-binding model. The Hamiltonian
employed by this package includes only the nearest-neighbor hopping term with
``t = 2.8``eV. The derivation of the formalism is available
[here](https://arxiv.org/pdf/2007.06984.pdf).

# Example:

Pick a few graphene coordinates:
```
a1 = graphene_A(0, 17)
a2 = graphene_A(1, 1)
a3 = graphene_B(3, -2)
a4 = graphene_A(2, 7)
a5 = graphene_A(0, 1)
```

Next, initialize a couple of impurities and couple them to the coordinates:

```
ϵ1 = 1.7
ϵ2 = -0.2

V1 = 0.4
V2 = 2+0.9im
V3 = -2.1

imp1 = ImpurityState(ϵ1, [Coupling(V1, a4), Coupling(V2, a3)])
imp2 = ImpurityState(ϵ2, [Coupling(V3, a5)])
```

Create a few coupling integrals `c` and generate a coupling dictionary `pert`:
```
c1 = ComplexF64(3)
c2 = ComplexF64(5 + 1im)
c3 = ComplexF64(9 - 2im)

pert = new_perturbation()
```

Add the `c`'s to the dictionary and initialize the system:
```
pert = add_perturbation(pert, a1, a1, c1)
pert = add_perturbation(pert, a2, a1, c2)
pert = add_perturbation(pert, a3, a1, c3)

μ = 0.2
T = 0.0

my_system = mk_GrapheneSystem(μ, T, [imp1, imp2], pert)
```

`my_system` can now be used to calculate the relevant Green's functions
([`δG_R`](@ref), [`G_R`](@ref), [`δΓ`](@ref), [`Γ`](@ref)), from which
quantities of interest can be computed.
"""
module GrapheneQFT

using MLStyle
using Cubature
using QuadGK
using LinearAlgebra

export graphene_A,
    graphene_B,
    Coupling,
    ImpurityState,
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
