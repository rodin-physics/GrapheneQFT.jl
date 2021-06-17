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
V2 = 2.0
V3 = -7.1

imp1 = ImpurityState(ϵ1, [(V1, a4), (V2, a3)])
imp2 = ImpurityState(ϵ2, [(V3, a5)])
```

Create a few coupling integrals `c` and use them to create some perturbations
```
c1 = 3.0
c2 = 5.0
c3 = -9.0

p1 = (a1, a1, c1)
p2 = (a2, a1, c2)
p3 = (a3, a1, c3)
```

Finally, initialize the system:
```
μ = 0.2
T = 0.0

my_system = mkGrapheneSystem(μ, T, [imp1, imp2], [p1, p2, p3])
```

`my_system` can now be used to calculate the relevant Green's functions
([`δG_R`](@ref), [`G_R`](@ref), [`δΓ`](@ref), [`Γ`](@ref)), from which
quantities of interest can be computed.
"""
module GrapheneQFT

using Cubature
using Interpolations
using LinearAlgebra
using MLStyle
using QuadGK

export GrapheneCoord,
    graphene_A,
    graphene_B,
    graphene_neighbors,
    graphene_multiple_neighbors,
    crystal_to_cartesian,
    Coupling,
    ImpurityState,
    GrapheneSystem,
    mkGrapheneSystem,
    δG_R,
    G_R,
    δΓ,
    Γ,
    Ψ_pz,
    coulomb_potential_pz,
    coulomb_potential_pz_interp,
    coulomb_energy_pz_pz

include("computed_quantities.jl")
include("orbitals.jl")

end # module
