using Test, QuadGK, GrapheneQFT

function random_atom()
      u = rand(-20:20)
      v = rand(-20:20)
      return GrapheneQFT.GrapheneCoord(
            u,
            v,
            rand() < 0.5 ? GrapheneQFT.A : GrapheneQFT.B,
      )
end

r1 = rand(-20:20, 10)
r2 = rand(-20:20, 10)
@test map((x, y) -> graphene_A(x, y), r1, r2) ==
      map((x, y) -> GrapheneQFT.GrapheneCoord(x, y, GrapheneQFT.A), r1, r2)
@test map((x, y) -> graphene_B(x, y), r1, r2) ==
      map((x, y) -> GrapheneQFT.GrapheneCoord(x, y, GrapheneQFT.B), r1, r2)

a1 = graphene_A(0, 17)
a2 = graphene_A(1, 1)
a3 = graphene_B(3, -2)
a4 = graphene_A(2, 7)
a5 = graphene_A(0, 1)

ϵ1 = 1.7
ϵ2 = -0.2

V1 = 0.4
V2 = 2+0.9im
V3 = -2.1

imp1 = ImpurityState(ϵ1, [Coupling(V1, a4), Coupling(V2, a3)])
imp2 = ImpurityState(ϵ2, [Coupling(V3, a5)])

c1 = ComplexF64(3)
c2 = ComplexF64(5 + 1im)
c3 = ComplexF64(9 - 2im)

pert = new_perturbation()
pert = add_perturbation(pert, a1, a1, c1)
pert = add_perturbation(pert, a2, a1, c2)
pert = add_perturbation(pert, a3, a1, c3)

my_system = mk_GrapheneSystem(0.0, 0.0, [imp1, imp2], pert)

@test my_system.imps == [ϵ1, ϵ2]
@test my_system.scattering_atoms == [a5, a2, a4, a1, a3]

@test my_system.Δ == [
      0 0 0 0 0
      0 0 0 c2 0
      0 0 0 0 0
      0 conj(c2) 0 c1 conj(c3)
      0 0 0 c3 0
]
@test my_system.V == [
      0 V3
      0 0
      V1 0
      0 0
      V2 0
]
