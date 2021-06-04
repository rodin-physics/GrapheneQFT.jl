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
V2 = 2.0
V3 = -7.1

imp1 = ImpurityState(ϵ1, [Coupling(V1, a4), Coupling(V2, a3)])
imp2 = ImpurityState(ϵ2, [Coupling(V3, a5)])

c1 = 3.0
c2 = 5.0
c3 = -9.0

p1 = (a1, a1, c1)
p2 = (a2, a1, c2)
p3 = (a3, a1, c3)


my_system = mkGrapheneSystem(0.0, 0.0, [imp1, imp2], [p1, p2, p3])

@test my_system.imps == [ϵ1, ϵ2]
@test my_system.scattering_atoms == [a5, a2, a4, a1, a3]

@test my_system.Δ == [
      0 0 0 0 0
      0 0 0 c2 0
      0 0 0 0 0
      0 c2 0 c1 c3
      0 0 0 c3 0
]
@test my_system.V == [
      0 V3
      0 0
      V1 0
      0 0
      V2 0
]

rand_num = ((rand() - 1 / 2) + 1im * (rand() - 1 / 2)) * 20

@test map(
      (x, y) -> G_R(rand_num, x, y, my_system),
      repeat(my_system.scattering_atoms, 1, length(my_system.scattering_atoms)),
      permutedims(
            repeat(
                  my_system.scattering_atoms,
                  1,
                  length(my_system.scattering_atoms),
            ),
      ),
) |> real ≈
      map(
      (x, y) -> G_R(conj(rand_num), x, y, my_system),
      repeat(my_system.scattering_atoms, 1, length(my_system.scattering_atoms)),
      permutedims(
            repeat(
                  my_system.scattering_atoms,
                  1,
                  length(my_system.scattering_atoms),
            ),
      ),
) |> real

@test map(
      (x, y) -> G_R(rand_num, x, y, my_system),
      repeat(my_system.scattering_atoms, 1, length(my_system.scattering_atoms)),
      permutedims(
            repeat(
                  my_system.scattering_atoms,
                  1,
                  length(my_system.scattering_atoms),
            ),
      ),
) |> imag ≈
      -map(
      (x, y) -> G_R(conj(rand_num), x, y, my_system),
      repeat(my_system.scattering_atoms, 1, length(my_system.scattering_atoms)),
      permutedims(
            repeat(
                  my_system.scattering_atoms,
                  1,
                  length(my_system.scattering_atoms),
            ),
      ),
) |> imag
