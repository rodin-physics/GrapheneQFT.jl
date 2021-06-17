using Test, QuadGK, GrapheneQFT, Cubature
## Testing pristine graphene
a1 = graphene_A(0, 17)
a2 = graphene_A(1, 1)
a3 = graphene_B(3, -2)
a4 = graphene_A(2, 7)
a5 = graphene_A(0, 1)

@test a2 < a1
@test a3 > a1
@test a2 > a5
@test a1 > a5

@test sort(graphene_neighbors(a4)) ==
      [graphene_B(2, 7), graphene_B(3, 7), graphene_B(2, 8)]

@test sort(graphene_neighbors(a3)) ==
      [graphene_A(3, -3), graphene_A(2, -2), graphene_A(3, -2)]

@test crystal_to_cartesian(a1) == [
      GrapheneQFT.graphene_d1[1] * 0 + GrapheneQFT.graphene_d2[1] * 17,
      GrapheneQFT.graphene_d1[2] * 0 + GrapheneQFT.graphene_d2[2] * 17,
]
@test crystal_to_cartesian(a3) == [
      GrapheneQFT.graphene_d1[1] * 3 + GrapheneQFT.graphene_d2[1] * (-2),
      GrapheneQFT.graphene_d1[2] * 3 +
      GrapheneQFT.graphene_d2[2] * (-2) +
      GrapheneQFT.sublattice_shift,
]
##

##

##
r1 = rand(-20:20, 10)
r2 = rand(-20:20, 10)
@test map((x, y) -> graphene_A(x, y), r1, r2) ==
      map((x, y) -> GrapheneQFT.GrapheneCoord(x, y, GrapheneQFT.A), r1, r2)
@test map((x, y) -> graphene_B(x, y), r1, r2) ==
      map((x, y) -> GrapheneQFT.GrapheneCoord(x, y, GrapheneQFT.B), r1, r2)

ϵ1 = 1.7
ϵ2 = -0.2

V1 = 0.4
V2 = 2.0
V3 = -7.1

imp1 = ImpurityState(ϵ1, [(V1, a4), (V2, a3)])
imp2 = ImpurityState(ϵ2, [(V3, a5)])

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
scattering_pairs =
      [
            (x, y) for x in my_system.scattering_atoms,
            y in my_system.scattering_atoms
      ] |> vec

@test G_R(rand_num, scattering_pairs, my_system) |> real ≈
      G_R(conj.(rand_num), scattering_pairs, my_system) |> real

@test G_R(rand_num, scattering_pairs, my_system) |> imag ≈
      -G_R(conj.(rand_num), scattering_pairs, my_system) |> imag

@test abs(
      hcubature(
            r -> 2 * π * Ψ_pz(r[1], r[2]) .^ 2 * r[1]^2 * sin(r[2]),
            [0, 0],
            [40, π],
      )[1] - 1,
) < 1e-6

R_rand = 40 * rand()
τ_rand = π * rand()
exact_coulomb = coulomb_potential_pz(R_rand, τ_rand)[1]
@test abs(coulomb_potential_pz_interp(R_rand, τ_rand) - exact_coulomb) /
      exact_coulomb < 1e-2
