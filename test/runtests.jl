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

@test quadgk(
      x ->
            -G_R(
                  x + 1im * 1e-4,
                  graphene_A(0, 0),
                  graphene_A(0, 0),
                  new_graphene_system(),
            ) |> imag,
      -Inf,
      Inf,
)[1] ≈ π




#
#
#
# δG_R(1im, graphene_A(0, 0), graphene_A(0, 0), new_graphene_system())
# #
# new_graphene_system().scattering_atoms
# new_graphene_system().scattering_atoms


#
#
# zz =  new_graphene_system()
#
# zz.scattering_atoms
#
