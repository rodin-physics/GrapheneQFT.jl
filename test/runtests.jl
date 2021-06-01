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
            -GrapheneQFT.graphene_propagator(
                  graphene_A(0, 0),
                  graphene_A(0, 0),
                  x + 1im * 1e-4,
            ) |> imag,
      -Inf,
      Inf,
)[1] ≈ π


# r_array_1 = map(x -> random_atom(), 1:20)
# r_array_2 = map(x -> random_atom(), 1:20)
# r_complex = 30 .* (rand(20) .- 1 / 2) + 30im .* (rand(20) .- 1 / 2)
#
# @test map(
#       (x, y, z) -> graphene_propagator(x, y, z),
#       r_array_1,
#       r_array_2,
#       r_complex,
# ) == map(
#       (x, y, z) -> graphene_propagator(x, y, z),
#       r_array_2,
#       r_array_1,
#       r_complex,
# )
