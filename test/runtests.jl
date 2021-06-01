using Test, GrapheneQFT

r1 = rand(-20:20, 20)
r2 = rand(-20:20, 20)
@test map((x, y) -> graphene_A(x, y), r1, r2) ==
      map((x, y) -> GrapheneQFT.GrapheneCoord(x, y, GrapheneQFT.A), r1, r2)
@test map((x, y) -> graphene_B(x, y), r1, r2) ==
      map((x, y) -> GrapheneQFT.GrapheneCoord(x, y, GrapheneQFT.B), r1, r2)
# @test graphene_A(5, 1) == GrapheneQFT.GrapheneCoord(5, 1, GrapheneQFT.A)
# out = plusTwo(3)
#
# @test out == 5
