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

a1 = random_atom()
a2 = random_atom()
a3 = random_atom()

pert = new_perturbation()
pert = add_perturbation(pert, graphene_A(0, 0), graphene_A(0, 0), ComplexF64(3))
my_system = mk_GrapheneSystem(0.0, 0.0, GrapheneQFT.ImpurityState[], pert)

# @code_warntype δG_R(1im, graphene_A(0,0), graphene_A(0,0), z)
