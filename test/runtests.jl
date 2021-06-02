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
# #
# my_system = new_graphene_system()
# imp1 = new_impurity(1.0)
# imp2 = new_impurity(1.0)
#
# add_coupling!(imp1, 2, graphene_A(0, 0))
# add_coupling!(imp2, 2, graphene_A(2, 2))
#
# add_impurity!(my_system, imp1)
# add_impurity!(my_system, imp2)
#
# add_perturbation!(my_system, random_atom(), random_atom(), 6.1)
#
# @time GrapheneQFT.propagator_matrix(1im, my_system.scattering_atoms)
# @time δG_R(1im, graphene_A(0, 0), graphene_A(0, 0), my_system)
# @time Γ(1im, my_system)[1, 1]
#
# @time GrapheneQFT.Ω(1im, 12, 4)
#
# # remove_perturbation!(my_system, random_atom(), random_atom())
@time quadgk(
      x -> δG_R(x * 1im, graphene_A(0, 0), graphene_A(0, 0), my_system)|>real,
      0,
      Inf,
      rtol = 1e-2,
)
@time δG_R(2 * 1im, graphene_A(1, 1), graphene_A(1, 1), my_system)
# @code_warntype δG_R(2 * 1im, graphene_A(1, 1), graphene_A(1, 1), my_system)
# @time GrapheneQFT.Ω(1im, 12, 4)
# @time δG_R(2 * 1im, graphene_A(1, 1), graphene_A(1, 1), my_system)
@time GrapheneQFT.Ω(1im, 12, 4)
