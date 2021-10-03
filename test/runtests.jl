using Test, QuadGK, GrapheneQFT, LinearAlgebra

## Variables used for testing
a1 = GrapheneCoord(0, 17, A)
a2 = GrapheneCoord(1, 1, A)
a3 = GrapheneCoord(3, -2, B)
a4 = GrapheneCoord(2, 7, A)
a5 = GrapheneCoord(0, 1, A)

s1 = GrapheneState(a1, SpinUp)
s2 = GrapheneState(a2, SpinDown)
s3 = GrapheneState(a3, SpinUp)
s4 = GrapheneState(a4, SpinDown)
s5 = GrapheneState(a5, SpinUp)

ϵ1 = 1.7
ϵ2 = -0.2

V1 = 0.4
V2 = 2.0
V3 = -7.1

jx1 = 1.0
jy1 = 0.3
jz1 = 1.0

imp1 = ImpurityState(ϵ1, [(V1, a4), (V2, a3)])
imp2 = ImpurityState(ϵ2, [(V3, a5)])

sp1 = LocalSpin(jx1, jy1, jz1, a2)

c1 = 3.0 + 0im
c2 = 5.0 + 1im
c3 = -9.0 + 0im

h1 = Hopping(a1, a1, c1)
h2 = Hopping(a2, a1, c2)
h3 = Hopping(a3, a1, c3)
## Tests
@testset "Basic definitions" begin
    @test map(
        x -> GrapheneCoord(x.u, x.v, !x.sublattice),
        [a1, a2, a3, a4, a5],
    ) == [
        GrapheneCoord(0, 17, B),
        GrapheneCoord(1, 1, B),
        GrapheneCoord(3, -2, A),
        GrapheneCoord(2, 7, B),
        GrapheneCoord(0, 1, B),
    ]

    @test map(x -> GrapheneState(x.coord, !x.spin), [s1, s2, s3, s4, s5]) == [
        GrapheneState(a1, SpinDown),
        GrapheneState(a2, SpinUp),
        GrapheneState(a3, SpinDown),
        GrapheneState(a4, SpinUp),
        GrapheneState(a5, SpinDown),
    ]

    @test isless(a2, a1) == true
    @test isless(a3, a1) == false
    @test isless(a2, a5) == false
    @test isless(a1, a5) == false

    @test isless(s2, s1) == true
    @test isless(s3, s1) == false
    @test isless(s2, s5) == true
    @test isless(s1, s5) == false

    @test repr(a1) == "|0, 17, A⟩"
    @test repr(s3) == "|3, -2, B⟩⊗|SpinUp⟩"

end

@testset "Helper functions" begin
    @test sort(graphene_neighbors(a4)) == [
        GrapheneCoord(2, 7, B),
        GrapheneCoord(3, 7, B),
        GrapheneCoord(2, 8, B),
    ]

    @test sort(graphene_neighbors(a3)) == [
        GrapheneCoord(3, -3, A),
        GrapheneCoord(2, -2, A),
        GrapheneCoord(3, -2, A),
    ]
    #
    # @test vcat(a2, sort(graphene_neighbors(a2))) ==
    #       sort(graphene_multiple_neighbors(a2, 1))


    @test (
        crystal_to_cartesian(a1) .≈ [
            GrapheneQFT.graphene_d1[1] * 0 + GrapheneQFT.graphene_d2[1] * 17,
            GrapheneQFT.graphene_d1[2] * 0 + GrapheneQFT.graphene_d2[2] * 17,
        ]
    ) |> prod
    @test (
        crystal_to_cartesian(a3) .≈ [
            GrapheneQFT.graphene_d1[1] * 3 + GrapheneQFT.graphene_d2[1] * (-2),
            GrapheneQFT.graphene_d1[2] * 3 +
            GrapheneQFT.graphene_d2[2] * (-2) +
            GrapheneQFT.sublattice_shift,
        ]
    ) |> prod
end

@testset "No Impurity States" begin
    rand_num = ((rand() - 1 / 2) + 1im * (rand() - 1 / 2)) * 5
    my_system = mkGrapheneSystem(0.0, 0.0, [sp1, h1, h2, h3])
    @test my_system.T == 0.0
    @test my_system.μ == 0.0
    @test my_system.imps == Vector{Float64}[]
    @test my_system.scattering_states == vcat(
        [GrapheneState(x, SpinUp) for x in sort([a1, a2, a3])],
        [GrapheneState(x, SpinDown) for x in sort([a1, a2, a3])],
    )
    @test my_system.V == Array{Float64,2}(undef, 0, 0)
    @test my_system.Δ == [
        jz1 c2 0 (jx1-1im*jy1) 0 0
        conj(c2) c1 conj(c3) 0 0 0
        0 c3 0 0 0 0
        (jx1+1im*jy1) 0 0 (-jz1) c2 0
        0 0 0 conj(c2) c1 conj(c3)
        0 0 0 0 c3 0
    ]

    scattering_pair = [(s4, s5)]
    scattering_pair_rev = [(s5, s4)]

    @test G_R(rand_num, scattering_pair, my_system) |> real ≈
          G_R(conj.(rand_num), scattering_pair_rev, my_system) |> real

    @test G_R(rand_num, scattering_pair, my_system) |> imag ≈
          -G_R(conj.(rand_num), scattering_pair_rev, my_system) |> imag

    @test δG_R(rand_num, [(s1, s2), (s3, s4)], my_system) == G_R(rand_num, [(s1, s2), (s3, s4)], my_system)

end

@testset "Only Impurities" begin
    rand_num = ((rand() - 1 / 2) + 1im * (rand() - 1 / 2)) * 5
    my_system = mkGrapheneSystem(0.0, 0.0, Vector{Defect}([imp1, imp2]))
    @test my_system.T == 0.0
    @test my_system.μ == 0.0
    @test my_system.imps == [ϵ1, ϵ2]
    @test my_system.scattering_states == vcat(
        [GrapheneState(x, SpinUp) for x in sort([a3, a4, a5])],
        [GrapheneState(x, SpinDown) for x in sort([a3, a4, a5])],
    )
    @test my_system.V == [
        0 V3 0 0
        V1 0 0 0
        V2 0 0 0
        0 0 0 V3
        0 0 V1 0
        0 0 V2 0
    ]

    @test my_system.Δ == [
        0 0 0 0 0 0
        0 0 0 0 0 0
        0 0 0 0 0 0
        0 0 0 0 0 0
        0 0 0 0 0 0
        0 0 0 0 0 0
    ]

    scattering_pair = [(s1, s5)]
    scattering_pair_rev = [(s5, s1)]

    @test G_R(rand_num, scattering_pair, my_system) |> real ≈
          G_R(conj.(rand_num), scattering_pair_rev, my_system) |> real

    @test G_R(rand_num, scattering_pair, my_system) |> imag ≈
          -G_R(conj.(rand_num), scattering_pair_rev, my_system) |> imag

    @test G_R(rand_num, [(s1, s2)], my_system)[1] ≈ 0

    @test δΓ(rand_num, my_system) |> Diagonal ≈ conj.(δΓ(conj.(rand_num), my_system) |> Diagonal)

    @test Γ(rand_num, my_system) |> Diagonal ≈ conj.(Γ(conj.(rand_num), my_system) |> Diagonal)

    @test δρ_R_graphene(GrapheneState(a1, SpinUp), my_system) == δρ_R_graphene(GrapheneState(a1, SpinDown), my_system)

end

@testset "Impurities and Perturbation" begin
    rand_num = ((rand() - 1 / 2) + 1im * (rand() - 1 / 2)) * 5
    my_system = mkGrapheneSystem(0.0, 0.0, [sp1, h1, h2, h3, imp1, imp2])
    @test my_system.T == 0.0
    @test my_system.μ == 0.0
    @test my_system.imps == [ϵ1, ϵ2]
    @test my_system.scattering_states == vcat(
        [GrapheneState(x, SpinUp) for x in sort([a1, a2, a3, a4, a5])],
        [GrapheneState(x, SpinDown) for x in sort([a1, a2, a3, a4, a5])],
    )

    @test my_system.V == [
        0 V3 0 0
        0 0 0 0
        V1 0 0 0
        0 0 0 0
        V2 0 0 0
        0 0 0 V3
        0 0 0 0
        0 0 V1 0
        0 0 0 0
        0 0 V2 0
    ]

    @test my_system.Δ == [
        0 0 0 0 0 0 0 0 0 0
        0 jz1 0 c2 0 0 (jx1-1im*jy1) 0 0 0
        0 0 0 0 0 0 0 0 0 0
        0 conj(c2) 0 c1 conj(c3) 0 0 0 0 0
        0 0 0 c3 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 0 0
        0 (jx1+1im*jy1) 0 0 0 0 (-jz1) 0 c2 0
        0 0 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 conj(c2) 0 c1 conj(c3)
        0 0 0 0 0 0 0 0 c3 0
    ]

    scattering_pair = [(s1, s5), (s4, s2)]
    scattering_pair_rev = [(s5, s1), (s2, s4)]

    @test G_R(rand_num, scattering_pair, my_system) |> real ≈
          G_R(conj.(rand_num), scattering_pair_rev, my_system) |> real

    @test G_R(rand_num, scattering_pair, my_system) |> imag ≈
          -G_R(conj.(rand_num), scattering_pair_rev, my_system) |> imag

    @test δΓ(rand_num, my_system) |> Diagonal ≈ conj.(δΓ(conj.(rand_num), my_system) |> Diagonal)

end

@testset "No defects present" begin
    rand_num = ((rand() - 1 / 2) + 1im * (rand() - 1 / 2)) * 5
    my_system = mkGrapheneSystem(0.0, 0.0, Vector{Defect}([]))

    @test my_system.Δ == Array{GrapheneCoord}(undef, 0, 0)
    scattering_pair = [(s1, s5), (s4, s2)]
    @test δG_R(rand_num, scattering_pair, my_system) ==
          [(0.0 + 0.0im), (0.0 + 0.0im)]
    @test δρ_R_graphene(s2, my_system) == 0.0
    @test δρ_R_graphene(s3, my_system) == 0.0
end
