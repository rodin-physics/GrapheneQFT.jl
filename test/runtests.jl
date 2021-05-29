using Test, GrapheneQFT # This load both the test suite and our MyAwesomePackage

out = plusTwo(3)

@test out == 5
