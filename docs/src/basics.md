# Basic Usage
```@meta
CurrentModule = GrapheneQFT
DocTestSetup  = quote
    using GrapheneQFT, Plots
end
```
```@setup guide
using GrapheneQFT
```
```@setup guide2
using GrapheneQFT
using Plots
```
To view all functions and macros provided, run
```@repl guide
print(names(GrapheneQFT))
```
We can see that there are many functions provided and it can be easy to not know where to begin. An easy way to access basic documentation for each function is to enter help mode by typing `?` into the REPL and typing in the function name.
```
help?> GrapheneCoord
search: GrapheneCoord

 GrapheneCoord(u::Int, v::Int, sublattice::Sublattice)

 Lattice coordinate of a carbon atom.

 Each coordinate contains the sublattice index A or B, as well as the integer coefficients of the two basis vectors d\times(\pm 1 \hat{x} + \sqrt{3}\hat{y}) / 2 (u for +, v for
 -), with d = 2.46Å as the lattice constant.
```

## Coordinates and States in Graphene
We begin by learning how to define coordinates and states in graphene using the `GrapheneCoord` and `GrapheneState` structs respectively. Entering the help mode and searching for the documentation for `GrapheneCoord`, as we did above, shows that it is the lattice coordinate for a carbon atom in graphene with three fields: two integers as coefficients of the basis vectors and a sublattice index `A` or `B`. We can define coordinates and access their fields using `struct_name.field_name`:
```@repl coords
using GrapheneQFT#hide
c1 = GrapheneCoord(0, 0, A)
c2 = GrapheneCoord(2, -5, B)
c1.u
c2.v
c2.sublattice
```
To define a `GrapheneState`, we once again search for documentation by entering help mode to see the fields of this struct.
```
help?> GrapheneState
search: GrapheneState GrapheneSystem mkGrapheneSystem graphene_neighbors graphene_multiple_neighbors

  GrapheneState(coord::GrapheneCoord, spin::Spin)

  Quantum state of an electron in graphene, denoted by |u, v, L\rangle\otimes |\sigma\rangle in the drivation.

  The state is given by the GrapheneCoord of the orbital, as well as the electronic spin, which can take values SpinUp and SpinDown.
```
Following this, we see that the allowed values of the `spin` field are `SpinUp` and `SpinDown` and use them accordingly.
```@repl coords
s1 = GrapheneState(c1, SpinUp)
s2 = GrapheneState(c2, SpinDown)
```

## Defect Types
There are three kinds of possible defects that can be introduced within the GrapheneQFT framework: `ImpurityState`, `LocalSpin` and `Hopping`. The abstract `Defect` type encompasses all of these types. As usual, entering help mode and searching for documentation yields details on each defect type. The examples below give an overview of the available `Defect` types.
```@example guide
# Define coordinates of defect locations
c1 = GrapheneCoord(0, 0, A)
c2 = GrapheneCoord(0, 0, B)

# Define an ImpurityState with energy 0.2 eV coupled to c1 with energy 0.1 eV
imp1 = ImpurityState(0.2, [(0.1, c1)])

# Define a LocalSpin with a z spin component at c2 with coupling strength 0.1 eV
spin1 = LocalSpin(0.0, 0.0, 0.1, c2)

# Define a Hopping between c1 and c2 with modification of 0.2 eV in hopping energy
hop1 = Hopping(c1, c2, 0.2)

# Define an onsite energy of 0.05 eV on c1
ener1 = Hopping(c1, c1, 0.05)
nothing #hide
```

## Constructing a Graphene System
To begin calculating relevant quantities, we must first assemble the necessary elements and construct a `GrapheneSystem`. To do so, we make use of the `mkGrapheneSystem`, which, as the name suggests, constructs a `GrapheneSystem` by taking in the chemical potential and temperature of the system and a vector of `Defect`s. The easiest system to construct is, of course, one that is free of defects and is comprised of just pristine graphene.
```@example guide
# Define chemical potential and temperature
μ = 0.0
T = 0.0

# Define system
pristine_sys = mkGrapheneSystem(μ, T, Defect[])
```
Here, we have constructed a pristine system with chemical potential 0.1 eV and temperature at 0.0 K. Note that we have asserted the type of empty vector to be a vector of `Defect`s. As with all other structs, the fields of `GrapheneSystem` can be accessed in the usual way:
```@repl guide
pristine_sys = mkGrapheneSystem(0.1, 0.0, Defect[])#hide
pristine_sys.μ
pristine_sys.T
pristine_sys.Δ
pristine_sys.V
pristine_sys.scattering_states
pristine_sys.imps
```
Since a pristine system is not a very exciting one, we see that most fields are empty. However, the ability to access these properties will be useful in more complicated systems. When constructing a `GrapheneSystem`, the parametrization of defects is dependent on the nature of the problem and the types of defects present in the system we wish to investigate.
