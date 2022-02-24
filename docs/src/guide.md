# Getting Started
```@meta
CurrentModule = GrapheneQFT
DocTestSetup  = quote
    using GrapheneQFT
end
```
```@setup guide
using GrapheneQFT
```
```@setup guide2
using GrapheneQFT
using Plots
```
[GrapheneQFT.jl](https://github.com/rodin-physics/GrapheneQFT.jl) is an extendable module that facilitates the field theoretic treatment of impurities on monolayer graphene. This guide is intended to provide the user with basic examples necessary to get started with the package, as well as intermediate examples that extend on the provided set of functions. For users interested in the derivation of the field theoretic calculations, the [formalism](@ref Formalism) is covered here. For users purely interested in trying out the package, the list of available functions and their relevance to physical quantities is available [here](@ref api).

GrapheneQFT can be installed with the Julia package manager. In the Julia REPL, type `]` and run
```
pkg> add GrapheneQFT
```
While still in the package manager, the test suite can be run with `test GrapheneQFT` to ensure the package is working correctly. After installation, the package is imported in the usual way:
```@repl
using GrapheneQFT
```
We are now ready to explore GrapheneQFT. For the purposes of this guide, it is assumed that the user has a plotting toolset such as [Plots.jl](https://github.com/JuliaPlots/Plots.jl) installed.

# Basic Usage
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

# Examples
The following set of examples focus on extracting experimentally relevant quantities using the functionality provided by GrapheneQFT and visualizing the output.
## Spectral Function
### Single Line Plot
To calculate the spectral function for graphene, we recall that it is given by ``-2\textrm{Im}[G (\omega + i0^{+}, \mathbf{R})]``, where we use the `G_R` function to calculate the graphene Green's function. We first define the `GrapheneState` and `GrapheneSystem` we are interested in calculating the spectral function for, as well as a range of energies.
```@example spec_func
using GrapheneQFT#hide
# Define a state to calculate spectral function for
state = GrapheneState(GrapheneCoord(0,1,A), SpinUp)

# Define a GrapheneSystem
pristine_sys = mkGrapheneSystem(0.0, 0.0, Defect[])

# Define range of energies to calculate spectral function for
ωs = range(-9.5, 9.5, step = 0.05)
nothing#hide
```
Next, we can go ahead and calculate the spectral function by mapping the range of energies onto the `G_R` function.
```@example spec_func
# Calculate Green's function for chosen GrapheneState and GrapheneSystem
greens_func = map(ω -> G_R(ω + 1im * 1e-6, [(state, state)], pristine_sys)[1], ωs)

# Calculate spectral function
spec_func = (-2/π) .* imag.(greens_func)
nothing#hide
```
We can plot this using any backend of choice. Here, we use Plots.jl and recover the familiar graphene spectral function.
```@example spec_func
using Plots
plot(ωs, spec_func, color = "blue", label = "$state")
savefig("spec_func.png"); nothing #hide
```

![](spec_func.png)

### Multiple Line Plots
The astute reader would have noticed that `G_R` takes in a vector of `GrapheneState` Tuples, and the use of `[1]` in the script was to extract the first (and only) value in the result. The ability to generalize makes it easy to calculate `G_R` for multiple `GrapheneStates` simultaneously, with a few manipulations of the data. A more interesting system is used as an example below.
```@example guide
using Plots#hide
# Define coordinates and states
c1 = GrapheneCoord(0,0,A)
c2 = GrapheneCoord(1,2,B)

s1 = GrapheneState(c1, SpinUp)
s2 = GrapheneState(c1, SpinDown)
s3 = GrapheneState(c2, SpinUp)

# Define Defects
spin1 = LocalSpin(0.0, 0.0, 0.5, c1)
imp1 = ImpurityState(0.2,[(0.1, c2)])
hoppings = [Hopping(c1, x, 0.15) for x in graphene_neighbors(c1)]

# Define system with multiple Defects
test_sys = mkGrapheneSystem(0.0, 0.0, vcat(spin1, imp1, hoppings))

# Define range of energies to calculate spectral function for
ωs = range(-9.5, 9.5, step = 0.05)

# Calculate Green's function
greens_func_all = map(ω -> G_R(ω + 1im * 1e-6, [(s1, s1), (s2, s2), (s3, s3)], test_sys), ωs)

# Calculate and separate individual spectral functions
spec_func_all = (-2/π) .* imag.(greens_func_all)
spec_func1 = getindex.(spec_func_all,1)
spec_func2 = getindex.(spec_func_all,2)
spec_func3 = getindex.(spec_func_all,3)

# Plotting
plot(ωs, spec_func1, color = "blue", label = "$s1")
plot!(ωs, spec_func2, color = "red", label = "$s2")
plot!(ωs, spec_func3, color = "green", label = "$s3")
savefig("spec_func_multiple.png"); nothing #hide

```
![](spec_func_multiple.png)

### Spatial Map
In addition to a line plot of the spectral function across an energy range for a chosen `GrapheneState`, it is useful to plot spatial maps for a chosen energy. This particular example shows how to build functions on top of the existing functionality. Before the computation, we first define computation settings. These include a grid of `GrapheneCoord`s for both `A` and `B` sublattices and function to calculate the graphene spectral function for a `GrapheneState`.
```@example spatial_spec
using GrapheneQFT, Plots #hide
# Define system chemical potential and temperature
μ = 0.0
T = 0.0

# Define grid points and coordinate grids
nPts = 20
range_nPts = -nPts:1:nPts

coord_A = [GrapheneCoord(u,v,A) for u in range_nPts, v in range_nPts] |> vec
coord_B = [GrapheneCoord(u,v,B) for u in range_nPts, v in range_nPts] |> vec

# Define function
function spectral_graphene_spinUp(energy::Float64, coord::GrapheneCoord, sys::GrapheneSystem)
    state = GrapheneState(coord, SpinUp)
    res = (-2/π) * G_R(energy + 1im * 1e-6, [(state, state)], sys)
    return imag(res[1])
end
nothing#hide
```
We can now define the system and calculate the spectral function over the coordinate grids.
```julia
# Define coordinates and states
c1 = GrapheneCoord(0,0,A)
c1_neighbors = graphene_neighbors(c1)

# Define energy (eV) for spatial map
ω = 0.1

# Define GrapheneSystem with an ImpurityState at c1 and its neighbors
test_sys = mkGrapheneSystem(0.1, 0.0, Defect[ImpurityState(0.2, vcat((0.15, c1), [(0.1, x) for x in c1_neighbors]))])

# Calculation
resA = map(x -> spectral_graphene_spinUp(ω, x, test_sys), coord_A)
resB = map(x -> spectral_graphene_spinUp(ω, x, test_sys), coord_B)
signal = vcat(resA, resB)
```
To plot, we make use of the `crystal_to_cartesian` function.
```@example spatial_spec
c1 = GrapheneCoord(0,0,A) #hide
c1_neighbors = graphene_neighbors(c1) #hide
ω = 0.1 #hide
test_sys = mkGrapheneSystem(0.1, 0.0, Defect[ImpurityState(0.2, vcat((0.15, c1), [(0.1, x) for x in c1_neighbors]))]) #hide
resA = map(x -> spectral_graphene_spinUp(ω, x, test_sys), coord_A) #hide
resB = map(x -> spectral_graphene_spinUp(ω, x, test_sys), coord_B) #hide
signal = vcat(resA, resB) #hide
# Define Cartesian coordinates
A_latt = crystal_to_cartesian.(coord_A)
B_latt = crystal_to_cartesian.(coord_B)

full_lattice = vcat(A_latt, B_latt)
X_latt = first.(full_lattice)
Y_latt = last.(full_lattice)

scatter(X_latt, Y_latt,
        marker_z = signal, markersize = 3,
        xrange = (-40, 40), yrange = (-40, 40),
        label = nothing, clims = (0.00468, 0.0047),
        aspectratio = 1)
savefig("spec_func_spatial.png"); #hide
```
![](spec_func_spatial.png)
