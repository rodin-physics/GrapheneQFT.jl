# Examples
```@meta
CurrentModule = GrapheneQFT
DocTestSetup  = quote
    using GrapheneQFT
    using Plots
end
```
```@setup guide
using GrapheneQFT
```
```@setup guide2
using GrapheneQFT
using Plots
```
The following set of examples focus on extracting experimentally relevant quantities using the functionality provided by GrapheneQFT and visualizing the output. Here, we use Plots.jl although any plotting backend can be used. Note that the syntax for plotting differs across backends.
```@contents
Pages = ["examples.md"]
Depth = 3
```

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
We plot this and recover the familiar graphene spectral function.
```@example spec_func
using Plots
plot(ωs, spec_func, color = "blue", label = "$state")
```


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
```

### Spatial Map
In addition to a line plot of the spectral function across an energy range for a chosen `GrapheneState`, it is useful to plot spatial maps for a chosen energy. This particular example shows how to build functions on top of the existing functionality. Before the computation, we first define computation settings. These include a grid of `GrapheneCoord`s for both `A` and `B` sublattices and a function to calculate the graphene spectral function for a `GrapheneState`.
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
ω = 0.2

# Define GrapheneSystem with an ImpurityState at c1 and its neighbors
test_sys = mkGrapheneSystem(0.1, 0.0, Defect[ImpurityState(0.2, vcat((0.15, c1), [(0.1, x) for x in c1_neighbors]))])

# Calculation
resA = map(x -> spectral_graphene_spinUp(ω, x, test_sys), coord_A)
resB = map(x -> spectral_graphene_spinUp(ω, x, test_sys), coord_B)
signal = vcat(resA, resB)
```
To plot, we make use of the `crystal_to_cartesian` function to correctly place the carbon atoms.
```@example spatial_spec
c1 = GrapheneCoord(0,0,A) #hide
c1_neighbors = graphene_neighbors(c1) #hide
ω = 0.2 #hide
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

# Plotting
scatter(X_latt, Y_latt,
        marker_z = signal, color = :viridis,
        markersize = 4, markerstrokewidth = 0,
        xrange = (-30, 30), yrange = (-30, 30),
        label = nothing, clims = (0.009, 0.0098),
        xlabel = "Distance(Å)", ylabel = "Distance(Å)",
        guidefontsize = 8, aspectratio = 1)
```

## Free Energy
We use the `δF` function to calculate the variation in system free energy due to the defects.
### Line Plot
The script below plots how ``\delta F`` varies as the impurity-lattice coupling strength is tuned in a system with a single `ImpurityState`.
```@example guide2
# Define impurity energy and coordinate
ϵ = 0.2
c1 = GrapheneCoord(0,0,A)

# Define function to calculate δF for each coupling strength
function δF_imp_coupling(coup::Float64)
    imp = ImpurityState(ϵ, [(coup, c1)])
    test_sys = mkGrapheneSystem(0.0, 0.0, Defect[imp])
    return δF(test_sys)
end

# Define range of coupling strengths and map onto values
couplings = range(0.0, 0.5, step = 0.01)
res = map(δF_imp_coupling, couplings)

plot(couplings, res,
    xlabel = "Coupling Strength (eV)",
    ylabel = "δF (eV)", title = "Single ImpurityState",
    guidefontsize = 8, label = nothing)
```

### Heatmap
Similar to above, we can vary both the system's chemical potential and coupling strength to produce a heatmap of ``\delta F`` values.
```@example guide2
ϵ = 0.2 #hide
c1 = GrapheneCoord(0,0,A) #hide
# Define function to calculate δF for each coupling strength
function δF_imp_mu_coupling(mu::Float64, coup::Float64)
    imp = ImpurityState(ϵ, [(coup, c1)])
    test_sys = mkGrapheneSystem(mu, 0.0, Defect[imp])
    return δF(test_sys)
end

# Define range of coupling strengths and map onto values
couplings = range(0.0, 0.5, step = 0.01)
μ_values = range(0.0, 0.5, step = 0.005)

# Plotting
heatmap(μ_values, couplings, δF_imp_mu_coupling,
        xlabel = "μ (eV)", ylabel = "Coupling Strength (eV)",
        colorbar_title = "\nδF (eV)", guidefontsize = 8,
        right_margin = 2.5Plots.mm)
```

## Charge Density Modulation
The defect-induced variation to the charge density is calculated using the `δρ_R_graphene` function.
### Line Plot
The script below plots the variation in total ``\delta \rho_\mathbf{R}`` as a function of distance from the impurity in a single-impurity system. The total charge density is given by the sum of spin-resolved charge density values.
```@example guide2
# Define system parameters
ϵ = 0.2
c1 = GrapheneCoord(0,0,A)
imp = ImpurityState(0.2, [(0.1, c1)])
dist = -10:1:10

# Define GrapheneSystem
test_sys = mkGrapheneSystem(0.2, 0.0, Defect[imp])

# Calculation for spin-resolved charge density
res_SpinUp =  map(x -> δρ_R_graphene(GrapheneState(GrapheneCoord(x,0,A), SpinUp), test_sys), dist)
res_SpinDown =  map(x -> δρ_R_graphene(GrapheneState(GrapheneCoord(x,0,A), SpinDown), test_sys), dist)

# Plotting total charge density variation
plot(dist, (res_SpinUp .+ res_SpinDown);
    marker = (:circle, 5), title = "Total Charge Density Variation",
    xlabel = "Distance (d)", ylabel = "δρ (electrons/atom)",
    guidefontsize = 8, label = nothing)
```

## Miscellaneous
### Animations
It is often useful to create animations to visualize how a quantity changes with respect to some variable. One of the easiest ways to begin creating animations is to use the `@gif` and `@animate` macros in Plots.jl . The script below creates a GIF to show how the spectral function varies as the energy level of an `ImpurityState` is tuned.
```@example guide2
# Define defect as function that takes in the defect energy
imp_ener(ϵ) = ImpurityState(ϵ, [(0.8, GrapheneCoord(0,0,A))])

# Define range of energies for defect
ϵs = range(-5.0, 5.0, step = 0.05)

# Define range of energies to calculate spectral function for
ωs = range(-9.0, 9.0, step = 0.05)

# Define spectral function
function spectral_graphene_spinUp(energy::Float64, sys::GrapheneSystem)
    state = GrapheneState(GrapheneCoord(0,0,A), SpinUp)
    res = (-2/π) * G_R(energy + 1im * 1e-6, [(state, state)], sys)
    return imag(res[1])
end

# Define a plotting function for convenience
function plotter_func(ϵ::Float64)
    # Define GrapheneSystem
    test_sys = mkGrapheneSystem(0.1, 0.0, Defect[imp_ener(ϵ)])

    # Calculate and plot spectral function
    res = map(x -> spectral_graphene_spinUp(x, test_sys), ωs)
    plot(ωs, res,
        color = :blue,
        label = nothing,
        title = "ϵ = $ϵ eV",
        ylims = (-0.01, 0.9))
end

# Animation
anim = @animate for ii in ϵs
                    plotter_func(ii)
                end

gif(anim, fps = 20)
```

### Parallel Computing
The functionality of GrapheneQFT.jl lends itself well to parallel computing, especially with the use of parallel mapping. The module [`Distributed`](https://docs.julialang.org/en/v1/manual/distributed-computing/) is an implementation of distributed parallel computing that is part of Julia's standard library. The user manual shows how to use the module and benefit from parallel computing methods such as `pmap`.
