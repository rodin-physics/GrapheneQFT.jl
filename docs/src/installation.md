# Installation
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
This guide is intended to provide the user with basic examples necessary to get started with the package, as well as intermediate examples that extend on the provided set of functions. For users interested in the derivation of the field theoretic calculations, the [formalism](@ref Formalism) is covered here. For users purely interested in trying out the package, the list of available functions and their relevance to physical quantities is available [here](@ref api).


GrapheneQFT can be installed with the Julia package manager. In the Julia REPL, type `]` and run
```
pkg> add GrapheneQFT
```
While still in the package manager, the test suite can be run with `test GrapheneQFT` to ensure the package is working correctly. After installation, the package is imported in the usual way:
```@repl
using GrapheneQFT
```
We are now ready to explore GrapheneQFT. For the purposes of this guide, it is assumed that the user has a plotting toolset such as [Plots.jl](https://github.com/JuliaPlots/Plots.jl) installed.
