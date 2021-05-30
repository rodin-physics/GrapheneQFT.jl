var documenterSearchIndex = {"docs":
[{"location":"anotherPage.html#The-GrapheneQFT-Module","page":"Another page","title":"The GrapheneQFT Module","text":"","category":"section"},{"location":"anotherPage.html","page":"Another page","title":"Another page","text":"  GrapheneQFT","category":"page"},{"location":"anotherPage.html#GrapheneQFT","page":"Another page","title":"GrapheneQFT","text":"GrapheneQFT\n\nt = 28eV\n\n\n\n\n\n","category":"module"},{"location":"anotherPage.html#Module-Index","page":"Another page","title":"Module Index","text":"","category":"section"},{"location":"anotherPage.html","page":"Another page","title":"Another page","text":"Modules = [GrapheneQFT]\nOrder   = [:constant, :type, :function, :macro]","category":"page"},{"location":"anotherPage.html#Detailed-API","page":"Another page","title":"Detailed API","text":"","category":"section"},{"location":"anotherPage.html","page":"Another page","title":"Another page","text":"Modules = [GrapheneQFT]\nOrder   = [:constant, :type, :function, :macro]","category":"page"},{"location":"anotherPage.html#GrapheneQFT.GrapheneCoord","page":"Another page","title":"GrapheneQFT.GrapheneCoord","text":"Lattice coordinate of a carbon atom, generated using graphene_A or graphene_B. Each coordinate contains the sublattice index, as well as the integer coefficients of the two basis vectors dtimes(pm 1 hatx + sqrt3haty)  2, where d = 246Å is the graphene lattice constant.\n\n\n\n\n\n","category":"type"},{"location":"anotherPage.html#GrapheneQFT.Location","page":"Another page","title":"GrapheneQFT.Location","text":"Location(x::Float64, y::Float64, z::Float64)\n\nA structure describing a point in 3D space, with the lengths in Å\n\n\n\n\n\n","category":"type"},{"location":"anotherPage.html#GrapheneQFT.crystal_to_cartesian-Tuple{GrapheneQFT.GrapheneCoord}","page":"Another page","title":"GrapheneQFT.crystal_to_cartesian","text":"crystal_to_cartesian(coord::GrapheneCoord)\n\nConvert a crystal coordinate to a cartesian one.\n\nArguments\n\ncoord: a GrapheneCoord that needs to be converted to a 3D Location\n\nOutput\n\nLocation of the carbon atom\n\n\n\n\n\n","category":"method"},{"location":"anotherPage.html#GrapheneQFT.graphene_A-Tuple{Int64, Int64}","page":"Another page","title":"GrapheneQFT.graphene_A","text":"graphene_A(u::Int, v::Int)\n\nCreate a GrapheneCoord for an atom belonging to sublattice A at the unit cell (u, v)\n\nArguments\n\nu: coefficient of basis vector dtimes(1 hatx + sqrt3haty)  2\nv: coefficient of basis vector dtimes(-1 hatx + sqrt3haty)  2\n\nOutput\n\nGrapheneCoord of the carbon atom\n\n\n\n\n\n","category":"method"},{"location":"anotherPage.html#GrapheneQFT.graphene_B-Tuple{Int64, Int64}","page":"Another page","title":"GrapheneQFT.graphene_B","text":"graphene_B(u::Int, v::Int)\n\nCreate a GrapheneCoord for an atom belonging to sublattice B at the unit cell (u, v)\n\nArguments\n\nu: coefficient of basis vector dtimes(1 hatx + sqrt3haty)  2\nv: coefficient of basis vector dtimes(-1 hatx + sqrt3haty)  2\n\nOutput\n\nGrapheneCoord of the carbon atom\n\n\n\n\n\n","category":"method"},{"location":"anotherPage.html#GrapheneQFT.neighbors-Tuple{GrapheneQFT.GrapheneCoord}","page":"Another page","title":"GrapheneQFT.neighbors","text":"neighbors(atom::GrapheneCoord)\n\nDetermine the nearest neighbors of an atom.\n\nArguments\n\natom: GrapheneCoord whose nearest neighbors are to be determined\n\nOutput\n\nA Vector of GrapheneCoord containing the neighbors of atom\n\n\n\n\n\n","category":"method"},{"location":"anotherPage.html#GrapheneQFT.propagator-Tuple{GrapheneQFT.GrapheneCoord, GrapheneQFT.GrapheneCoord, Any}","page":"Another page","title":"GrapheneQFT.propagator","text":"propagator(a_l::GrapheneCoord, a_m::GrapheneCoord, z)\n\nThe propagator function picks out the correct element of the Ξ matrix based on the sublattices of the graphene coordinates\n\n\n\n\n\n","category":"method"},{"location":"index.html#GrapheneQFT.jl","page":"Index","title":"GrapheneQFT.jl","text":"","category":"section"},{"location":"index.html","page":"Index","title":"Index","text":"Documentation for GrapheneQFT.jl","category":"page"}]
}
