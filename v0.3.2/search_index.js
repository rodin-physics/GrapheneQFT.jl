var documenterSearchIndex = {"docs":
[{"location":"formalism.html#Formalism","page":"Formalism","title":"Formalism","text":"","category":"section"},{"location":"formalism.html#General-System","page":"Formalism","title":"General System","text":"","category":"section"},{"location":"formalism.html","page":"Formalism","title":"Formalism","text":"Instead of starting directly with the problem of defects in graphene, we start with a more general scenario which will make the derivation more transparent. Consider a pristine system described by a second-quantized Hamiltonian sum_jk c^dagger_j H_jkc_k, where j and k label the states of the system, c_j and c_j^dagger are the corresponding annihilation and creation operators, and H_jk gives the coupling between the states. Next, we introduce two types of perturbation to this system. First, we can modify the coupling between the states. Second, we introduce a second system described by a Hamiltonian sum_jk g^dagger_j h_jkg_k and couple it to the original one. The total Hamiltonian, then, can be written as","category":"page"},{"location":"formalism.html","page":"Formalism","title":"Formalism","text":"hatmathcalH = sum_jk c^dagger_j H_jkc_k\n    +\n    sum_jk c^dagger_j Delta_jkc_k\n    +\n    sum_jk left(c^dagger_j V_jkg_k + g^dagger_k V_jk^* c_jright)\n    +\n    sum_jk g^dagger_j h_jkg_k","category":"page"},{"location":"formalism.html","page":"Formalism","title":"Formalism","text":"where the second term in the first line describes the modified coupling in the original system and the first term in the second line gives the coupling between the two systems. The Hamiltonian can be written more compactly if we assemble the operators and the coupling elements into vectors and matrices, respectively:","category":"page"},{"location":"formalism.html","page":"Formalism","title":"Formalism","text":"hatmathcalH = mathbfc^dagger Hmathbfc + mathbfc^dagger Delta mathbfc\n    +left(mathbfc^dagger Vmathbfg + mathbfg^dagger V^dagger mathbfcright) + mathbfg^dagger hmathbfg","category":"page"},{"location":"formalism.html","page":"Formalism","title":"Formalism","text":"Next, we transcribe the Hamiltonian into the imaginary-time action","category":"page"},{"location":"formalism.html","page":"Formalism","title":"Formalism","text":"S = sum_nbarboldsymbolpsi_n overbraceleft(-iomega_n - mu + Hright)^-G^-1_iomega_n + muboldsymbolpsi_n + barboldsymbolpsi_nDelta boldsymbolpsi_n\n    +left(barboldsymbolpsi_n Vboldsymbolphi_n + barboldsymbolphi_n V^dagger boldsymbolpsi_nright) + barboldsymbolphi_n overbraceleft(-iomega_n - mu + hright)^-Gamma_iomega_n+mu^-1boldsymbolphi_n\n\t\n\t=\n\tsum_n\n\tbeginpmatrix\n\tbarboldsymbolpsi_n  barboldsymbolphi_n\n\tendpmatrix\n\tbeginpmatrix\n        -G_iomega_n+mu^-1 + Delta  V\n        \n        V^dagger  -Gamma_iomega_n+mu^-1\n    endpmatrix\n\tbeginpmatrix\n\tboldsymbolpsi_n  boldsymbolphi_n\n\tendpmatrix\n\t","category":"page"},{"location":"formalism.html","page":"Formalism","title":"Formalism","text":"where omega_n are the fermionic Matsubara frequencies, mu is the chemical potential, and boldsymbolphi_n and boldsymbolpsi_n (barboldsymbolphi_n and barboldsymbolpsi_n) are vectors of Grassmann numbers corresponding to mathbfg and mathbfc (mathbfg^dagger and mathbfc^dagger). We also identify G_z and Gamma_z as the Green's functions for the two systems.","category":"page"},{"location":"formalism.html","page":"Formalism","title":"Formalism","text":"Exponentiating -S and integrating over all the Grassmann variables yields the partition function","category":"page"},{"location":"formalism.html","page":"Formalism","title":"Formalism","text":"mathcalZ\n    =prod_n\n    leftbeta\n    beginpmatrix\n        -G_iomega_n+mu^-1 + Delta  V\n        \n        V^dagger  -Gamma_iomega_n+mu^-1\n    endpmatrix\n    right","category":"page"},{"location":"formalism.html","page":"Formalism","title":"Formalism","text":"where beta = 1  (k_BT). We can identify the full system Green's function from the determinant in mathcalZ","category":"page"},{"location":"formalism.html","page":"Formalism","title":"Formalism","text":"mathbfG_z\n    =\n    beginpmatrix\n        G_z^-1 - Delta  -V\n        \n        -V^dagger  Gamma_z^-1\n    endpmatrix^-1\n    =\n    beginpmatrix\n         mathcalG_z  mathcalG_z V Gamma_z\n        \n         Gamma_zV^dagger mathcalG_z  Gamma_z+ Gamma_zV^dagger mathcalG_z V Gamma_z\n    endpmatrix","category":"page"},{"location":"formalism.html","page":"Formalism","title":"Formalism","text":"where","category":"page"},{"location":"formalism.html","page":"Formalism","title":"Formalism","text":"mathcalG_z\n    =left(G_z^-1 - Delta - VGamma_z V^daggerright)^-1\n    = G_z+G_zleft(Delta + VGamma_z V^daggerright)\n\tleft1-G_zleft(Delta + VGamma_z V^daggerright)right^-1G_z ","category":"page"},{"location":"formalism.html","page":"Formalism","title":"Formalism","text":"is the full Green's function of the bulk system including the effects of defects and perturbations. Similarly, the bottom right block in mathbfG_z corresponds to the full Green's function of the impurity states including their coupling to the perturbed bulk system:","category":"page"},{"location":"formalism.html","page":"Formalism","title":"Formalism","text":"Gamma_z + Gamma_z V^dagger mathcalG_z V Gamma_z\n= Gamma_z + Gamma_z V^dagger G_z left1 - (Delta + V Gamma_z V^dagger) G_z right^-1 V Gamma_z  ","category":"page"},{"location":"formalism.html","page":"Formalism","title":"Formalism","text":"From the partition function, we can also write down the Helmholtz free energy F = -beta^-1 ln mathcalZ:","category":"page"},{"location":"formalism.html","page":"Formalism","title":"Formalism","text":"F = -beta^-1sum_n lnleftbeta\n    beginpmatrix\n        -G_iomega_n+mu^-1 + Delta  V\n        \n        V^dagger  -Gamma_iomega_n+mu^-1\n    endpmatrix\n    right\n    \n    = -beta^-1sum_n lnBiggbeta\n    beginpmatrix\n        -G_iomega_n+mu^-1  0\n        \n        0  -Gamma_iomega_n+mu^-1\n    endpmatrix\n    left\n    1\n    +\n    beginpmatrix\n        -G_iomega_n+mu  0\n        \n        0  -Gamma_iomega_n+mu\n    endpmatrix\n    beginpmatrix\n        Delta  V\n        \n        V^dagger  0\n    endpmatrix\n    right\n    Bigg","category":"page"},{"location":"formalism.html","page":"Formalism","title":"Formalism","text":"Removing the part of F corresponding to the free energy of the two isolated system in the absence of any perturbation yields the defect- and coupling-induced modification to F","category":"page"},{"location":"formalism.html","page":"Formalism","title":"Formalism","text":"delta F\n    = -beta^-1sum_n lnleft\n    beginpmatrix\n       1- G_iomega_n+muDelta  -G_iomega_n+muV\n        \n       - Gamma_iomega_n+muV^dagger  1\n    endpmatrix\n    right\n    \n    = -beta^-1sum_n lnleft\n      1- G_iomega_n+muleft(Delta+VGamma_iomega_n+muV^daggerright)\n    right","category":"page"},{"location":"formalism.html#Graphene","page":"Formalism","title":"Graphene","text":"","category":"section"},{"location":"formalism.html","page":"Formalism","title":"Formalism","text":"Having performed the required manipulations for a general system, we now focus on graphene. Because the electronic properties of graphene are dominated by the carbon pi orbitals, the electronic states in a pristine system can be described by mathbfrLrangleotimessigmarangle, where mathbfr is the coordinate of the unit cell hosting the orbital, L is the sublattice of the atom, and sigma is the spin of the electron. This is the basis corresponding to the operators c_k in the derivation above.","category":"page"},{"location":"formalism.html","page":"Formalism","title":"Formalism","text":"The next required ingredient is the graphene Green's function G_z. Since G_z = (z - H)^-1, it is a matrix whose elements give the propagation amplitudes between graphene states. The matrix elements can be calculated using","category":"page"},{"location":"formalism.html","page":"Formalism","title":"Formalism","text":"langlemathbfr Lotimeslanglesigma(z - hatH)^-1mathbfrLrangleotimessigmarangle = frac1Nsum_mathbfqqlanglemathbfq L otimes langle sigmae^imathbfqcdotmathbfr(z - hatH)^-1e^-imathbfqcdotmathbfrmathbfqLrangleotimessigmarangle","category":"page"},{"location":"formalism.html","page":"Formalism","title":"Formalism","text":"where we Fourier-transformed the real-space states so that the sums over the momenta mathbfq and mathbfq run over the entire Brillouin zone. Because the momentum-space Hamiltonian is diagonal in mathbfq and sigma, we can write","category":"page"},{"location":"formalism.html","page":"Formalism","title":"Formalism","text":"langlemathbfr Lotimeslanglesigma(z - hatH)^-1mathbfrLrangleotimessigmarangle = frac1Nsum_mathbfqlanglemathbfq Le^imathbfqcdotmathbfr(z - hatH)^-1e^-imathbfqcdotmathbfrmathbfqLrangle delta_sigmasigma\n\n= langle Lfrac1Nsum_mathbfqe^imathbfqcdotleft(mathbfr-mathbfrright)(z - H_mathbfq)^-1Lrangle delta_sigmasigma","category":"page"},{"location":"formalism.html","page":"Formalism","title":"Formalism","text":"Using","category":"page"},{"location":"formalism.html","page":"Formalism","title":"Formalism","text":"H_mathbfq =beginpmatrix 0  -t f_mathbfq  -t f^*_mathbfq  0endpmatrix ","category":"page"},{"location":"formalism.html","page":"Formalism","title":"Formalism","text":"where t = 28eV is the nearest-neighbor hopping integral, f_mathbfq = 1 + e^imathbfqcdotmathbfd_1+e^imathbfqcdotmathbfd_2, and mathbfd_12=dleft(pm1 sqrt3right)2 are the lattice vectors, we have","category":"page"},{"location":"formalism.html","page":"Formalism","title":"Formalism","text":"frac1Nsum_mathbfq left(z -  H_mathbfq right)^-1 e^ileft(mathbfr_k - mathbfr_jright)cdotmathbfq  = frac1Nsum_mathbfq\nbeginpmatrix\nz - tf_mathbfq\n\n- tf_mathbfq^*  z\nendpmatrixfrac1z^2 - t^2 leftf_mathbfqright^2 e^imathbfr_kjcdotmathbfq","category":"page"},{"location":"formalism.html","page":"Formalism","title":"Formalism","text":"To perform the summation over mathbfq, we first introduce","category":"page"},{"location":"formalism.html","page":"Formalism","title":"Formalism","text":"Omega^uv_z =\n\tfrac1Nsum_mathbfqinmathrmBZ\n\tfrac\n\t\te^imathbfqcdot left(umathbfd_1 + vmathbfd_2right)\n\t\n\tz^2 - t^2left f_mathbfqright^2","category":"page"},{"location":"formalism.html","page":"Formalism","title":"Formalism","text":"with umathbfd_1 + vmathbfd_2 = fracd2left(u - v sqrt3left(u+vright)right). Using mathbfqcdot left(umathbfd_1 + vmathbfd_2right)  = fracd2leftleft(u - vright)q_x + sqrt3left(u+vright)q_yright and turning the momentum sum into an integral yields","category":"page"},{"location":"formalism.html","page":"Formalism","title":"Formalism","text":"Omega^uvleft(zright)\n\t= frac1left(2piright)^2oint dx oint dy\n\tfrac\n\te^i leftleft(u - vright)x + left(u+vright)yright\n\tz^2 - t^2left(1 + 4cos^2 x + 4 cos xcos y right)","category":"page"},{"location":"formalism.html","page":"Formalism","title":"Formalism","text":"From","category":"page"},{"location":"formalism.html","page":"Formalism","title":"Formalism","text":"oint dtheta frace^ilthetaW-costheta = 2pi fracleft(W - sqrtW - 1sqrtW + 1right)^lsqrtW - 1sqrtW + 1","category":"page"},{"location":"formalism.html","page":"Formalism","title":"Formalism","text":"we get","category":"page"},{"location":"formalism.html","page":"Formalism","title":"Formalism","text":"Omega^uv_z = frac12pifrac14t^2\n\toint dx frace^ileft(u - vright)xcos xfracleft(W - sqrtW - 1sqrtW + 1right)^u+vsqrtW - 1sqrtW + 1\n\t\n\tW = fracfracz^2t^2-14cos x-cos x","category":"page"},{"location":"formalism.html","page":"Formalism","title":"Formalism","text":"Finally, for mathbfr = umathbfd_1 + vmathbfd_2, we have","category":"page"},{"location":"formalism.html","page":"Formalism","title":"Formalism","text":"frac1Nsum_mathbfq left(z -  H_mathbfq right)^-1 e^imathbfrcdotmathbfq  \n\t=\n\tbeginpmatrix\n\t\tzOmega^uv_z\n\t\t\n\t\t- tleftOmega^uv_z + Omega^uv_+z right\n\t\t\n\t\t- tleftOmega^uv_z + Omega^uv_-zright\n\t\t\n\t\tzOmega^uv_z\n\tendpmatrix\n\t\n\tOmega^uv_pmz\n\t=\n\t frac12pifrac14t^2\n\toint dx 2e^ileft(u - vright)xfracleft(W - sqrtW - 1sqrtW + 1right)^u+vpm 1sqrtW - 1sqrtW + 1","category":"page"},{"location":"formalism.html","page":"Formalism","title":"Formalism","text":"Using the multiples of the basis vectors to describe the electronic states (mathbfrLrangleotimessigmaranglerightarrow u v Lrangleotimessigmarangle), the matrix elements langle u v Lotimeslanglesigma(z - hatH)^-1u v Lrangleotimessigmarangle of the Green's function become","category":"page"},{"location":"formalism.html","page":"Formalism","title":"Formalism","text":"langle Lbeginpmatrix\n\t\tzOmega^u-uv-v_z\n\t\t\n\t\t- tleftOmega^u-uv-v_z + Omega^u-uv-v_+z right\n\t\t\n\t\t- tleftOmega^u-uv-v_z + Omega^u-uv-v_-zright\n\t\t\n\t\tzOmega^u-uv-v_z\n\tendpmatrixLrangle delta_sigmasigma","category":"page"},{"location":"formalism.html","page":"Formalism","title":"Formalism","text":"Depending on the sublattices of the two states, we pick out the appropriate element of the 2times 2 matrix above.","category":"page"},{"location":"formalism.html#Graphene-Green's-Function","page":"Formalism","title":"Graphene Green's Function","text":"","category":"section"},{"location":"formalism.html","page":"Formalism","title":"Formalism","text":"Following the general derivation above, the individual elements of the full real-space Green's function are","category":"page"},{"location":"formalism.html","page":"Formalism","title":"Formalism","text":"mathcalG_z^jk\n=  G_z^jk+sum_lmG_z^jlleftleft left(Delta +  VGamma_zV^daggerright)^-1-\nG_zright^-1right_lmG_z^mk","category":"page"},{"location":"formalism.html","page":"Formalism","title":"Formalism","text":"where j, k, l, and m are compact labels for the graphene states u v Lrangleotimessigmarangle. As written, the sum over l and m includes all the states in the system. However, the expression can be made considerably simpler. The term Delta + VGamma_z V^dagger is a square matrix including all the system states. Let us rearrange the elements in the matrix to give it a block diagonal form:","category":"page"},{"location":"formalism.html","page":"Formalism","title":"Formalism","text":"Delta + VGamma_zV^dagger rightarrow beginpmatrix000Aendpmatrix","category":"page"},{"location":"formalism.html","page":"Formalism","title":"Formalism","text":"where we created a square matrix A that contains only the perturbed states. Consequently,","category":"page"},{"location":"formalism.html","page":"Formalism","title":"Formalism","text":"left(Delta + VGamma_zV^daggerright)^-1 rightarrow beginpmatrixinfty 00A^-1endpmatrix","category":"page"},{"location":"formalism.html","page":"Formalism","title":"Formalism","text":"and","category":"page"},{"location":"formalism.html","page":"Formalism","title":"Formalism","text":"left left(Delta +  VGamma_zV^daggerright)^-1-\nG_zright^-1 rightarrow beginpmatrix0 00 leftA^-1 - tildeG_zright^-1endpmatrix","category":"page"},{"location":"formalism.html","page":"Formalism","title":"Formalism","text":"where tildeG_z is a portion of G_z that only includes the perturbed states. Hence, the lm summation only needs to run over these, too:","category":"page"},{"location":"formalism.html","page":"Formalism","title":"Formalism","text":"mathcalG_z^jk\n=  G_z^jk+sum_lminmathrmpertG_z^jlleftleft left(tildeDelta +  tildeVGamma_ztildeV^daggerright)^-1-\ntildeG_zright^-1right_lmG_z^mk","category":"page"},{"location":"formalism.html","page":"Formalism","title":"Formalism","text":"Here, the tilde over the matrices indicates that only the elements corresponding to the perturbed states are retained. Note that a particular state is included in tildeV and tildeDelta even if it is perturbed by only one of the terms!","category":"page"},{"location":"formalism.html#Impurity-Green's-Function","page":"Formalism","title":"Impurity Green's Function","text":"","category":"section"},{"location":"formalism.html","page":"Formalism","title":"Formalism","text":"For the full impurity Green's function, we follow the discussion above to obtain","category":"page"},{"location":"formalism.html","page":"Formalism","title":"Formalism","text":"Gamma_z + Gamma_z tildeV^dagger tildeG_z left1 - (tildeDelta + tildeV Gamma_z tildeV^dagger) tildeG_z right^-1 tildeV Gamma_z  ","category":"page"},{"location":"api.html#API","page":"API","title":"API","text":"","category":"section"},{"location":"api.html","page":"API","title":"API","text":"Modules = [GrapheneQFT]\nOrder   = [:constant, :type, :function, :macro]","category":"page"},{"location":"api.html#Pristine-Graphene","page":"API","title":"Pristine Graphene","text":"","category":"section"},{"location":"api.html","page":"API","title":"API","text":"  GrapheneCoord\n  GrapheneState\n  Defect\n  ImpurityState\n  LocalSpin\n  Hopping\n  GrapheneSystem\n  graphene_neighbors\n  graphene_multiple_neighbors\n  crystal_to_cartesian\n  mkGrapheneSystem\n  peierls_phase","category":"page"},{"location":"api.html#GrapheneQFT.GrapheneCoord","page":"API","title":"GrapheneQFT.GrapheneCoord","text":"GrapheneCoord(u::Int, v::Int, sublattice::Sublattice)\n\nLattice coordinate of a carbon atom.\n\nEach coordinate contains the sublattice index A or B, as well as the integer coefficients of the two basis vectors dtimes(pm 1 hatx + sqrt3haty)  2 (u for +, v for -), with d = 246Å as the lattice constant.\n\n\n\n\n\n","category":"type"},{"location":"api.html#GrapheneQFT.GrapheneState","page":"API","title":"GrapheneQFT.GrapheneState","text":"GrapheneState(coord::GrapheneCoord, spin::Spin)\n\nQuantum state of an electron in graphene, denoted by  u v Lrangleotimes sigmarangle in the drivation.\n\nThe state is given by the GrapheneCoord of the orbital, as well as the electronic spin, which can take values SpinUp and SpinDown.\n\n\n\n\n\n","category":"type"},{"location":"api.html#GrapheneQFT.Defect","page":"API","title":"GrapheneQFT.Defect","text":"Defect\n\nAbstract type for all defect types.\n\n\n\n\n\n","category":"type"},{"location":"api.html#GrapheneQFT.ImpurityState","page":"API","title":"GrapheneQFT.ImpurityState","text":"ImpurityState(ϵ::Float64, coupling::Vector{Tuple{Float64,GrapheneCoord}}) <: Defect\n\nAn impurity state of energy ϵ (in eV) coupled to the graphene system. The tuples in the coupling field contain all the coupling energies (in eV) and the  corresponding GrapheneCoord's.\n\n\n\n\n\n","category":"type"},{"location":"api.html#GrapheneQFT.LocalSpin","page":"API","title":"GrapheneQFT.LocalSpin","text":"LocalSpin(x::Float64, y::Float64, z::Float64, coord::GrapheneCoord) <: Defect\n\nA local spin with components x, y, and z located at coord. See GrapheneCoord.\n\n\n\n\n\n","category":"type"},{"location":"api.html#GrapheneQFT.Hopping","page":"API","title":"GrapheneQFT.Hopping","text":"Hopping(c1::GrapheneCoord, c2::GrapheneCoord, Δ::ComplexF64) <: Defect\n\nHopping modification (in eV) between two GrapheneCoord's. If c1==c2,  this quantity corresponds to a local energy modification.\n\n\n\n\n\n","category":"type"},{"location":"api.html#GrapheneQFT.GrapheneSystem","page":"API","title":"GrapheneQFT.GrapheneSystem","text":"GrapheneSystem(\n    μ::Float64,\n    T::Float64,\n    Δ::Array{ComplexF64,2},\n    V::Array{Float64,2},\n    scattering_states::Vector{GrapheneState},\n    imps::Vector{Float64},\n)\n\nA structure describing the perturbed graphene system.\n\nSee also mkGrapheneSystem.\n\n\n\n\n\n","category":"type"},{"location":"api.html#GrapheneQFT.graphene_neighbors","page":"API","title":"GrapheneQFT.graphene_neighbors","text":"graphene_neighbors(atom::GrapheneCoord)\n\nDetermine the nearest neighbors of an atom and return a vector of the corresponding GrapheneCoord's.\n\n\n\n\n\n","category":"function"},{"location":"api.html#GrapheneQFT.graphene_multiple_neighbors","page":"API","title":"GrapheneQFT.graphene_multiple_neighbors","text":"graphene_multiple_neighbors(atom::GrapheneCoord, n::Int)\n\nReturn an array of GrapheneCoord's obtained by iteratively running graphene_neighbors n times, applying it only to the newly-added GrapheneCoord's from the past iteration. The entries in the result are unique and sorted by their distance from atom.\n\nArguments\n\natom: GrapheneCoord from which the iteration begins\nn: number of iterations of graphene_neighbors.\n\n\n\n\n\n","category":"function"},{"location":"api.html#GrapheneQFT.crystal_to_cartesian","page":"API","title":"GrapheneQFT.crystal_to_cartesian","text":"crystal_to_cartesian(coord::GrapheneCoord)\n\nConvert a GrapheneCoord to a cartesian point with lengths in Å, where GrapheneCoord(0, 0, A) is at the origin.\n\n\n\n\n\n","category":"function"},{"location":"api.html#GrapheneQFT.mkGrapheneSystem","page":"API","title":"GrapheneQFT.mkGrapheneSystem","text":"mkGrapheneSystem(\n    μ::Float64,\n    T::Float64,\n    defects::Vector{Defect},\n)\n\nConstruct GrapheneSystem.\n\nArguments\n\nμ: chemical potential\nT: temperature\ndefects: a list of [Defect]'s\n\nWhen supplying []Hopping](@ref) in defects, ensure that each coordinate pair appears only once as repeated pairs with different couplings will cause eariler values to be overwritten. The order of the coordinates does not matter.\n\nThe function constructs a GrapheneSystem with the Δ and V matrices.  In addition, a list of all GrapheneState that are perturbed (scattering_states field in GrapheneCoord) and a list of impurity energies (imps field in GrapheneCoord) are included.\n\n\n\n\n\n","category":"function"},{"location":"api.html#GrapheneQFT.peierls_phase","page":"API","title":"GrapheneQFT.peierls_phase","text":"peierls_phase(vec_pot, a1::GrapheneCoord, a2::GrapheneCoord)\n\nCalculate the phase used in the Peierls substitution to include the effects of the magnetic field.\n\nTo make the units work out better, the magnetic field mathbfB(x y z) = Phi_0mathcalV mathbff(x y z), where mathcalV is the area of the graphene unit cell in Å², Phi_0 = h  2e is the magnetic flux quantum, and mathbff(x y z) = nabla times mathbfg(xyz) is a dimensionless vector function. Note that mathbff = 1 produces a field of about 40000 T.\n\nFor the vector potential, we have mathbfA(x y z) =Phi_0  mathcalV mathbfg(x y z), where mathbfg(x y z) has the units of Å. Using the definition of the Peierls phase, one gets\n\nphi = -fracpi Phi_0 int mathbfAcdot dmathbfl =\n-fracpi mathcalV int  mathbfg(x y z) cdot dmathbfl\nrightarrow\n-fracpi mathcalV int  mathbfg_xy(x y) cdot dmathbfl\n\nThe last step follows from the fact that the graphene system resides in the xy plane, so one needs to retain only the x and y componends of mathbfg, as denoted by the subscript xy.\n\nArguments\n\nvec_pot(x,y): mathbfg_xy(x y) in Å with x and y in Å.\na1: GrapheneCoord of the \"from\" atom.\na2: GrapheneCoord of the \"to\" atom.\n\nvec_pot(x,y) needs to return a tuple corresponding to the vector potential in     x and y directions.\n\n\n\n\n\n","category":"function"},{"location":"api.html#Green's-Functions","page":"API","title":"Green's Functions","text":"","category":"section"},{"location":"api.html","page":"API","title":"API","text":"  δG_R\n  G_R\n  δΓ\n  Γ","category":"page"},{"location":"api.html#GrapheneQFT.δG_R","page":"API","title":"GrapheneQFT.δG_R","text":"δG_R(z::ComplexF64,\n pairs::Vector{Tuple{GrapheneState,GrapheneState}},\n s::GrapheneSystem)\n\nThe correction to the real-space graphene Green's function in the presence of defects as a function of complex energy z.\n\nThe function returns a vector of ComplexF64 for each pair of     GrapheneState's in pairs.\n\nArguments\n\nz: complex energy\npairs: pairs of GrapheneState's for which δG_R is calculated\ns: GrapheneSystem for which δG_R is calculated\n\n\n\n\n\n","category":"function"},{"location":"api.html#GrapheneQFT.G_R","page":"API","title":"GrapheneQFT.G_R","text":"G_R(z::ComplexF64,\npairs::Vector{Tuple{GrapheneState,GrapheneState}},\ns::GrapheneSystem)\n\nThe full real-space graphene Green's function in the presence of defects as a function of complex energy z.\n\nArguments\n\nz: complex energy\npairs: pairs of GrapheneState's for which G_R is calculated\ns: GrapheneSystem for which G_R is calculated\n\n\n\n\n\n","category":"function"},{"location":"api.html#GrapheneQFT.δΓ","page":"API","title":"GrapheneQFT.δΓ","text":"δΓ(z::ComplexF64, s::GrapheneSystem)\n\nThe correction to the impurity Green's function due to the impurities' interaction with graphene.\n\nArguments\n\nz: complex energy\ns: GrapheneSystem for which δΓ is calculated\n\n\n\n\n\n","category":"function"},{"location":"api.html#GrapheneQFT.Γ","page":"API","title":"GrapheneQFT.Γ","text":"Γ(z::ComplexF64, s::GrapheneSystem)\n\nThe full impurity Green's function with the correction due to the impurities' interaction with graphene.\n\nArguments\n\nz: complex energy\ns: GrapheneSystem for which Γ is calculated\n\n\n\n\n\n","category":"function"},{"location":"index.html#GrapheneQFT.jl","page":"Home","title":"GrapheneQFT.jl","text":"","category":"section"},{"location":"index.html","page":"Home","title":"Home","text":"  GrapheneQFT","category":"page"},{"location":"index.html#GrapheneQFT","page":"Home","title":"GrapheneQFT","text":"This package provides provides a set of functions to facilitate the field-theoretic treatment of monolayer graphene using the tight-binding model. The Hamiltonian employed by this package includes only the nearest-neighbor hopping term with t = 28eV. The derivation of the formalism is available here.\n\n\n\n\n\n","category":"module"}]
}
