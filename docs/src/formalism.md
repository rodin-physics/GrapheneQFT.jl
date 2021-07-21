# Formalism

## Two-system Hamiltonian

To model graphene coupled to external electronic states and influenced by an external perturbation, it is helpful to start with the following generic Hamiltonian

$$\hat{\mathcal{H}} = \sum_{ab} c^\dagger_a H_{ab} +\Delta_{ab}c_b + \sum_{ab} d^\dagger_a h_{ab} d_b+ \sum_{ab} d^\dagger_a V_{ab} c_b + c^\dagger_b V_{ab}^* d_a
	\,.$$

Here, we have introduced two fermionic systems with the corresponding second-quatized operators $c$ and $d$, governed by the Hamiltonians $H$ and $h$, respectively. The subscripts $a$ and $b$ label the individual states in these systems, while $V_{ab}$ represent the coupling between them.

The expression can be made more compact by defining $\mathbf{c}$ and $\mathbf{d}$ as column vectors of $c_a$ and $d_a$, respectively:

$$\hat{\mathcal{H}} = 	\mathbf{c}^\dagger H\mathbf{c}
	+ 	\mathbf{d}^\dagger h\mathbf{d}
	+ \mathbf{d}^\dagger V \mathbf{c} +  \mathbf{c}^\dagger V^\dagger \mathbf{d}
	\,.$$

At this point, we identify $h$ as the Hamiltonian describing the external electronic states and $H = H_0 + \Delta$ as a combination of the pristine graphene Hamiltonian $H_0$ and a perturbation $\Delta$. Below, we will use these Hamiltonians to extract the relevant Green's functions and compute the Helmholtz free energy.

## Partition Function

Because $\hat{\mathcal{H}}$ is normal-ordered, we can translate it into the imaginary-time
action

$$S = \sum_{\omega_n}	\bar\psi_n \left[\left(-i\omega_n  - \mu\right)\mathbf{1}+ H\right]\psi_n
	+ \bar\phi_n\left[\left(-i\omega_n  - \mu\right)\mathbf{1}+ h \right]\phi_n
	+ \bar\phi_nV \psi_n +  \bar\psi_n V^\dagger \phi_n
	\\
	= \sum_{\omega_n}	\bar\psi_n \left(-G^{-1}_n\right)\psi_n
		+ \bar\phi_n\left(-\Gamma^{-1}_n \right)\phi_n
		+ \bar\phi_nV \psi_n +  \bar\psi_n V^\dagger \phi_n
		\\
	= \sum_{\omega_n}
	\begin{pmatrix}
		\bar\psi_n & \bar\phi_n
	\end{pmatrix}	 
	\begin{pmatrix}
		-G^{-1}_n	& V^\dagger
			\\
			V &-\Gamma^{-1}_n
	\end{pmatrix}
	\begin{pmatrix}
		\psi_n \\ \phi_n
	\end{pmatrix}
	\,,$$

where $\omega_n$ are the fermionic Matsubara frequencies, $\psi_n$ and $\phi_n$ are vectors of the Grassmann fields corresponding to $c$ and $d$ operators, respectively.

Exponentiating $-S$ and integrating over all the Grassmann fields yields the partition function

$$\mathcal{Z} = \prod_n \det \left[\beta\begin{pmatrix}
	-G^{-1}_n	& V^\dagger
		\\
		V &-\Gamma^{-1}_n
\end{pmatrix}\right]\,.$$

## Green's Functions

## Free Energy

## Graphene Propagator
