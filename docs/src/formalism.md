# Formalism

## Hamiltonian

To model graphene coupled to external electronic states and influenced by an external perturbation, it is helpful to start with the following generic Hamiltonian

$$\hat{\mathcal{H}} = \sum_{ab} c^\dagger_a  H_{ab} c_b + \sum_{ab} g^\dagger_a h_{ab} g_b+ \sum_{ab} g^\dagger_a V_{ab} c_b + c^\dagger_b V_{ab}^* g_a
	\,.$$

Here, we have introduced two fermionic systems with the corresponding second-quatized operators $c$ and $g$, governed by the Hamiltonians $H$ and $h$, respectively, while $V_{ab}$ represent the coupling between them. The subscripts $a$ and $b$ label the individual states in these systems.

The expression can be made more compact by defining $\mathbf{c}$ and $\mathbf{d}$ as column vectors of $c_a$ and $d_a$:

$$\hat{\mathcal{H}} = 	\mathbf{c}^\dagger H\mathbf{c}
	+ 	\mathbf{g}^\dagger h\mathbf{g}
	+ \mathbf{g}^\dagger V \mathbf{c} +  \mathbf{c}^\dagger V^\dagger \mathbf{g}
	\,.$$

At this point, we identify $h$ as the Hamiltonian describing the external electronic states and $H = H_0 + \Delta$ as a combination of the pristine graphene Hamiltonian $H_0$ and a perturbation $\Delta$.

Pristine graphene Hamiltonian can be written using the tight-binding formalism in momentum space as

$$\hat{H}_0 = \sum_\mathbf{q} \begin{pmatrix}a_\mathbf{q}^\dagger & b_\mathbf{q}^\dagger \end{pmatrix} \begin{pmatrix} 0 & -t f_\mathbf{q} \\ -t f^*_\mathbf{q} & 0\end{pmatrix} \begin{pmatrix}a_\mathbf{q}\\ b_\mathbf{q}\end{pmatrix} = \sum_\mathbf{q}c^\dagger_\mathbf{q}H_\mathbf{q} c_\mathbf{q}\,,$$

where $a_\mathbf{q}$ and $b_\mathbf{q}$ are the annihilation operators for the $p_z$ orbitals of the two sublattices, $t = 2.8$eV is the nearest-neighbor hopping integral, $f_\mathbf{q} = 1 + e^{i\mathbf{q}\cdot\mathbf{d}_1}+e^{i\mathbf{q}\cdot\mathbf{d}_2}$, and $\mathbf{d}_{1/2}=d\left(\pm1, \sqrt{3}\right)/2$ are the lattice vectors. On the other hand, coupling to the external states and the perturbation $\Delta$ are easier to write in real space. Hence, using $c^\dagger_\mathbf{q} = N^{-1/2}\sum_\mathbf{R}c_\mathbf{R}^\dagger e^{i\mathbf{R}\cdot\mathbf{q}} = \mathbf{c}^\dagger \Theta_\mathbf{q}^\dagger$, where $N$ is the number of unit cells in the system and $\Theta_\mathbf{q}^\dagger$ is a column vector of $e^{i\mathbf{R}\cdot\mathbf{q}}$, we have

$$\hat{H}_0 =\mathbf{c}^\dagger \left(   \sum_\mathbf{q}  \Theta_\mathbf{q}^\dagger H_\mathbf{q} \Theta_\mathbf{q} \right)\mathbf{c} =\mathbf{c}^\dagger \Theta^\dagger H_0^\mathbf{Q} \Theta \mathbf{c} = \mathbf{c}^\dagger H_0 \mathbf{c}\,.$$

Here, $\mathbf{c}$ is a column vector of $c_\mathbf{R}$, $\Theta$ is a column vector of $\Theta_\mathbf{q}$, and $H_0^\mathbf{Q}$ is a block-diagonal matrix of $H_\mathbf{q}$.


## Partition Function

Because $\hat{\mathcal{H}}$ is normal-ordered, we can translate it into the imaginary-time
action

$$S = \sum_{\omega_n}	\bar\psi_n \left[\left(-i\omega_n  - \mu\right)\mathbf{1}+ H\right]\psi_n
	+ \bar\phi_n\left[\left(-i\omega_n  - \mu\right)\mathbf{1}+ h \right]\phi_n
	+ \bar\phi_nV \psi_n +  \bar\psi_n V^\dagger \phi_n
	\\
	= \sum_{\omega_n}	\bar\psi_n \left(-G^{-1}_{i\omega_n+\mu}\right)\psi_n
		+ \bar\phi_n\left(-\Gamma^{-1}_{0,i\omega_n+\mu} \right)\phi_n
		+ \bar\phi_nV \psi_n +  \bar\psi_n V^\dagger \phi_n
		\\
	= \sum_{\omega_n}
	\begin{pmatrix}
		\bar\psi_n & \bar\phi_n
	\end{pmatrix}	 
	\begin{pmatrix}
		-G^{-1}_{i\omega_n+\mu}	& V^\dagger
			\\
			V &-\Gamma^{-1}_{0,i\omega_n+\mu}
	\end{pmatrix}
	\begin{pmatrix}
		\psi_n \\ \phi_n
	\end{pmatrix}
	\,,$$

where $\omega_n$ are the fermionic Matsubara frequencies, $\mu$ is the chemical potential, and $\psi_n$ and $\phi_n$ are vectors of the Grassmann fields corresponding to $c$ and $d$ operators, respectively.

Exponentiating $-S$ and integrating over all the Grassmann fields yields the partition function

$$\mathcal{Z} = \prod_{\omega_n} \det \left[\beta\begin{pmatrix}
	-G^{-1}_{i\omega_n+\mu}	& V^\dagger
		\\
		V &-\Gamma^{-1}_{0,i\omega_n+\mu}
\end{pmatrix}\right]\,,$$

where $\beta$ is $1 / (k_B T)$.

## Green's Functions

The partition function can be also written as

$$\mathcal{Z} =
	\prod_{\omega_n }
	\left|-\beta \Gamma^{-1}_{0,i\omega_n+\mu}\right|
	\left|-\beta\left( G_{i\omega_n+\mu} ^{-1} -  V\Gamma_{0,i\omega_n+\mu}V^\dagger\right)\right|
		\\
		=
    \prod_{\omega_n }
    \left|-\beta G^{-1}_{i\omega_n+\mu}\right|
    \left|-\beta\left(\Gamma^{-1}_{0,i\omega_n+\mu} - V^\dagger G_{i\omega_n+\mu}V\right)\right|\,.$$

We identify the quantity in the parentheses of the first line as the inverse of the  full graphene Green's function (including the effects of the external states and a perturbation). The quantity in the parentheses of the second line is the inverse of the full Green's function of the external states.

### Graphene Green's Function

Defining the pristine graphene Green's function in real space as $G_{0,z} = \left(z -  H_0 \right)^{-1}$ allows us write the full real-space Green's function as $\mathcal{G}_z^{-1} =  G_{0,z}^{-1} - \Delta -  V\Gamma_{0,z}V^\dagger$. Inverting this expression yields

$$\mathcal{G}_z =  \left[G_{0,z}^{-1} - \Delta -  V\Gamma_{0,z}V^\dagger\right]^{-1}
=  \left[1 - G_{0,z}\left(\Delta +  V\Gamma_{0,z}V^\dagger\right)\right]^{-1}G_{0,z}
\\
=  G_{0,z}+G_{0,z}\left(\Delta +  V\Gamma_{0,z}V^\dagger\right)\left[ 1-
G_{0,z}\left(\Delta +  V\Gamma_{0,z}V^\dagger\right)\right]^{-1}G_{0,z}
\\
=  G_{0,z}+G_{0,z}\left[ \left(\Delta +  V\Gamma_{0,z}V^\dagger\right)^{-1}-
G_{0,z}\right]^{-1}G_{0,z}\,.$$

We can obtain $G_{0,z}$ from the momentum-space graphene Hamiltonian as follows:

$$G_{0,z} = \left(z -  H_0 \right)^{-1} = \left(z - \Theta^\dagger H_0^\mathbf{Q}\Theta \right)^{-1} = \Theta\left(z -  H_0^\mathbf{Q} \right)^{-1}\Theta^\dagger\,.$$

The $2\times 2$ elements of $G_{0,z}$ are computed from

$$G_{0,z,2\times 2}^{jk} = \sum_{lm}\Theta_{jl}\left[\left(z -  H_0^\mathbf{Q} \right)^{-1}\right]_{lm}\Theta^\dagger_{mk}
\\
\rightarrow \sum_{\mathbf{q}}\Theta_{\mathbf{R}_j\mathbf{q}}\left[\left(z -  H_0^\mathbf{Q} \right)^{-1}\right]_\mathbf{q}\Theta^\dagger_{\mathbf{q}\mathbf{R}_k} = \frac{1}{N}\sum_\mathbf{q} \left[\left(z -  H_0^\mathbf{Q} \right)^{-1}\right]_\mathbf{q} e^{i\left(\mathbf{R}_k - \mathbf{R}_j\right)\cdot\mathbf{q}}
\\
= \frac{1}{N}\sum_\mathbf{q}
\begin{pmatrix}
z& - tf_\mathbf{q}
\\
- tf_\mathbf{q}^* & z
\end{pmatrix}\frac{1}{z^2 - t^2 \left|f_\mathbf{q}\right|} e^{i\mathbf{R}_{kj}\cdot\mathbf{q}}= G_{0,z}\left(\mathbf{R}_{kj}\right)\,.$$

To compute $G_{0,z}(\mathbf{R})$, we first introduce

$$\Omega^{u,v}_z =
	\frac{1}{N}\sum_{\mathbf{q}\in\mathrm{BZ}}
	\frac{
		e^{i\mathbf{q}\cdot \left(u\mathbf{d}_1 + v\mathbf{d}_2\right)}
	}
	{z^2 - t^2\left| f_{\mathbf{q}}\right|^2}$$

with $u\mathbf{d}_1 + v\mathbf{d}_2 = \frac{d}{2}\left(u - v, \sqrt{3}\left(u+v\right)\right)$ and $t = 2.8$eV as the nearest-neighbor hopping energy. Using $\mathbf{q}\cdot \left(u\mathbf{d}_1 + v\mathbf{d}_2\right)  = \frac{d}{2}\left[\left(u - v\right)q_x + \sqrt{3}\left(u+v\right)q_y\right]$ and turning the momentum sum into an integral yields

$$\Omega^{u,v}\left(z\right)
	= \frac{1}{\left(2\pi\right)^2}\oint dx \oint dy
	\frac
	{e^{i \left[\left(u - v\right)x + \left(u+v\right)y\right]}}
	{z^2 - t^2\left(1 + 4\cos^2 x + 4 \cos x\cos y \right)}\,.$$

From

$$\oint d\theta \frac{e^{il\theta}}{W-\cos\theta} = 2\pi \frac{\left(W - \sqrt{W - 1}\sqrt{W + 1}\right)^{|l|}}{\sqrt{W - 1}\sqrt{W + 1}}\,,$$

we get

$$\Omega^{u,v}_z = \frac{1}{2\pi}\frac{1}{4t^2}
	\oint dx \frac{e^{i\left(u - v\right)x}}{\cos x}\frac{\left(W - \sqrt{W - 1}\sqrt{W + 1}\right)^{|u+v|}}{\sqrt{W - 1}\sqrt{W + 1}}\,,
	\\
	W = \frac{\frac{z^2}{t^2}-1}{4\cos x}-\cos x\,.$$

Finally, $G_{0,z}({\mathbf{R}})$ for $\mathbf{R} = u\mathbf{d}_1 + v\mathbf{d}_2$ can be written as

$$G_{0,z}({\mathbf{R}})
	=
	\begin{pmatrix}
		z\Omega^{u,v}_z
		&
		- t\left[\Omega^{u,v}_z + \Omega^{u,v}_{+,z} \right]
		\\
		- t\left[\Omega^{u,v}_z + \Omega^{u,v}_{-,z}\right]
		&
		z\Omega^{u,v}_z
	\end{pmatrix}\,,
	\\
	\Omega^{u,v}_{\pm,z}
	=
	 \frac{1}{2\pi}\frac{1}{4t^2}
	\oint dx \,2e^{i\left(u - v\right)x}\frac{\left(W - \sqrt{W - 1}\sqrt{W + 1}\right)^{|u+v\pm 1|}}{\sqrt{W - 1}\sqrt{W + 1}}\,.$$

The individual elements of the full real-space Green's function are

$$\mathcal{G}_z^{jk}
=  G_{0,z}^{jk}+\sum_{lm}G_{0,z}^{jl}\left\{\left[ \left(\Delta +  V\Gamma_{0,z}V^\dagger\right)^{-1}-
G_{0,z}\right]^{-1}\right\}_{lm}G_{0,z}^{mk}\,,$$

where $G_{0,z}^{jk}$ are the appropriate elements of the $2\times 2$ $G_{0,z}(\mathbf{R}_{kj})$ blocks.

As written, the sum over $l$ and $m$ includes all the atoms in the system. However, the expression can be made considerably simpler. The term $\Delta + V\Gamma_{0,z}V^\dagger$ is a square matrix including all the system atoms. Let us rearrange the elements in the matrix to give it a block diagonal form:

$$\Delta + V\Gamma_{0,z}V^\dagger \rightarrow \begin{pmatrix}0&0\\0&F\end{pmatrix}\,,$$

where we created a square matrix $F$ that contains only the perturbed atoms. Consequently,

$$\left(\Delta + V\Gamma_{0,z}V^\dagger\right)^{-1} \rightarrow \begin{pmatrix}\infty &0\\0&F^{-1}\end{pmatrix}$$

and

$$\left[ \left(\Delta +  V\Gamma_{0,z}V^\dagger\right)^{-1}-
G_{0,z}\right]^{-1} \rightarrow \begin{pmatrix}0 &0\\0& \left[F^{-1} - \tilde{G}_{0,z}\right]^{-1}\end{pmatrix}\,,$$

where $\tilde{G}_{0,z}$ is a portion of $G_{0,z}$ that only includes the perturbed atoms. Hence, the $lm$ summation only needs to run over the perturbed atoms, too:

$$\mathcal{G}_z^{jk}
=  G_{0,z}^{jk}+\sum_{lm\in\mathrm{pert}}G_{0,z}^{jl}\left\{\left[ \left(\tilde{\Delta} +  \tilde{V}\Gamma_{0,z}\tilde{V}^\dagger\right)^{-1}-
\tilde{G}_{0,z}\right]^{-1}\right\}_{lm}G_{0,z}^{mk}\,.$$

Here, the tilde over the matrices indicates that only the elements corresponding to the perturbed atoms are retained. Note that a particular atom is included in $\tilde{V}$ *and* $\tilde{\Delta}$ even if it is perturbed by only one of the terms!

### Impurity Green's Function

To calculate the full impurity Green's function, denoted by $\Gamma_z$, we write

$$\Gamma_z = \left[\Gamma^{-1}_{0,z} - V^\dagger \left(z - H_0 - \Delta\right)^{-1}V\right]^{-1} = \left[1 - \Gamma_{0,z}V^\dagger \left(z - H_0 - \Delta\right)^{-1}V\right]^{-1}\Gamma_{0,z}
\\
= \Gamma_{0,z} + \Gamma_{0,z}V^\dagger\Lambda_z\left(1-V\Gamma_{0,z}V^\dagger \Lambda_z\right)^{-1}V\Gamma_{0,z}\,,$$

where

$$\Lambda_z = \left(z - H_0 - \Delta\right)^{-1} = \left[1 - \left(z - H_0\right)^{-1}\Delta\right]^{-1} \left(z - H_0\right)^{-1} =\left(1 - G_{0,z}\Delta\right)^{-1} G_{0,z}\,.$$

Similarly to the discussion above, we only need to include the graphene atoms that are perturbed in the $\Lambda$ and $V$ matrices. This yields

$$\Gamma_z
= \Gamma_{0,z} + \Gamma_{0,z}\tilde{V}^\dagger\tilde{\Lambda}_z\left(1-\tilde{V}\Gamma_{0,z}\tilde{V}^\dagger \tilde{\Lambda}_z\right)^{-1}\tilde{V}\Gamma_{0,z}\,,$$

where

$$\tilde{\Lambda}_z =\left(1 - \tilde{G}_{0,z}\tilde{\Delta}\right)^{-1}\tilde{G}_{0,z}\,.$$

## Free Energy




## Occupation Number
