# Formalism

## General System

Instead of starting directly with the problem of defects in graphene, we start with a more general scenario which will make the derivation more transparent. Consider a pristine system described by a second-quantized Hamiltonian $\sum_{jk} c^\dagger_j H_{jk}c_k$, where $j$ and $k$ label the states of the system, $c_j$ and $c_j^\dagger$ are the corresponding annihilation and creation operators, and $H_{jk}$ gives the coupling between the states. Next, we introduce two types of perturbation to this system. First, we can modify the coupling between the states. Second, we introduce a second system described by a Hamiltonian $\sum_{jk} g^\dagger_j h_{jk}g_k$ and couple it to the original one. The total Hamiltonian, then, can be written as

$$\hat{\mathcal{H}} = \sum_{jk} c^\dagger_j H_{jk}c_k 
    + 
    \sum_{jk} c^\dagger_j \Delta_{jk}c_k
    +
    \sum_{jk} \left(c^\dagger_j V_{jk}g_k + g^\dagger_k V_{jk}^* c_j\right)
    +
    \sum_{jk} g^\dagger_j h_{jk}g_k\,,$$
where the second term in the first line describes the modified coupling in the original system and the first term in the second line gives the coupling between the two systems. The Hamiltonian can be written more compactly if we assemble the operators and the coupling elements into vectors and matrices, respectively:

$$\hat{\mathcal{H}} = \mathbf{c}^\dagger H\mathbf{c} + \mathbf{c}^\dagger \Delta \mathbf{c}
    +\left(\mathbf{c}^\dagger V\mathbf{g} + \mathbf{g}^\dagger V^\dagger \mathbf{c}\right) + \mathbf{g}^\dagger h\mathbf{g}\,.$$

Next, we transcribe the Hamiltonian into the imaginary-time action

$$S = \sum_{n}\bar{\boldsymbol{\psi}}_n \overbrace{\left(-i\omega_n - \mu + H\right)}^{-G^{-1}_{i\omega_n + \mu}}\boldsymbol{\psi}_n + \bar{\boldsymbol{\psi}}_n\Delta \boldsymbol{\psi}_n
    +\left(\bar{\boldsymbol{\psi}}_n V\boldsymbol{\phi}_n + \bar{\boldsymbol{\phi}}_n V^\dagger \boldsymbol{\psi}\right) + \bar{\boldsymbol{\phi}}_n \overbrace{\left(-i\omega_n - \mu + h\right)}^{-\Gamma_{i\omega_n+\mu}^{-1}}\boldsymbol{\phi}_n
	\\
	=
	\sum_{n}
	\begin{pmatrix}
	\bar{\boldsymbol{\psi}}_n & \bar{\boldsymbol{\phi}}_n
	\end{pmatrix}
	\begin{pmatrix}
        -G_{i\omega_n+\mu}^{-1} + \Delta & V
        \\
        V^\dagger & -\Gamma_{i\omega_n+\mu}^{-1}
    \end{pmatrix}
	\begin{pmatrix}
	\boldsymbol{\psi}_n \\ \boldsymbol{\phi}_n
	\end{pmatrix}
	\,,$$

where $\omega_n$ are the fermionic Matsubara frequencies, $\mu$ is the chemical potential, and $\boldsymbol{\phi}_n$ and $\boldsymbol{\psi}_n$ ($\bar{\boldsymbol{\phi}}_n$ and $\bar{\boldsymbol{\psi}}_n$) are vectors of Grassmann numbers corresponding to $\mathbf{g}$ and $\mathbf{c}$ ($\mathbf{g}^\dagger$ and $\mathbf{c}^\dagger$). We also identify $G_{z}$ and $\Gamma_z$ as the Green's functions for the two systems.

Exponentiating $-S$ and integrating over all the Grassmann variables yields the partition function

$$\mathcal{Z} 
    =\prod_n
    \left|\beta
    \begin{pmatrix}
        -G_{i\omega_n+\mu}^{-1} + \Delta & V
        \\
        V^\dagger & -\Gamma_{i\omega_n+\mu}^{-1}
    \end{pmatrix}
    \right|\,,$$
where $\beta = 1 / (k_BT)$.

From this expression, we can write down the Helmholtz free energy $F = -\beta^{-1} \ln \mathcal{Z}$:

$$F = -\beta^{-1}\sum_n \ln\left|\beta
    \begin{pmatrix}
        -G_{i\omega_n+\mu}^{-1} + \Delta & V
        \\
        V^\dagger & -\Gamma_{i\omega_n+\mu}^{-1}
    \end{pmatrix}
    \right|
    \\
    = -\beta^{-1}\sum_n \ln\Bigg|\beta
    \begin{pmatrix}
        -G_{i\omega_n+\mu}^{-1} & 0
        \\
        0 & -\Gamma_{i\omega_n+\mu}^{-1}
    \end{pmatrix}
    \left[
    1
    +
    \begin{pmatrix}
        -G_{i\omega_n+\mu} & 0
        \\
        0 & -\Gamma_{i\omega_n+\mu}
    \end{pmatrix}
    \begin{pmatrix}
        \Delta & V
        \\
        V^\dagger & 0
    \end{pmatrix}
    \right]
    \Bigg|\,.$$

Removing the part of $F$ corresponding to the free energy of the two isolated system in the absence of any perturbation yields the defect- and coupling-induced modification to $F$

$$\delta F 
    = -\beta^{-1}\sum_n \ln\left|
    \begin{pmatrix}
       1- G_{i\omega_n+\mu}\Delta & -G_{i\omega_n+\mu}V
        \\
       - \Gamma_{i\omega_n+\mu}V^\dagger & 1
    \end{pmatrix}
    \right|
    \\
    = -\beta^{-1}\sum_n \ln\left|
      1- G_{i\omega_n+\mu}\left(\Delta+V\Gamma_{i\omega_n+\mu}V^\dagger\right)
    \right|\,.$$

We can also identify the full system Green's function from the determinant in $\mathcal{Z}$

$$\mathbf{G}_{z} 
    =
    \begin{pmatrix}
        G_{z}^{-1} - \Delta & -V
        \\
        -V^\dagger & \Gamma_{z}^{-1}
    \end{pmatrix}^{-1}
    =
    \begin{pmatrix}
         \mathcal{G}_z &\ \mathcal{G}_z V \Gamma_{z}
        \\
         \Gamma_{z}V^\dagger \mathcal{G}_z & \Gamma_{z}+ \Gamma_{z}V^\dagger \mathcal{G}_z V \Gamma_{z}
    \end{pmatrix}\,,$$
where

$$\mathcal{G}_z 
    =\left(G_{z}^{-1} - \Delta - V\Gamma_z V^\dagger\right)^{-1} 
    = G_{z}+G_{z}\left(\Delta + V\Gamma_z V^\dagger\right)
	\left[1-G_{z}\left(\Delta + V\Gamma_z V^\dagger\right)\right]^{-1}G_{z} \,,$$
is the full Green's function of the bulk system including the effects of defects and perturbations. Similarly, the bottom right block in $\mathbf{G}_z$ corresponds to the full Green's function of the impurity states including their coupling to the perturbed bulk system.

## Graphene Propagator

Having performed the required manipulations for a general system, we now focus on graphene. Because the electronic properties of graphene are dominated by the carbon $\pi$ orbitals, the electronic states in a pristine system can be described by $|\mathbf{r},L\rangle\otimes|\sigma\rangle$, where $\mathbf{r}$ is the coordinate of the unit cell hosting the orbital, $L$ is the sublattice of the atom, and $\sigma$ is the spin of the electron. This is the basis corresponding to the operators $c_k$ in the derivation above.

The next required ingredient is the graphene Green's function $G_z$. Since $G_z = (z - H)^{-1}$, it is a matrix whose elements give the propagation amplitudes between graphene states. The matrix elemnts can be calculated using

$$\langle\mathbf{r}, L|\otimes\langle\sigma|(z - \hat{H})^{-1}|\mathbf{r}',L'\rangle\otimes|\sigma'\rangle = \frac{1}{N}\sum_{\mathbf{qq}'}\langle\mathbf{q}, L| \otimes \langle \sigma|e^{i\mathbf{q}\cdot\mathbf{r}}(z - \hat{H})^{-1}e^{-i\mathbf{q}'\cdot\mathbf{r}'}|\mathbf{q}',L'\rangle\otimes|\sigma'\rangle,$$

where we Fourier-transformed the real-space states so that the sums over the momenta $\mathbf{q}$ and $\mathbf{q}'$ run over the entire Brillouin zone. Because the momentum-space Hamiltonian is diagonal in $\mathbf{q}$ and $\sigma$, we can write

$$\langle\mathbf{r}, L|\otimes\langle\sigma|(z - \hat{H})^{-1}|\mathbf{r}',L'\rangle\otimes|\sigma'\rangle = \frac{1}{N}\sum_{\mathbf{q}}\langle\mathbf{q}, L|e^{i\mathbf{q}\cdot\mathbf{r}}(z - \hat{H})^{-1}e^{-i\mathbf{q}\cdot\mathbf{r}'}|\mathbf{q},L'\rangle \delta_{\sigma\sigma'} 
\\ 
= \langle L|\frac{1}{N}\sum_{\mathbf{q}}e^{i\mathbf{q}\cdot\left(\mathbf{r}-\mathbf{r}'\right)}(z - H_\mathbf{q})^{-1}|L'\rangle \delta_{\sigma\sigma'}.$$

Using

$$H_\mathbf{q} =\begin{pmatrix} 0 & -t f_\mathbf{q} \\ -t f^*_\mathbf{q} & 0\end{pmatrix} ,$$

where $t = 2.8$eV is the nearest-neighbor hopping integral, $f_\mathbf{q} = 1 + e^{i\mathbf{q}\cdot\mathbf{d}_1}+e^{i\mathbf{q}\cdot\mathbf{d}_2}$, and $\mathbf{d}_{1/2}=d\left(\pm1, \sqrt{3}\right)/2$ are the lattice vectors, we have


$$\frac{1}{N}\sum_\mathbf{q} \left(z -  H_\mathbf{q} \right)^{-1} e^{i\left(\mathbf{r}_k - \mathbf{r}_j\right)\cdot\mathbf{q}}  = \frac{1}{N}\sum_\mathbf{q}
\begin{pmatrix}
z& - tf_\mathbf{q}
\\
- tf_\mathbf{q}^* & z
\end{pmatrix}\frac{1}{z^2 - t^2 \left|f_\mathbf{q}\right|^2} e^{i\mathbf{r}_{kj}\cdot\mathbf{q}}\,.$$

To perform the summation over $\mathbf{q}$, we first introduce

$$\Omega^{u,v}_z =
	\frac{1}{N}\sum_{\mathbf{q}\in\mathrm{BZ}}
	\frac{
		e^{i\mathbf{q}\cdot \left(u\mathbf{d}_1 + v\mathbf{d}_2\right)}
	}
	{z^2 - t^2\left| f_{\mathbf{q}}\right|^2}$$

with $u\mathbf{d}_1 + v\mathbf{d}_2 = \frac{d}{2}\left(u - v, \sqrt{3}\left(u+v\right)\right)$. Using $\mathbf{q}\cdot \left(u\mathbf{d}_1 + v\mathbf{d}_2\right)  = \frac{d}{2}\left[\left(u - v\right)q_x + \sqrt{3}\left(u+v\right)q_y\right]$ and turning the momentum sum into an integral yields

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

Finally, for $\mathbf{r} = u\mathbf{d}_1 + v\mathbf{d}_2$, we have

$$\frac{1}{N}\sum_\mathbf{q} \left(z -  H_\mathbf{q} \right)^{-1} e^{i\mathbf{r}\cdot\mathbf{q}}  
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

Using the multiples of the basis vectors to describe the electronic states ($|\mathbf{r},L\rangle\otimes|\sigma\rangle\rightarrow |u, v, L\rangle\otimes|\sigma\rangle$), the matrix elements of the Green's function become

$$\langle u, v, L|\otimes\langle\sigma|(z - \hat{H})^{-1}|u', v' ,L'\rangle\otimes|\sigma'\rangle 
\\
= \langle L|\begin{pmatrix}
		z\Omega^{u-u',v-v'}_z
		&
		- t\left[\Omega^{u-u',v-v'}_z + \Omega^{u-u',v-v'}_{+,z} \right]
		\\
		- t\left[\Omega^{u-u',v-v'}_z + \Omega^{u-u',v-v'}_{-,z}\right]
		&
		z\Omega^{u-u',v-v'}_z
	\end{pmatrix}|L'\rangle \delta_{\sigma\sigma'}.$$







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


Using the fact that $\mathbf{r}$ can be written as $u\mathbf{d}_1 + v\mathbf{d}_2$, where $\mathbf{d}_{1/2}=d\left(\pm1, \sqrt{3}\right)/2$ are the lattice vectors, the electronic basis becomes $|u,v,L\rangle\otimes|\sigma\rangle$. 

The next required ingredient is the graphene Green's function $G_z$. To obtain this quantity, first recall that  Using $a^\dagger_{\mathbf{q}\sigma} = N^{-1/2}\sum_\mathbf{R}a_{\mathbf{R}\sigma}^\dagger e^{i\mathbf{R}\cdot\mathbf{q}}$ (same for the other sublattice), where $N$ is the number of unit cells in the system, we have

$$\hat{H}_0 = \frac{1}{N}\sum_{\mathbf{RR}'}\sum_{\mathbf{qq}'}\sum_{\sigma\sigma'} \begin{pmatrix}a_{\mathbf{R}\sigma}^\dagger & b_{\mathbf{R}\sigma}^\dagger \end{pmatrix} \begin{pmatrix} 0 & -t f_\mathbf{q} \\ -t f^*_\mathbf{q} & 0\end{pmatrix} e^{i\mathbf{q}\cdot\mathbf{R}-i\mathbf{q}'\cdot\mathbf{R}'} \begin{pmatrix}a_{\mathbf{R}'\sigma'}\\ b_{\mathbf{R}'\sigma'}\end{pmatrix}\delta_{\mathbf{qq}'}\delta_{\sigma\sigma'}
\\
= \sum_{\mathbf{RR}'}\sum_{\sigma} \begin{pmatrix}a_{\mathbf{R}\sigma}^\dagger & b_{\mathbf{R}\sigma}^\dagger \end{pmatrix}\left[\frac{1}{N}\sum_{\mathbf{q}} \begin{pmatrix} 0 & -t f_\mathbf{q} \\ -t f^*_\mathbf{q} & 0\end{pmatrix} e^{i\mathbf{q}\cdot\left(\mathbf{R}-\mathbf{R}'\right)} \right]\begin{pmatrix}a_{\mathbf{R}'\sigma}\\ b_{\mathbf{R}'\sigma}\end{pmatrix}.$$

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


## Adding Spin Index
To include spin effects, we consider a similar formalism as above, with an additional spin term $J$ in the Hamiltonian.

$$\hat{\mathcal{H}} = 	\mathbf{c}^\dagger(H_0 + \Delta + J)\mathbf{c}
	+ 	\mathbf{g}^\dagger h\mathbf{g}
	+ \mathbf{g}^\dagger V \mathbf{c} +  \mathbf{c}^\dagger V^\dagger \mathbf{g} + H_{spin}
	\,.$$
We define $\Upsilon = \Delta + J$. The creation and annihilation operator vectors now include the spin index and we now write the pristine graphene Hamiltonian using the tight-binding formalism in momentum space as

$$\hat{H}_0 = \sum_\mathbf{q} \begin{pmatrix} a_{\uparrow, \mathbf{q}}^\dagger & b_{\uparrow, \mathbf{q}}^\dagger & a_{\downarrow, \mathbf{q}}^\dagger & b_{\downarrow, \mathbf{q}}^\dagger \end{pmatrix} \begin{pmatrix} 0 & -t f_\mathbf{q} & 0 & 0 \\ -t f^*_\mathbf{q} & 0 & 0 & 0 \\ 0 & 0 & 0 & -t f_{\mathbf{q}} \\ 0 & 0 -t f^*_{\mathbf{q}} \end{pmatrix} \begin{pmatrix} a_{\uparrow, \mathbf{q}} \\ b_{\uparrow, \mathbf{q}} \\ a_{\downarrow, \mathbf{q}} \\ b_{\downarrow, \mathbf{q}} \end{pmatrix} = \sum_\mathbf{q}c^\dagger_\mathbf{q}H_\mathbf{q} c_\mathbf{q}\,,$$

where all the terms are the same as before and $a_{\uparrow, \mathbf{q}}, a_{\downarrow, \mathbf{q}}$ and $b_{\uparrow, \mathbf{q}}, b_{\downarrow, \mathbf{q}}$ correspond to the different spins of the $p_z$ orbitals in the two sublattices. We can think of $J$ as an extension of the perturbation term $\Delta$. $\Delta$ allows for coupling between orbitals of the same spin; with the consideration of the spin index, $\Delta$ takes on a block-diagonal form. In the term $\Upsilon$, $J$ now accounts for the off-diagonal terms that allow coupling between orbitals of different spins.

Since the inclusion of the spin index does not change the form of the imaginary-time action, it follows that the partition function and by extension, the form of the full Green's functions will be almost identical.

### Graphene Green's Function

The full real-space graphene Green's function is given by $\mathcal{G}_z^{-1} = G_{0,z}^{-1} - \Upsilon - V \Gamma_{0,z}V^\dagger$. Inverting the expression gives us

$$\mathcal{G}_z = [G_{0,z}^{-1} - (\Upsilon + V \Gamma_{0,z} V^\dagger)]^{-1}
= G_{0,z} + G_{0,z} [(\Upsilon + V \Gamma_{0,z} V^\dagger)^{-1} - G_{0,z}]^{-1} G_{0,z} \\
= G_{0,z} + G_{0,z} [(\Delta + J + V \Gamma_{0,z} V^\dagger)^{-1} - G_{0,z}]^{-1} G_{0,z}\,.$$

### Impurity Green's Function

The full impurity Green's function, denoted by $\Gamma_z$, is given by

$$\Gamma_z = [\Gamma_{0,z}^{-1} - V^\dagger (z - H_0 - \Upsilon)^{-1} V]^{-1}
= \Gamma_{0,z} + \Gamma_{0,z} V^\dagger (\Lambda_z^{-1} - V \Gamma_{0,z} V^\dagger)^{-1} V \Gamma_{0,z}\,,$$

where, similar to the previous formalism,

$$\Lambda_z = (z - H_0 - \Upsilon)^{-1} = [1 - (z - H_0)^{-1} \Upsilon]^{-1} (z - H_0)^{-1} \\
= [1 - G_{0,z}(\Delta + J)]^{-1} G_{0,z} \, .$$

As discussed above, we only need to include atoms that are perturbed in the $\Lambda_z$ and/or the $V$ matrices. This yields

$$\Gamma_z = \Gamma_{0,z} + \Gamma_{0,z} \tilde{V}^\dagger (\tilde{\Lambda}_z^{-1} - \tilde{V} \Gamma_{0,z} \tilde{V}^\dagger)^{-1} \tilde{V} \Gamma_{0,z}\,,$$

where
$$\Lambda_z = [1 - \tilde{G}_{0,z}(\tilde{\Delta} + \tilde{J})]^{-1} \tilde{G}_{0,z} \, .$$

The term $\tilde{J}$ is a diagonal matrix composed of $\hat{J}_j$, where the index $j$ runs over perturbed atoms.

To work out the form of $\hat{J}_j$, we use the Ruderman–Kittel–Kasuya–Yosida (RKKY) formalism to obtain

$$\hat{J}_j = J_j \times \begin{pmatrix}
S_{jz} & 0 & S_{j+} & 0 \\
0 & S_{jz} & 0 & S_{j+} \\
S_{j-} & 0  & -S_{jz} & 0 \\
0 & S_{j-} & 0 & -S_{jz} \end{pmatrix} \, ,$$

where $S_{j\pm} = S_{jx} \pm i S_{jy}$.
