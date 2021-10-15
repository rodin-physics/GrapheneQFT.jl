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
    +\left(\bar{\boldsymbol{\psi}}_n V\boldsymbol{\phi}_n + \bar{\boldsymbol{\phi}}_n V^\dagger \boldsymbol{\psi}_n\right) + \bar{\boldsymbol{\phi}}_n \overbrace{\left(-i\omega_n - \mu + h\right)}^{-\Gamma_{i\omega_n+\mu}^{-1}}\boldsymbol{\phi}_n
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
We can identify the full system Green's function from the determinant in $\mathcal{Z}$

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
is the full Green's function of the bulk system including the effects of defects and perturbations. Similarly, the bottom right block in $\mathbf{G}_z$ corresponds to the full Green's function of the impurity states including their coupling to the perturbed bulk system:

$$\Gamma_z + \Gamma_z V^\dagger \mathcal{G}_z V \Gamma_z
= \Gamma_z + \Gamma_z V^\dagger G_z \left[1 - (\Delta + V \Gamma_z V^\dagger) G_z \right]^{-1} V \Gamma_z \, .$$
From the partition function, we can also write down the Helmholtz free energy $F = -\beta^{-1} \ln \mathcal{Z}$:

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


## Graphene

Having performed the required manipulations for a general system, we now focus on graphene. Because the electronic properties of graphene are dominated by the carbon $\pi$ orbitals, the electronic states in a pristine system can be described by $|\mathbf{r},L\rangle\otimes|\sigma\rangle$, where $\mathbf{r}$ is the coordinate of the unit cell hosting the orbital, $L$ is the sublattice of the atom, and $\sigma$ is the spin of the electron. This is the basis corresponding to the operators $c_k$ in the derivation above.

The next required ingredient is the graphene Green's function $G_z$. Since $G_z = (z - H)^{-1}$, it is a matrix whose elements give the propagation amplitudes between graphene states. The matrix elements can be calculated using

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

Using the multiples of the basis vectors to describe the electronic states ($|\mathbf{r},L\rangle\otimes|\sigma\rangle\rightarrow |u, v, L\rangle\otimes|\sigma\rangle$), the matrix elements $\langle u, v, L|\otimes\langle\sigma|(z - \hat{H})^{-1}|u', v' ,L'\rangle\otimes|\sigma'\rangle$ of the Green's function become

$$\langle L|\begin{pmatrix}
		z\Omega^{u-u',v-v'}_z
		&
		- t\left[\Omega^{u-u',v-v'}_z + \Omega^{u-u',v-v'}_{+,z} \right]
		\\
		- t\left[\Omega^{u-u',v-v'}_z + \Omega^{u-u',v-v'}_{-,z}\right]
		&
		z\Omega^{u-u',v-v'}_z
	\end{pmatrix}|L'\rangle \delta_{\sigma\sigma'}.$$

Depending on the sublattices of the two states, we pick out the appropriate element of the $2\times 2$ matrix above.

### Graphene Green's Function

Following the general derivation above, the individual elements of the full real-space Green's function are

$$\mathcal{G}_z^{jk}
=  G_z^{jk}+\sum_{lm}G_z^{jl}\left\{\left[ \left(\Delta +  V\Gamma_zV^\dagger\right)^{-1}-
G_z\right]^{-1}\right\}_{lm}G_z^{mk}\,,$$

where $j$, $k$, $l$, and $m$ are compact labels for the graphene states $|u, v ,L\rangle\otimes|\sigma\rangle$. As written, the sum over $l$ and $m$ includes all the states in the system. However, the expression can be made considerably simpler. The term $\Delta + V\Gamma_z V^\dagger$ is a square matrix including all the system states. Let us rearrange the elements in the matrix to give it a block diagonal form:

$$\Delta + V\Gamma_zV^\dagger \rightarrow \begin{pmatrix}0&0\\0&A\end{pmatrix}\,,$$

where we created a square matrix $A$ that contains only the perturbed states. Consequently,

$$\left(\Delta + V\Gamma_zV^\dagger\right)^{-1} \rightarrow \begin{pmatrix}\infty &0\\0&A^{-1}\end{pmatrix}$$

and

$$\left[ \left(\Delta +  V\Gamma_zV^\dagger\right)^{-1}-
G_z\right]^{-1} \rightarrow \begin{pmatrix}0 &0\\0& \left[A^{-1} - \tilde{G}_z\right]^{-1}\end{pmatrix}\,,$$

where $\tilde{G}_z$ is a portion of $G_z$ that only includes the perturbed states. Hence, the $lm$ summation only needs to run over these, too:

$$\mathcal{G}_z^{jk}
=  G_z^{jk}+\sum_{lm\in\mathrm{pert}}G_z^{jl}\left\{\left[ \left(\tilde{\Delta} +  \tilde{V}\Gamma_z\tilde{V}^\dagger\right)^{-1}-
\tilde{G}_z\right]^{-1}\right\}_{lm}G_z^{mk}\,.$$

Here, the tilde over the matrices indicates that only the elements corresponding to the perturbed states are retained. Note that a particular state is included in $\tilde{V}$ *and* $\tilde{\Delta}$ even if it is perturbed by only one of the terms!

### Impurity Green's Function

For the full impurity Green's function, we follow the discussion above to obtain

$$\Gamma_z + \Gamma_z \tilde{V}^\dagger \tilde{G}_z \left[1 - (\tilde{\Delta} + \tilde{V} \Gamma_z \tilde{V}^\dagger) \tilde{G}_z \right]^{-1} \tilde{V} \Gamma_z \, .$$
