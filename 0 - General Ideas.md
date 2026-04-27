# Presentation
   Presentation takes 20 minutes, 18 for present, 2 for questions, sign up for one of three sessions.
   They'll be 1-4pm on April 29, 30, and May 1




# Ideas
- Visualizing Spin Liquid States in J1-J2 models
	- Visualizing Quantum Magnetic Phases in the 2D J1-J2 model
	- Doing this to see the differences between Hartree Fock, Hartree Fock Bogoliubov
- **Visualizations:**
	- Real Space Spin Texture
		- **Neel Phase:** $J_2/J_1 < 0.4$ Checkerboard Pattern
		- **Stripe Phase:** $J_2/J_1 > 0.6$ Stripes obviously
		- **QSL:** $0.4 < J_2/J_1 < 0.6$: No order Strange stuff
	- Magnetization
	- Pairing Strength (hfb only)
	- Energy surface as a function of symmetry-breaking parameter (pairing strength) to show how the minimum energy shifts as you project back to a fixed particle number
- **Simulation**:
	- **HF**: As standard, gets the $\rho$ matrix and eigenvectors, from which we can extract information
		- **Restricted**: Results are either Collinear in $\hat{z}$ or coplanar in the $\hat{x}\hat{y}$ plane (i think)
		- **Unrestricted**: Results are 3D vectors
		- **General**: Results are Complex Vectors (3D flags, I think) 
	- **HFB**: Allows breaking Number Symmetry by including a pairing channel. Gives $\rho$, $\kappa$, and quasiparticle eigenvectors


# Visualizations
###   Getting the Spin Vectors from $\rho$
   The density matrix $\rho_{ij}$ contains the local correlations. Spin arrows would be from diagonal elements like $$\langle S_i^a\rangle = \Tr{\rho_{ii}\sigma^a}\implies \vec{S}_i = (\Tr{\rho_{ii}\sigma^x},\ \Tr{\rho_{ii}\sigma^y},\ \Tr{\rho_{ii}\sigma^z})$$which gives you the 3D vector for a **Quiver Plot**. 

###   Getting Bond information from $\kappa$
   The pairing density $\kappa_{ij}$ tells how how "entangled" or "paired" two sites are in a way that breaks the $U(1)$ symmetry. $\kappa_{ij}$ is a complex phase where:
   - **Magnitude**: $|\kappa_{ij}|$ represents the strength of the "bond"
   - **Phase**: $\text{arg}(\kappa_{ij})$ represents the $U(1)$ phase. 
   We can plot these as bonds between sites in which:
   - **Thickness/Brightness** is magnitude $|\kappa_{ij}|$
   - **Color** is phase $\theta = \text{arg}(\kappa_{ij})$

###   Quasiparticles $(u,v)$
   When diagonalizing the HFB hamiltonian you get eigenvectors of $(u,v)$ which are the bogoliubov coefficients. These define quasiparticles $$\gamma_{k}^\dagger =\sum\limits_i\left(u_{ki}c_i^\dagger + v_{ki}c_i\right)$$We can visualize the "spatial profile" of a single quasiparticle by picking the eigenvector with the lowest energy and plotting $$|u_{ki}|^2 + |v_{ki}|^2$$as a heatmap across the lattice. This would display where a "hole" or "excitation" is likely to live in the frustrated system. Possibly revealing if the excitation is localized or delocalized.


# Evolutions
###   Parameter Sweep
   1) Choose parameter, starting and ending values, and step size
   2) Converge your solution for the starting parameter using random input matrices
   3) Alter the parameter, and use your previous solution as input to the next solution
   4) Repeat 3 until complete
###   Temperature (Finite-T HF/HFB)
   HF/HFB are **mean field theories**, their whole job is to average out the random fluctuations. Adding temperature will not give wiggling or fluctuating states over time, it will show melting. As temperature increase, the magnetic order will be destroyed. I could derive the curie temperature for my particular system and then show that magnetic order is destroyed there.
   We must calculate the **Thermal Occupation Number** $f_k$ $$\beta = \frac{1}{k_bT}\qquad f_k = \frac{1}{e^{\beta E_k} + 1}$$When we calculate the density matrix from our eigenvalues where we replace the standard $T=0$ formula by $$\rho = v^*v^T \to \rho_{ij} = \sum\limits_{k}\left[v_{ki}^*(1-f_k)v_{kj} + u^*_{ki}f_ku_{kj}\right]$$this will display like:
   - **Arrows (Magnetization)**: As $T$ goes up the thermal fluctuations destroy the magnetic order. 3D arrows will get shorter as the magnitude of $\langle S_i\rangle$ shrinks.
   - **Bonds (Pairing)**: The pairing gap $\Delta$ is extremely sensitive and as $T$ increases the colored bonds will get dimmer and thinner. *This temperature sensitive is why superconductors (many of which we believe use cooper pairs as their mechanism) can only exist at extremely low temperatures!*
   - **Critical Point $T_c$**: At a certain temperature we hit a phase transition at which the arrows and bonds completely disappear. This means we have a disordered paramagnet.

###   Time Domain (TDHF, TDHFB)
   Assuming HFB at this point, the Louiville-von Neumann equation ($\hbar = 1$) is $$\dd{\mathcal{R}}{t} = -i[\mathcal{H}(\mathcal{R}(t)), \mathcal{R}(t)]$$Since $\mathcal{H}$ depends on $\mathcal{R}$ this becomes a nonlinear ODE. We can simulate the dynamics using RK4 or Crank-Nicholson. The critical catch is that you have to rebuild $\mathcal{H}$ at every intermediate Runge-Kutta step. This where my developed superoperator method will do well. This method analyzes the hamiltonian, builds a list of all terms in the 2-body potential, then maps these terms with some index reordering into a large $N^2 \times N^2$ sparse matrix $\mathcal{G}$. Due to the particular nature of the index reordering we can simply apply a sparse matrix multiplication on a flattened $\rho$ or $\kappa$ to build the 2-body fock and pairing potentials $\Gamma$ and $\Delta$. This allows creation of the fock matrix in approximately $O(N)$ operations due to the extreme sparsity of the resulting supermatrix $\mathcal{G}$ for these nearest-neighbor lattice models. 
   I'd need derivative function $$f(\mathcal{R}) = -i[\mathcal{H}(\mathcal{R}), \mathcal{R}]$$ and my timesteps from $t\to t+dt$ would be
   1) Assign starting randomized $\mathcal{R}_n$
   2) Build $\mathcal{H}_1$ from $\mathcal{R}_n$
	   1) Calculate $k_1 = dt\cdot f(\mathcal{R}_n)$
   2) Form intermediate step $\mathcal{R}_{tmp} = \mathcal{R}_n + 0.5 k_1$
	   1) Build $\mathcal{H}_2(\mathcal{R}_{tmp})$ 
	   2) Calculate $k_2 = dt \cdot f(\mathcal{R}_{tmp})$
   3) Form $\mathcal{R}_{tmp2} = \mathcal{R}_n + 0.5 k_2$
	   1) Build $\mathcal{H}_3(\mathcal{R}_{tmp2})$ 
	   2) Calculate $k_3 = dt\cdot f(\mathcal{R}_{tmp2})$
   4) Form $\mathcal{R}_{tmp3} = \mathcal{R}_n + k_3$
	   1) Build $\mathcal{H}_4(\mathcal{R}_{tmp3})$
	   2) Calculate $k_4 = dt \cdot f(\mathcal{R}_{tmp3})$
   5) Take step $\mathcal{R}_{n+1} = \mathcal{R}_n + \frac{1}{6}(k_1 + 2k_2 + 2k_3 + k_4)$
   6) **Check** $\mathcal{R_n}$ for unitarity $\mathcal{R}^2 = \mathcal{R}$ and purify it. Monitor $\Tr{\mathcal{R}}$ to ensure we don't explode!
   7) Repeat 2-5

   This process for static parameters will simply arrive at a static solution and stay there. To see dynamics I'll need to:
   8) Find the static ground state $\mathcal{R}_0$ for specific parameters (all $J_1$ or all $J_2$)
   9) At $t=0$ abruptly change the parameters, turn on strong frustration
   10) Use $\mathcal{R}_0$ as the initial condition but evolve it using the new Hamiltonian's Superoperators

   Hopefully this will look like coherent, deterministic oscillations, like a fluid or waves.
   - **Arrows**: The classical expectation values of the spins will precess. Magnons should propagate across the grid.
   - **Bonds**: It'll oscillate. They should pulse and change color as the phase of the superconducting spin-liquid order parameters winds in time. It should display a time dependent breakdown and attempted recovery of the $U(1)$ symmetry.
 