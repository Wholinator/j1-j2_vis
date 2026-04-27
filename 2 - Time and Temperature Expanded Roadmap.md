# Post-Implementation Roadmap: HFB → TDHFB → FTHFB and Visualization

This roadmap picks up after the static-HFB code from `hfb_j1j2_plan.md` is
implemented, validated, and producing trustworthy 16×16 J1-J2 ground states. It
covers (1) hardening the static solver, (2) extending to time-dependent HFB
following Ebata et al., (3) extending to finite-temperature HFB following
Goodman, (4) building the Python visualization layer, and (5) the long-term
science program. Each phase has explicit entry criteria, deliverables, exit
criteria, and risks.

The reading anchors are:

- **Goodman (1981)**, "Finite-Temperature HFB Theory," Nucl. Phys. A352, 30.
  Defines the FTHFB equations, proves they have the same form as $T=0$ HFB
  with only the densities reweighted by Fermi-Dirac occupations, derives the
  grand-potential variational principle, and shows the key fact that
  $R^2 \ne R$ for $T>0$ removes the canonical-basis arbitrariness.
- **Ebata et al. (2010)**, "Canonical-basis time-dependent Hartree-Fock-Bogoliubov
  theory and linear-response calculations," Phys. Rev. C 82, 034306. Gives the
  TDHFB equation $i\dot R = [\mathcal{H}, R]$, its component form for $\rho$
  and $\kappa$, the conservation laws, and the linear-response recipe (kick
  the ground state, propagate, Fourier-transform).

---

## Phase 0: Static-HFB hardening (entry phase)

**Entry criteria.** All ten validation tests in §12 of `hfb_j1j2_plan.md`
pass. The two-site exact-spectrum test passes to $10^{-12}$. A 16×16 run from
at least three independent seeds converges to the same energy to $10^{-6}$ per
site for at least one $J_2/J_1$ value.

**Goals.** Make the static code reliable enough to be the foundation for
TDHFB and FTHFB. Bugs that survive into Phase 2 will be very hard to find
because TDHFB unitary evolution preserves errors in the initial state.

### 0.1 Persistent state and checkpoint format

Define a single binary serialization format used by the static solver, the
TDHFB propagator, and the FTHFB solver. HDF5 is the right choice: it handles
complex matrices natively, is readable from Python via `h5py` without extra
dependencies, and supports compression and partial reads.

```
hfb_state.h5
├── /metadata
│   ├── lattice_Lx (int)
│   ├── lattice_Ly (int)
│   ├── M (int)                       # 2 * Lx * Ly
│   ├── J1 (float), J2 (float)
│   ├── J_convention (string)         # "OrderedMatrixLiteral" | "SymmetricPhysicalBonds"
│   ├── beta (float)                  # 0 = T=0
│   ├── mu (float)
│   ├── N_target (float)
│   ├── git_commit (string)
│   ├── solver_version (string)
│   ├── timestamp_utc (string)
│   └── seed_label (string)           # which initialization seed produced this
├── /density
│   ├── rho (complex M×M)
│   └── kappa (complex M×M)
├── /generalized_density
│   └── R (complex 2M×2M)             # may be derived; store for convenience
├── /bdg
│   ├── eigenvalues (real 2M)
│   ├── U (complex M×2M)              # top half of eigenvector matrix W
│   └── V (complex M×2M)              # bottom half of W
├── /occupations
│   └── f (real 2M)                   # Fermi-Dirac occupations (T=0: 0/1)
├── /fields
│   ├── h_eff (complex M×M)           # h0 + Gamma - mu*I
│   ├── Gamma (complex M×M)
│   └── Delta (complex M×M)
├── /observables
│   ├── energy (float)
│   ├── entropy (float)               # 0 at T=0
│   ├── grand_potential (float)
│   ├── particle_number (float)
│   ├── site_density (real L)         # <n_i> per site
│   ├── site_magnetization (real L)   # <S^z_i> per site
│   └── bond_pairing (real, structured)  # <c_{i up} c_{j dn}> on each bond
└── /scf_history
    ├── energy_per_iter (real)
    ├── residual_per_iter (real)
    ├── mu_per_iter (real)
    └── number_per_iter (real)
```

C++ side: a thin wrapper around the HDF5 C++ API. Implement `save_state` and
`load_state` so a propagation can resume from any saved checkpoint.

### 0.2 Order parameter machinery

Even before TDHFB, the static code must compute and report:

- **Singlet pair amplitude on each bond.** For each bond $(i,j)$, the $S_z=0$
  singlet amplitude is $\Delta^{\text{singlet}}_{ij} = \tfrac{1}{\sqrt{2}} (\langle c_{i\uparrow} c_{j\downarrow}\rangle - \langle c_{i\downarrow} c_{j\uparrow}\rangle) = \tfrac{1}{\sqrt{2}} (\kappa_{j\downarrow, i\uparrow} - \kappa_{j\uparrow, i\downarrow})$.
- **$d$-wave vs $s$-wave decomposition** on NN bonds: project the NN bond
  amplitudes onto $\cos k_x \pm \cos k_y$ form factors.
- **Site magnetization** $\langle S^z_i \rangle = \tfrac{1}{2}(\rho_{i\uparrow,i\uparrow} - \rho_{i\downarrow,i\downarrow})$.
- **Static spin structure factor** $S(\mathbf{q}) = \tfrac{1}{L}\sum_{ij} e^{i\mathbf{q}\cdot(\mathbf{r}_i - \mathbf{r}_j)} \langle S^z_i S^z_j\rangle$, evaluated using Wick decomposition since we have $\rho, \kappa$.

These will be the variables we watch oscillate in TDHFB and decay with
temperature in FTHFB. Build them now while the static code is the only thing
that can be wrong.

### 0.3 Multi-seed driver and phase-diagram scan

Wrap the SCF in a driver that:

1. Sweeps $J_2/J_1 \in [0, 1]$ in steps of 0.05.
2. At each point runs all five seed types from §11.3 of the static plan.
3. Records converged energy, $\mu$, order parameters, and convergence status.
4. Produces an HDF5 file per converged solution and a summary CSV/JSON of the
   ensemble.

This scan is the science-grade product of Phase 0 and the input to all later
phases. Time-evolution tests start from one of these checkpoints.

**Exit criteria for Phase 0.**

- Phase scan complete across $J_2/J_1 \in [0, 1]$ at 16×16.
- For each $J_2/J_1$ point: the lowest-energy seed agrees with at least one other seed to $10^{-6}$ per site, or the energy gap between competing minima is documented.
- Round-trip checkpoint test: load a saved state, recompute energy/fields, agree with stored values to $10^{-10}$.
- Energy cross-check (Wick vs trace) passes for every saved state.

---

## Phase 1: Static-HFB extensions

These are static-HFB features that are useful on their own and are
prerequisites for the dynamical extensions. Pick them off in any order.

### 1.1 Local single-occupancy Lagrange multipliers

Even though we're not enforcing local single occupancy as a hard constraint,
adding the *option* to penalize deviations is cheap and useful for comparing
to constrained mean-field treatments. Add:

$$h_{i\sigma, i\sigma} \mapsto h_{i\sigma, i\sigma} + \lambda_i$$

with the soft update rule

$$\lambda_i \leftarrow \lambda_i + \eta_\lambda (\langle n_{i\uparrow} + n_{i\downarrow}\rangle - 1)$$

on every SCF step. Default $\eta_\lambda = 0$ (off). When on, monitor
$\max_i |\langle n_i\rangle - 1|$ as a convergence criterion in addition to
the existing ones.

### 1.2 $S_z$-block-reduced BdG solver

When the initialization preserves $S_z$ (always, in our default workflow), the
$1024 \times 1024$ BdG matrix reduces to a $1024 \times 1024$ Hermitian matrix
that decouples into two equivalent $1024 \times 1024$ blocks via the Nambu
mixing of $(c_{i\uparrow}, c_{j\downarrow}^\dagger)$. The reduction gives
roughly a factor of 4 in diagonalization cost and a factor of 2 in memory.

Implementation: a separate `SzBlockBdgSolver` class that takes
$h_{\uparrow\uparrow}$, $h_{\downarrow\downarrow}$, and
$\Delta_{\uparrow\downarrow}$ as inputs and returns half the spectrum plus
the corresponding eigenvectors. The full-BdG path remains available as a
ground truth for cross-validation.

**Test.** For a converged state from Phase 0, run both solvers and check the
extracted $\rho, \kappa$ agree to $10^{-10}$.

### 1.3 Anderson mixer history persistence

To allow restart-and-continue runs to avoid re-warming the mixer, save the
Anderson history vectors to the checkpoint (as `/mixer/X_history`,
`/mixer/F_history`).

### 1.4 Diagonalization performance

Switch the production path from `Eigen::SelfAdjointEigenSolver` to LAPACK
`zheevd`. Benchmark on the 16×16 BdG matrix; expect 3-10× speedup with MKL.
Add a `--diagonalizer` CLI flag for A/B comparison during development.

**Exit criteria for Phase 1.** All extensions implemented, each gated by a
test that compares against the Phase-0 baseline. No regressions on the
Phase-0 scan.

---

## Phase 2: Time-Dependent HFB

This phase implements full TDHFB following Ebata et al. §II.B. We explicitly
**do not** adopt the canonical-basis approximation (Eq. 28 of the paper) —
that approximation forces $\Delta$ to be diagonal in the canonical basis,
which for our problem means BCS-like singlet pairing only and would lose
much of the physics we want to study. We propagate the full $\rho, \kappa$.

### 2.1 The equations to solve

In the same conventions as the static plan ($\rho_{pq} = \langle c_q^\dagger c_p\rangle$,
$\kappa_{pq} = \langle c_q c_p\rangle$, $h$ Hermitian, $\Delta^T = -\Delta$):

$$i \frac{\partial \rho}{\partial t} = [h, \rho] + \kappa \Delta^* - \Delta \kappa^*$$

$$i \frac{\partial \kappa}{\partial t} = h \kappa + \kappa h^* + \Delta (1 - \rho^*) - \rho \Delta$$

These are Eqs. (8) and (9) of Ebata et al. Equivalently, the generalized
density $R$ obeys $i \partial_t R = [\mathcal{H}, R]$, which is the form we
will integrate.

The Hamiltonian $\mathcal{H}(t) = \mathcal{H}[R(t)]$ depends on $R$ at the
*same* time, so this is a nonlinear ODE. Integrators must be chosen accordingly.

### 2.2 Integrator selection

Simple forward Euler (which Ebata et al. used for their Cb form, Eq. 72-74)
is unstable for the full TDHFB matrix equation at the time steps we want.
Three options, in increasing order of cost and quality:

**(a) Crank-Nicolson with a fixed-point inner loop.** Implicit, unitary in
the linear sense, second-order accurate. The inner loop converges in 2-3
iterations for small steps. Good default.

$$R^{n+1} = R^n - i \, \mathrm{d}t \cdot [\bar{\mathcal{H}}, \bar R], \quad \bar X = \tfrac{1}{2}(X^n + X^{n+1})$$

**(b) Magnus expansion (4th order).** $R^{n+1} = e^{-i\Omega} R^n e^{i\Omega}$
with $\Omega$ a Magnus combination of $\mathcal{H}$ at intermediate times.
Exactly preserves $R^\dagger = R$, $R$ block structure, and (for time-independent
$\mathcal{H}$) energy. Best for long-time runs.

**(c) Predictor-corrector with adaptive step.** Lighter weight than Magnus,
heavier than Crank-Nicolson. Useful if $\mathcal{H}(t)$ has stiff features.

**Recommendation:** implement (a) first as the workhorse, then (b) for
long-time linear-response runs. The Cb-TDHFB paper's simple Euler (Eq. 72)
is not adequate for full TDHFB at the matrix level.

### 2.3 Conservation diagnostics

At every step, monitor:

- **Total energy** $E(t)$: should be conserved to $10^{-8}$ (relative) for
  Crank-Nicolson with reasonable step sizes, $10^{-12}$ for Magnus 4.
- **Particle number** $N(t) = \mathrm{Tr}\,\rho(t)$: conserved to roundoff.
- **Idempotency at $T=0$**: $\|R(t)^2 - R(t)\|_F$ should remain at $10^{-10}$
  level. Drift indicates integrator error.
- **Hermiticity** $\|\rho - \rho^\dagger\|$ and antisymmetry
  $\|\kappa + \kappa^T\|$.

If any of these drift, the integrator step is too big or the symmetry
projection between steps is needed. Symmetry projection (the
$\rho \leftarrow \tfrac{1}{2}(\rho + \rho^\dagger)$ trick) can restore
constraints but masks integrator quality, so use it sparingly during
development and report drift before projection.

### 2.4 Linear response: the kick-and-propagate recipe

Following Ebata et al. §V.D:

1. Load a converged static HFB state.
2. Apply a one-body instantaneous kick:
   $\rho \to e^{i\eta F} \rho \, e^{-i\eta F}$,
   $\kappa \to e^{i\eta F} \kappa \, e^{i\eta \tilde F}$
   for a chosen Hermitian one-body operator $F$ and small $\eta$ (linear regime).
3. Propagate from $t=0$ to $t = T_{\max}$, recording $\langle F(t)\rangle$ or
   another observable of interest.
4. Compute the response function via Fourier transform with Gaussian or
   exponential damping (Eq. 76 of Ebata et al.):

$$\chi_{FF}(\omega) = -\frac{1}{\pi \eta} \, \mathrm{Im} \int_0^\infty dt \, e^{i\omega t - \Gamma t / 2} [\langle F\rangle(t) - \langle F\rangle(0)]$$

For J1-J2 the natural choices for $F$ are:

- **Magnetic-field probes**: $F_z = \sum_i e^{i\mathbf{q}\cdot\mathbf{r}_i} S^z_i$ for the dynamical spin structure factor at wavevector $\mathbf{q}$.
- **Pairing probes**: $F = \sum_{ij} f_{ij} (c_{i\uparrow} c_{j\downarrow} + \mathrm{h.c.})$ to drive Higgs-like amplitude oscillations.
- **Bond probes**: $F = \sum_{ij} g_{ij} \mathbf{S}_i \cdot \mathbf{S}_j$ to drive density-of-states-like responses.

The output spectrum gives the mean-field RPA approximation to the
corresponding dynamical susceptibility.

### 2.5 Validation tests for TDHFB

1. **Trivial dynamics from a stationary state.** Start from a converged HFB
   solution with no perturbation. Propagate for $T = 100/J_1$. All observables
   should be constant to $10^{-8}$. This catches integrator bugs and any
   residual non-stationarity in the static solution.

2. **Free-fermion limit.** Set $V = 0$, add a small NN hopping $t_h$ in $h^0$.
   Kick with a one-body density operator. The propagation reduces to TDHF on
   a band; compare the response to direct band-structure calculation.

3. **Energy conservation.** With Crank-Nicolson at $\mathrm{d}t = 0.01/J_1$,
   energy should be conserved to $10^{-6}$ over $T = 50/J_1$. Magnus 4 should
   reach $10^{-10}$.

4. **Number conservation.** $|N(t) - N(0)| < 10^{-10}$ throughout.

5. **Idempotency at $T=0$.** $\|R(t)^2 - R(t)\|_F < 10^{-9}$ over the run.

6. **Goldstone-mode test (Ebata §III.D).** For a paired ground state, the
   pairing-rotation transformation $R \to e^{i\theta N} R \, e^{-i\theta N}$
   should be a zero-energy mode: applying an infinitesimal kick of this form
   should produce dynamics at $\omega = 0$. The peak in the response at zero
   frequency confirms the Nambu-Goldstone mode is correctly reproduced.

7. **Small-amplitude → linear response sanity.** Reduce the kick amplitude
   $\eta$ by 10× and check the response spectrum scales as $\eta^2$ (intensity)
   while peak positions are unchanged. Departure from linearity sets the
   maximum usable $\eta$.

8. **Time-reversal test.** Propagate forward to $T$, then backward to $0$;
   $R(2T)$ should match $R(0)$ to integrator precision.

### 2.6 Performance and parallelism

Per Crank-Nicolson step: 2-3 BdG-style matrix exponentials (or equivalent
linear solves) plus field rebuilds. Costs roughly $10 \times$ a static SCF
step. A typical linear-response run with $T = 50/J_1$ and
$\mathrm{d}t = 0.02/J_1$ is 2500 steps, so $\sim 25{,}000$ static-equivalent
operations per run. Tractable on a workstation overnight.

For Magnus 4, each step is roughly $4\times$ Crank-Nicolson. Use Magnus only
when conservation matters more than wall time.

**Parallelism.** Different perturbations $F$ and different ground states
(different $J_2/J_1$) are embarrassingly parallel. Use a process pool at the
driver level. Within a single propagation, MKL/OpenBLAS already parallelizes
the dense linear algebra.

**Exit criteria for Phase 2.**

- All eight TDHFB validation tests pass.
- Linear-response spectrum from a known-paired Phase-0 ground state has its
  expected Nambu-Goldstone peak at $\omega = 0$.
- Magnus-4 energy drift below $10^{-9}$ over $T = 100/J_1$.
- Documented spectrum file format (HDF5 with $\langle F(t)\rangle$, $\chi(\omega)$,
  metadata identifying the source ground state and kick operator).

---

## Phase 3: Finite-Temperature HFB

Following Goodman §3-4. Goodman's central result is that **FTHFB has the same
matrix form as $T=0$ HFB**. The fields $\Gamma, \Delta$ are built with the
same formulas, the BdG eigenvalue problem is the same, and only the
construction of $\rho, \kappa$ from the eigenvectors changes.

### 3.1 The equations to solve

From Goodman Eqs. (3.19), (3.20), (3.27), (4.20):

$$\rho = U f \tilde U^* + V^\dagger (1-f) V$$

$$\kappa = U f \tilde V^* + V^\dagger (1-f) U$$

$$f_i = \frac{1}{1 + e^{\beta E_i}}, \qquad S = -k_B \sum_i [f_i \ln f_i + (1-f_i)\ln(1-f_i)]$$

where $E_i$ are the (positive) BdG quasiparticle energies and $(U, V)$ are
their eigenvectors. The grand potential is $\Omega = E - TS - \mu N$, and
the FTHFB self-consistency requires *both* the orbitals and the occupations
to converge (Goodman §4 "second condition is new").

### 3.2 Implementation: minimal delta from $T=0$

The static code already implements `density_finiteT(H, beta)` as an
alternative to `density_zeroT(H)`. The remaining changes are:

1. **Entropy and grand-potential reporting.** Compute and store $S$ and
   $\Omega$ alongside $E$. Convergence is now on $\Omega$, not $E$, because
   $\Omega$ is what's being minimized. Energy is no longer monotone during
   FTHFB SCF; $\Omega$ is.

2. **Pairing convergence diagnostic.** At finite $T$ the order parameter can
   genuinely vanish at a critical $\beta_c$. Plot $\|\Delta\|$ vs $\beta$
   from a $\beta$ sweep and identify $\beta_c$ as a smooth-but-rapid drop
   (the mean-field "phase transition," Goodman Fig. 1).

3. **Occupation convergence.** Track $\max_i |f_i^{(n+1)} - f_i^{(n)}|$ as
   an additional convergence criterion. This is automatic if the residual
   on $R$ is converged, but reporting it explicitly helps diagnose stalled
   SCF runs near $\beta_c$.

4. **Annealing schedule.** For large lattices and small gaps, start at high
   $T$ (e.g., $\beta J_1 = 1$), converge, lower $T$ by a factor of $\sim 1.3$,
   converge again using the previous solution as the seed, until reaching
   the target $\beta$. This is much more reliable than starting at low $T$.

5. **$T = 0$ limit.** At very large $\beta$ (e.g., $\beta J_1 > 10^4$),
   the smooth $f$ becomes numerically degenerate with the $T=0$ projector.
   Detect this and either switch to the projector path or increase the
   step-function threshold in `density_finiteT`. The validation suite must
   include: at $\beta = 10^6$, FTHFB and $T=0$ HFB agree to $10^{-9}$.

### 3.3 Validation tests for FTHFB

1. **High-$\beta$ recovery of $T=0$.** Run FTHFB at $\beta J_1 = 10^4$ from
   a $T=0$ checkpoint. The converged state should match the $T=0$ state to
   $10^{-7}$ in $\rho, \kappa$ and to $10^{-9}$ in energy.

2. **Entropy bounds.** $S \ge 0$ at every iteration. $S \le M \ln 2$ (the
   maximum-entropy bound for $M$ two-level systems).

3. **Grand-potential monotonicity.** With small enough mixing $\Omega$
   should be non-increasing across SCF steps once linear-mixing warmup is
   done. Use this to tune the mixing parameters and the annealing schedule.

4. **Goodman degenerate-shell sanity check.** For the analytically solvable
   half-filled degenerate-shell pairing model (Goodman §6), the code should
   reproduce $\Delta(T)$ from his Eq. (6.2) to four digits and the critical
   temperature $k_B T_c = \Delta_0 / 2$ to better than 1%. This is a
   non-trivial test of the SCF logic at finite $T$.

5. **Specific-heat anomaly.** Compute $C_v = \beta^2 \partial^2 (\beta \Omega) / \partial \beta^2$
   numerically by central differences across a $\beta$ sweep. Around
   the pairing transition, $C_v$ should show the mean-field jump
   characteristic of a second-order transition (or the smooth crossover for
   a finite-system calculation). The presence and shape of this feature is
   a sensitive test of the FTHFB convergence quality.

6. **Order-parameter melting curve.** Sweep $\beta$ and plot the singlet pair
   amplitude vs $T$. Compare qualitatively to Goodman Fig. 1 and document
   $T_c$ as a function of $J_2/J_1$.

### 3.4 Time-dependent finite-temperature HFB (stretch goal)

Once both Phase 2 and Phase 3 are working, combining them gives TDHFB at
finite $T$: propagate $R(t)$ as in Phase 2, but starting from a finite-$T$
self-consistent state (where $R^2 \ne R$). The same equations apply; the
diagnostics differ (no idempotency; entropy is conserved by the TDHFB
equations). The dynamical response at finite $T$ gives access to thermal
broadening and order-parameter relaxation.

**Exit criteria for Phase 3.**

- Six FTHFB validation tests pass.
- $T$-sweep across $J_2/J_1 \in [0.0, 0.5, 1.0]$ produces $\Delta(T)$ curves
  with identifiable mean-field $T_c$.
- High-$\beta$ recovery of $T=0$ ground state passes for at least three
  $J_2/J_1$ values.

---

## Phase 4: Python visualization and analysis layer

The C++ side handles SCF, time evolution, and observable computation. The
Python side handles everything users actually look at: plots, diagnostics,
exploratory analysis, and report generation. The contract between the two
is the HDF5 schema in §0.1.

### 4.1 Repository layout

```
hfb_j1j2/
├── ... (existing C++ source)
└── pyhfb/
    ├── pyhfb/
    │   ├── __init__.py
    │   ├── io.py              # HDF5 readers/writers; thin wrappers
    │   ├── observables.py     # derived quantities: S(q), correlations, etc.
    │   ├── plots/
    │   │   ├── lattice.py     # site-resolved plots: magnetization, density
    │   │   ├── bond.py        # bond-resolved: pairing, current
    │   │   ├── spectrum.py    # BdG spectrum, response spectra
    │   │   ├── phase.py       # phase-diagram and T-sweep plots
    │   │   └── dynamics.py    # animations of TDHFB observables
    │   ├── tests/
    │   └── cli.py             # typer-based command-line driver
    ├── pyproject.toml
    └── README.md
```

Build with `pip install -e ./pyhfb`. Pin to Python 3.11+ for typing
ergonomics; depend on `numpy`, `scipy`, `h5py`, `matplotlib`, `pandas`,
`typer`. Use `seaborn` for statistical plots if needed.

### 4.2 Core data class

```python
@dataclass
class HFBState:
    metadata: dict
    rho: np.ndarray            # complex (M, M)
    kappa: np.ndarray          # complex (M, M)
    R: np.ndarray              # complex (2M, 2M), Hermitian
    eigenvalues: np.ndarray    # real (2M,)
    U: np.ndarray              # complex (M, 2M)
    V: np.ndarray              # complex (M, 2M)
    f: np.ndarray              # real (2M,)
    h_eff: np.ndarray
    Gamma: np.ndarray
    Delta: np.ndarray
    energy: float
    entropy: float
    grand_potential: float
    particle_number: float
    site_density: np.ndarray
    site_magnetization: np.ndarray
    bond_pairing: np.ndarray   # structured: bond index -> amplitude

    @classmethod
    def load(cls, path: pathlib.Path) -> "HFBState": ...

    @property
    def L(self) -> int:
        return self.metadata["lattice_Lx"] * self.metadata["lattice_Ly"]

    def site_index(self, x: int, y: int) -> int: ...
    def spin_orbital(self, site: int, spin: int) -> int: ...
```

Also define `HFBTrajectory` for TDHFB output (sequence of states or, more
efficiently, time series of observables) and `HFBTSweep` for finite-$T$
sweeps.

### 4.3 Required visualizations

**Static / single-state plots.**

1. **Lattice heatmap of site magnetization.** $\langle S^z_i \rangle$ on the
   $L_x \times L_y$ grid, diverging colormap centered at zero. Identifies
   Néel vs stripe vs disordered states at a glance.

2. **Bond-resolved pairing amplitude.** Each NN and NNN bond drawn with line
   thickness proportional to $|\Delta^{\text{singlet}}_{ij}|$ and color
   indicating sign/phase. The $d$-wave vs $s$-wave character is visually
   obvious in this view.

3. **BdG quasiparticle spectrum.** Eigenvalues plotted with horizontal
   ticks; mark the gap, the chemical potential, and the highest occupied
   state. Useful for spotting near-zero modes that might cause SCF
   instability.

4. **Static spin structure factor $S(\mathbf{q})$.** Heatmap on the
   Brillouin zone. The position of the peak distinguishes phases:
   $(\pi, \pi)$ for Néel, $(\pi, 0)$ for stripe.

5. **Density-of-states (DOS).** Smoothed histogram of BdG eigenvalues. The
   pairing gap appears as a hard minimum at $\omega = 0$.

**Phase-diagram plots.**

6. **Energy vs $J_2/J_1$ for each seed.** Lowest-energy envelope is the
   ground state; crossings indicate mean-field phase transitions.

7. **Order parameters vs $J_2/J_1$.** Néel order parameter, stripe order
   parameter, $d$-wave amplitude, and $s$-wave amplitude on the same axes.
   The competition between phases shows up as one rising while another
   falls.

**TDHFB plots.**

8. **Real-time observables.** $\langle F(t)\rangle$ for the chosen probe.
   Useful diagnostic plot during development.

9. **Response spectrum $\chi(\omega)$.** Imaginary part vs $\omega$ with
   smoothing. Mark expected Goldstone mode at $\omega = 0$ and any gap
   features.

10. **Dynamical structure factor $S(\mathbf{q}, \omega)$.** Heatmap of
    response strength on the $(\mathbf{q}, \omega)$ plane along
    high-symmetry lines of the Brillouin zone. This is the headline
    deliverable for the TDHFB phase: directly comparable to inelastic
    neutron scattering experiments.

**FTHFB plots.**

11. **Order parameter vs $T$.** Singlet amplitude, magnetization, etc.,
    versus $k_B T / J_1$. Mark $T_c$. Compare to Goodman Fig. 1 for
    qualitative shape.

12. **Free-energy and specific-heat curves.** $\Omega(T)$ and $C_v(T)$ vs
    $T$. Anomaly at $T_c$ should be visible.

13. **Phase diagram in $(J_2/J_1, T)$ plane.** Heatmap of order parameter
    or contour plot of $T_c$ vs $J_2/J_1$.

### 4.4 CLI interface

```
$ pyhfb summarize state.h5
$ pyhfb plot site-magnetization state.h5 --output mag.png
$ pyhfb plot pairing-bonds state.h5 --output pairing.png --component dwave
$ pyhfb plot spectrum state.h5 --output spectrum.png
$ pyhfb plot phase-scan results/scan_*.h5 --output phase.png
$ pyhfb plot response trajectory.h5 --probe Sz_pi_pi --output response.png
$ pyhfb plot dynamical-sf trajectory_*.h5 --output sqw.png
$ pyhfb plot order-vs-temperature tsweep_*.h5 --output melting.png
$ pyhfb animate magnetization trajectory.h5 --output mag.mp4 --fps 30
```

### 4.5 Validation of the visualization layer

1. **Round-trip test.** Save a state, load it via `pyhfb.io`, recompute
   energy and observables, compare to stored values: $10^{-12}$ agreement.

2. **Schema test.** Reject HDF5 files missing required fields with clear
   error messages.

3. **Plot regression test.** Use `pytest-mpl` to compare generated plots
   against committed PNG baselines. Fail on $> 1\%$ pixel difference.
   Catches matplotlib backend changes and numerical regressions in
   observables.

**Exit criteria for Phase 4.** All 13 plot types produce reasonable output
on the Phase-0 scan, the Phase-2 trajectories, and the Phase-3 sweeps.
Round-trip, schema, and plot regression tests all pass.

---

## Phase 5: Long-term science program

With all infrastructure in place, the actual physics:

### 5.1 Mean-field phase diagram of J1-J2 at $T=0$

Sweep $J_2/J_1 \in [0, 1]$ at 16×16 with at least 5 random seeds plus the
four structured seeds (uniform, Néel, stripe, $d$-wave). Identify:

- The Néel-to-disordered transition near $J_2/J_1 \approx 0.4$.
- The disordered-to-stripe transition near $J_2/J_1 \approx 0.6$.
- Any paired/spin-liquid mean-field minima in the intermediate region.

Compare to QMC and DMRG results in the literature. HFB will give an upper
bound on the ground-state energy; the gap is the mean-field error.

### 5.2 Finite-size scaling

Repeat the scan at 8×8, 12×12, 16×16, 20×20. Watch how the phase boundaries
and the order-parameter magnitudes scale. Mean-field correlations are
infinite-range; the finite-size scaling will be qualitatively different
from QMC.

### 5.3 Dynamical signatures

For each $J_2/J_1$ in the converged scan, run TDHFB with kicks at
characteristic wavevectors:

- $\mathbf{q} = (\pi, \pi)$: Néel response.
- $\mathbf{q} = (\pi, 0)$: stripe response.
- $\mathbf{q} = (0, 0)$: uniform mode, includes Goldstone.
- $\mathbf{q} = (\pi/2, \pi/2)$: probe of intermediate physics.

Plot $S(\mathbf{q}, \omega)$ along high-symmetry paths. The presence of
Goldstone modes, gap structure, and continuum features is the science
output.

### 5.4 Finite-$T$ phase diagram

Sweep $\beta$ at the same $J_2/J_1$ points. Locate the mean-field
transitions in the $(J_2/J_1, T)$ plane. Compute $T_c$ as a function of
$J_2/J_1$ and compare to expectations from the $T=0$ order-parameter
magnitudes (mean field predicts $T_c \propto |\Delta_0|$).

### 5.5 Systematic improvements (future work)

Beyond the scope of this plan but worth flagging:

- **Number projection.** HFB breaks particle-number symmetry; projection
  restores it and lowers the energy. Implement via the Lipkin-Nogami
  approximation as a first step, exact projection later.
- **Spin projection.** Even if $S_z$ is preserved, total $S^2$ is not.
  Spin projection lowers energy further.
- **Beyond mean field.** Random phase approximation built on HFB ground
  state, then second RPA, then symmetry-projected HFB. Each level removes
  more of the mean-field error.
- **Different lattices.** Triangular, kagome, honeycomb. The infrastructure
  generalizes; only the lattice and bond definitions change.

---

## Master schedule and dependencies

```
Phase 0 (hardening)                  
  └── HDF5 schema, observables, multi-seed scanner

Phase 1 (extensions)                
  ├── 1.1 local Lagrange multipliers (optional)
  ├── 1.2 Sz block solver
  ├── 1.3 mixer persistence
  └── 1.4 LAPACK switch

Phase 2 (TDHFB)                     
  ├── 2.1-2.3 integrators and diagnostics
  ├── 2.4 linear response
  └── 2.5 validation

Phase 3 (FTHFB)                     
  ├── 3.1-3.2 entropy, grand potential, annealing
  └── 3.3 validation

Phase 4 (Python viz)                
  ├── 4.1-4.2 io and core class
  ├── 4.3 plots
  └── 4.4-4.5 CLI and tests

Phase 5 (science)                   
```

Total infrastructure time: 8-12 weeks. Phase 4 can start in parallel once
the HDF5 schema is locked, even before Phase 2 is complete — early plots
of Phase-0 ground states will catch issues with the static solver before
they propagate.

---

## Common pitfalls and how to avoid them

1. **Storing $R$ but not $U, V$.** $R$ is enough for static observables
   but loses the eigenstructure needed to restart TDHFB cleanly or to
   compute quasiparticle-resolved quantities. Always save eigenvectors.

2. **Forgetting to save metadata.** A folder of HDF5 files with no record
   of which $J_2/J_1$ or $\beta$ each one represents is useless within
   weeks. The metadata block in §0.1 is mandatory, including a git commit
   hash so you can reproduce the run.

3. **TDHFB integrator drift masquerading as physics.** If energy drifts
   monotonically over a long run, it's the integrator, not the physics.
   Always plot $E(t)$, not just the final value. Magnus 4 should give
   essentially flat energy; if it doesn't, the field-build code has a bug.

4. **FTHFB SCF stalling near $T_c$.** The order parameter is genuinely
   trying to vanish, and the SCF iteration is sensitive to mixing. Annealing
   from high $T$ down through $T_c$ is much more reliable than direct
   convergence at a single $T$ near the transition.

5. **Confusing visualization with validation.** A pretty heatmap is not a
   correctness check. Plot regression tests (`pytest-mpl`) catch subtle
   numerical drift; visual inspection alone does not.

6. **Coupling Python to C++ too tightly.** Resist the urge to call the
   C++ solver from Python via pybind11 in Phase 4. The HDF5 boundary is
   slower per-call but enormously more robust: you can run C++ on a
   cluster, copy files locally, and analyze with Python interactively.
   Coupling them tightly creates a build-system mess and a hard dependency
   on identical compiler/library versions.

7. **Not budgeting for re-runs.** Bugs found in Phase 2 may invalidate
   Phase 0 outputs. Keep the Phase-0 scanner runnable from a single command
   so re-runs are cheap.
