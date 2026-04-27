# Abrikosov-Fermion HFB for the J1–J2 Model: Final Implementation Plan

## 0. Scope and design goals

This plan specifies a C++ implementation of unrestricted Hartree-Fock-Bogoliubov (HFB) mean-field theory for the spin-1/2 J1-J2 Heisenberg model on a 16×16 square lattice, using the Abrikosov (parton) fermion representation. Number symmetry is allowed to break; local single-occupancy is **not** enforced (only the global parton number). The plan is the consolidation of two design reviews and supersedes both.

**Targets.**
- Lattice: 16×16, periodic boundaries (configurable).
- Spin orbitals: $M = 2L = 512$. Full BdG matrix: $2M \times 2M = 1024 \times 1024$.
- Diagonalization: dense Hermitian via LAPACK (`zheevd` or `zheevr`); Eigen for prototyping.
- Field construction: sparse antisymmetrized two-body tensor with direct loop.
- Mixing: Anderson/Pulay from the start (with linear-mixing warmup).
- Validation: mandatory two-site exact-spectrum test plus a full diagnostic suite.

---

## 1. Conventions (fixed once, used everywhere)

### 1.1 Indexing

Spatial sites: $i = x + L_x y$ with $0 \le x < L_x$, $0 \le y < L_y$. Spin: $\sigma \in \{0,1\} = \{\uparrow, \downarrow\}$. Spin orbital index (spin-major):

$$p(i,\sigma) = i + \sigma L$$

```cpp
enum Spin : int { UP = 0, DOWN = 1 };

inline int so(int site, Spin spin, int L) {
    return site + static_cast<int>(spin) * L;
}
```

### 1.2 Densities

$$\rho_{pq} \equiv \langle c_q^\dagger c_p \rangle, \qquad \kappa_{pq} \equiv \langle c_q c_p \rangle$$

These satisfy

$$\rho^\dagger = \rho, \qquad \kappa^T = -\kappa$$

### 1.3 Generalized density and BdG matrix

Nambu spinor $\Psi = (c_1, \ldots, c_M, c_1^\dagger, \ldots, c_M^\dagger)^T$. Generalized density:

$$R = \begin{pmatrix} \rho & \kappa \\ -\kappa^* & I - \rho^* \end{pmatrix}$$

BdG/quasiparticle Hamiltonian:

$$\mathcal{H} = \begin{pmatrix} h - \mu I & \Delta \\ -\Delta^* & -(h - \mu I)^* \end{pmatrix}, \quad h = h^\dagger, \quad \Delta^T = -\Delta$$

### 1.4 Two-body operator and antisymmetrized tensor

Every term is stored as a normal-ordered monomial

$$C \, c_p^\dagger c_q^\dagger c_s c_r$$

(operator order: $p, q, s, r$). The compressed antisymmetrized two-body Hamiltonian is

$$H_2 = \frac{1}{4} \sum_{pqrs} V_{pqrs} \, c_p^\dagger c_q^\dagger c_s c_r$$

with

$$V_{pqrs} = -V_{qprs} = -V_{pqsr} = V_{qpsr}$$

### 1.5 Mean fields

$$\Gamma_{pq}[\rho] = \sum_{rs} V_{prqs} \, \rho_{sr}$$

$$\Delta_{pq}[\kappa] = \frac{1}{2} \sum_{rs} V_{pqrs} \, \kappa_{rs}$$

These follow from $\delta E / \delta \rho_{qp}$ and $\delta E / \delta \kappa_{qp}^*$ with the conventions above.

### 1.6 Wick contractions (canonical, verified)

Under the stated conventions:

$$\boxed{\langle c_p^\dagger c_q^\dagger c_s c_r \rangle = \rho_{rp} \rho_{sq} - \rho_{sp} \rho_{rq} + \kappa_{pq}^* \kappa_{rs}}$$

**Derivation (the sign that everything depends on).** Three Wick pairings:

1. $\langle c_p^\dagger c_r \rangle \langle c_q^\dagger c_s \rangle = \rho_{rp} \rho_{sq}$ (using $\langle c_b^\dagger c_a \rangle = \rho_{ab}$).
2. $-\langle c_p^\dagger c_s \rangle \langle c_q^\dagger c_r \rangle = -\rho_{sp} \rho_{rq}$.
3. $+\langle c_p^\dagger c_q^\dagger \rangle \langle c_s c_r \rangle$.

For the anomalous piece: $\langle c_s c_r \rangle$ matches $\kappa_{ab} = \langle c_b c_a \rangle$ with $b = s$, $a = r$, giving $\kappa_{rs}$. And $\langle c_p^\dagger c_q^\dagger \rangle = \overline{\langle c_q c_p \rangle}$; the inner expectation matches $\kappa_{ab} = \langle c_b c_a \rangle$ with $b = q$, $a = p$, giving $\kappa_{pq}$. So $\langle c_p^\dagger c_q^\dagger \rangle = \kappa_{pq}^*$, and the term is $+\kappa_{pq}^* \kappa_{rs}$.

This is the **plus sign** convention used throughout. It is consistent with $\Delta_{pq} = \frac{1}{2} \sum V_{pqrs} \kappa_{rs}$.

---

## 2. Hamiltonian input: bond-counting policy

### 2.1 The user's expression

$$H = \sum_{i,j} J_{ij} \left[ c_{i\uparrow}^\dagger c_{j\downarrow}^\dagger c_{j\uparrow} c_{i\downarrow} + \tfrac{1}{4} c_{i\uparrow}^\dagger c_{j\uparrow}^\dagger c_{j\uparrow} c_{i\uparrow} + \tfrac{1}{4} c_{i\downarrow}^\dagger c_{j\downarrow}^\dagger c_{j\downarrow} c_{i\downarrow} - \tfrac{1}{2} c_{i\uparrow}^\dagger c_{j\downarrow}^\dagger c_{j\downarrow} c_{i\uparrow} \right]$$

This is naturally an **ordered-pair** sum. Loop over all nonzero $(i,j)$ entries of the input matrix as written.

### 2.2 The factor-of-two trap

If the user wants the physical Heisenberg Hamiltonian

$$H = J_1 \sum_{\langle ij \rangle_1} \mathbf{S}_i \cdot \mathbf{S}_j + J_2 \sum_{\langle ij \rangle_2} \mathbf{S}_i \cdot \mathbf{S}_j$$

(unordered bonds, each bond counted once), and feeds in a symmetric matrix with $J_{ij} = J_{ji} = J_n$, the ordered double-sum will double-count each bond. The two valid input modes:

```cpp
enum class JConvention {
    OrderedMatrixLiteral,    // sum exactly as supplied
    SymmetricPhysicalBonds   // symmetric matrix; halve each entry internally
};
```

For `SymmetricPhysicalBonds`, set $J_{ij} = J_{ji} = J_n / 2$ on each bond pair. The internal solver always uses the ordered-pair monomial expansion. Document this prominently in the API:

```cpp
// Build a symmetric J matrix encoding the physical Heisenberg Hamiltonian
//   H = J1 * sum_{<ij>_1} S_i.S_j + J2 * sum_{<ij>_2} S_i.S_j
// using the ordered-pair convention. Each bond pair (i,j),(j,i) is set
// to J_n / 2 so that the ordered double-sum reproduces H exactly.
Eigen::MatrixXd build_J1J2_matrix(int Lx, int Ly,
                                   double J1, double J2,
                                   bool periodic = true);
```

### 2.3 Sanity assertions on the input

```cpp
assert(J.rows() == J.cols());
assert(J.rows() == L);
for (int i = 0; i < L; ++i) {
    assert(J(i,i) == 0.0);              // no self-coupling
}
// For SymmetricPhysicalBonds mode:
assert((J - J.transpose()).norm() < 1e-12);
```

---

## 3. Monomial representation

Each ordered pair $(i,j)$ with $J_{ij} \ne 0$ generates **four** monomials:

```cpp
struct Monomial {
    int p, q;                    // creators in order c_p^\dagger c_q^\dagger
    int r, s;                    // annihilators: operator is c_p^\dag c_q^\dag c_s c_r
    std::complex<double> coeff;
};

void add_J_pair_monomials(int i, int j, double Jij, int L,
                           std::vector<Monomial>& mons) {
    using cplx = std::complex<double>;

    auto add = [&](int p, int q, int r, int s, cplx C) {
        assert(p != q);
        assert(r != s);
        mons.push_back({p, q, r, s, C});
    };

    // c_{i up}^dag c_{j dn}^dag c_{j up} c_{i dn}        coeff = J
    add(so(i,UP,L), so(j,DOWN,L), so(i,DOWN,L), so(j,UP,L), Jij);

    // (1/4) c_{i up}^dag c_{j up}^dag c_{j up} c_{i up}  coeff = J/4
    add(so(i,UP,L), so(j,UP,L), so(i,UP,L), so(j,UP,L), 0.25 * Jij);

    // (1/4) c_{i dn}^dag c_{j dn}^dag c_{j dn} c_{i dn}  coeff = J/4
    add(so(i,DOWN,L), so(j,DOWN,L), so(i,DOWN,L), so(j,DOWN,L), 0.25 * Jij);

    // (-1/2) c_{i up}^dag c_{j dn}^dag c_{j dn} c_{i up} coeff = -J/2
    add(so(i,UP,L), so(j,DOWN,L), so(i,UP,L), so(j,DOWN,L), -0.5 * Jij);
}
```

**Index decoding** (operator $c_p^\dagger c_q^\dagger c_s c_r$):
- For the spin-flip term $c_{i\uparrow}^\dagger c_{j\downarrow}^\dagger c_{j\uparrow} c_{i\downarrow}$: $p = i\!\uparrow$, $q = j\!\downarrow$, then the rightmost annihilator is $c_{i\downarrow}$ so $r = i\!\downarrow$, and the inner annihilator is $c_{j\uparrow}$ so $s = j\!\uparrow$.

Loop:

```cpp
for (int i = 0; i < L; ++i) {
    for (int j = 0; j < L; ++j) {
        if (i == j) continue;
        double Jij = J(i, j);
        if (Jij == 0.0) continue;
        add_J_pair_monomials(i, j, Jij, L, monomials);
    }
}
```

---

## 4. Sparse antisymmetrized $V$ tensor

### 4.1 Expansion

For each monomial $C \, c_p^\dagger c_q^\dagger c_s c_r$ (with $p \ne q$, $r \ne s$), add four sparse entries:

$$V_{pqrs} \mathrel{+}= C, \quad V_{qprs} \mathrel{+}= -C, \quad V_{pqsr} \mathrel{+}= -C, \quad V_{qpsr} \mathrel{+}= C$$

```cpp
struct VEntry {
    int a, b, c, d;
    std::complex<double> val;
};

void add_antisymmetrized(int p, int q, int r, int s,
                          std::complex<double> C,
                          std::vector<VEntry>& V) {
    assert(p != q && r != s);
    V.push_back({p, q, r, s,  C});
    V.push_back({q, p, r, s, -C});
    V.push_back({p, q, s, r, -C});
    V.push_back({q, p, s, r,  C});
}
```

### 4.2 Compression (mandatory, before SCF)

```cpp
void compress_V(std::vector<VEntry>& V, double drop_tol = 1e-14) {
    std::sort(V.begin(), V.end(), [](const VEntry& x, const VEntry& y) {
        return std::tie(x.a, x.b, x.c, x.d) < std::tie(y.a, y.b, y.c, y.d);
    });
    std::vector<VEntry> out;
    for (const auto& e : V) {
        if (!out.empty() && out.back().a == e.a && out.back().b == e.b
            && out.back().c == e.c && out.back().d == e.d) {
            out.back().val += e.val;
        } else {
            out.push_back(e);
        }
    }
    V.clear();
    for (auto& e : out) {
        if (std::abs(e.val) > drop_tol) V.push_back(e);
    }
}
```

### 4.3 Size estimate (16×16 J1-J2)

Unordered NN bonds: $2L = 512$. Unordered NNN bonds: $2L = 512$. Ordered pairs (both directions): $4L = 1024$ each, total $2048$. With 4 monomials per pair and 4 antisymmetrized entries per monomial: $\sim 32{,}768$ raw entries before compression. Compression typically reduces this by a factor of 2-4. The number of distinct $(a,b,c,d)$ keys is comfortably under $10^5$ — trivial compared to the $1024 \times 1024$ diagonalization.

---

## 5. Field construction

```cpp
void build_fields(const std::vector<VEntry>& V,
                   const Eigen::MatrixXcd& rho,
                   const Eigen::MatrixXcd& kappa,
                   Eigen::MatrixXcd& Gamma,
                   Eigen::MatrixXcd& Delta) {
    const int M = static_cast<int>(rho.rows());
    Gamma.setZero(M, M);
    Delta.setZero(M, M);

    for (const auto& e : V) {
        // Gamma_{a,c} += V_{abcd} * rho_{d,b}    (a=p, b=r, c=q, d=s)
        Gamma(e.a, e.c) += e.val * rho(e.d, e.b);
        // Delta_{a,b} += (1/2) V_{abcd} * kappa_{c,d}  (a=p, b=q, c=r, d=s)
        Delta(e.a, e.b) += 0.5 * e.val * kappa(e.c, e.d);
    }

    // Verify symmetries BEFORE cleanup. Bugs hide if you symmetrize blindly.
    const double g_err = (Gamma - Gamma.adjoint()).norm();
    const double d_err = (Delta + Delta.transpose()).norm();
    if (g_err > 1e-9 || d_err > 1e-9) {
        std::ostringstream oss;
        oss << "Field symmetry violation: |Gamma - Gamma^dag| = " << g_err
            << ", |Delta + Delta^T| = " << d_err;
        throw std::runtime_error(oss.str());
    }

    // Clean roundoff-level noise only.
    Gamma = 0.5 * (Gamma + Gamma.adjoint());
    Delta = 0.5 * (Delta - Delta.transpose());
}
```

The pre-cleanup check is non-negotiable: it has caught real bugs in independent HFB codes during development. After full validation passes, you may relax the throw to a warning.

---

## 6. BdG matrix construction

```cpp
Eigen::MatrixXcd build_bdg(const Eigen::MatrixXcd& h_eff,   // h0 + Gamma - mu*I
                            const Eigen::MatrixXcd& Delta) {
    const int M = static_cast<int>(h_eff.rows());
    Eigen::MatrixXcd H(2*M, 2*M);
    H.topLeftCorner(M, M)     =  h_eff;
    H.topRightCorner(M, M)    =  Delta;
    H.bottomLeftCorner(M, M)  = -Delta.conjugate();
    H.bottomRightCorner(M, M) = -h_eff.conjugate();

    const double herm_err = (H - H.adjoint()).norm();
    if (herm_err > 1e-9) {
        throw std::runtime_error("BdG matrix not Hermitian: err = "
                                  + std::to_string(herm_err));
    }
    return H;
}
```

Do not symmetrize the BdG matrix during development. Only after validation should you optionally apply $H \leftarrow \tfrac{1}{2}(H + H^\dagger)$ to suppress accumulated roundoff.

---

## 7. Density extraction

### 7.1 Zero temperature: lowest-$M$ projector

Always take the lowest $M$ eigenvectors. Sign-based selection of "negative-energy" states is fragile near gap closings.

```cpp
struct HFBState {
    Eigen::MatrixXcd rho;
    Eigen::MatrixXcd kappa;
};

HFBState density_zeroT(const Eigen::MatrixXcd& H) {
    const int twoM = static_cast<int>(H.rows());
    const int M = twoM / 2;

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> es(H);
    if (es.info() != Eigen::Success) {
        throw std::runtime_error("BdG diagonalization failed");
    }

    Eigen::MatrixXcd Wlow = es.eigenvectors().leftCols(M);
    Eigen::MatrixXcd R = Wlow * Wlow.adjoint();

    HFBState out;
    out.rho   = R.topLeftCorner(M, M);
    out.kappa = R.topRightCorner(M, M);
    out.rho   = 0.5 * (out.rho + out.rho.adjoint());
    out.kappa = 0.5 * (out.kappa - out.kappa.transpose());
    return out;
}
```

### 7.2 Finite temperature: Fermi-Dirac

Recommended whenever the BdG gap is small or to smooth SCF convergence:

```cpp
HFBState density_finiteT(const Eigen::MatrixXcd& H, double beta) {
    const int twoM = static_cast<int>(H.rows());
    const int M = twoM / 2;

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> es(H);
    if (es.info() != Eigen::Success) {
        throw std::runtime_error("BdG diagonalization failed");
    }

    Eigen::VectorXd f(twoM);
    for (int a = 0; a < twoM; ++a) {
        double x = beta * es.eigenvalues()(a);
        if      (x >  40.0) f(a) = 0.0;
        else if (x < -40.0) f(a) = 1.0;
        else                f(a) = 1.0 / (1.0 + std::exp(x));
    }
    Eigen::MatrixXcd R = es.eigenvectors() * f.asDiagonal()
                       * es.eigenvectors().adjoint();

    HFBState out;
    out.rho   = R.topLeftCorner(M, M);
    out.kappa = R.topRightCorner(M, M);
    out.rho   = 0.5 * (out.rho + out.rho.adjoint());
    out.kappa = 0.5 * (out.kappa - out.kappa.transpose());
    return out;
}
```

Practical schedule: start at $\beta = 50$, anneal up by factors of $\sim 1.5$ once the residual stabilizes, finish at $\beta \gtrsim 500$ or switch to $T = 0$ for the last few iterations.

---

## 8. Energy evaluation

### 8.1 Primary: Wick from monomials

$$E = \sum_\alpha C_\alpha \left[ \rho_{r_\alpha p_\alpha} \rho_{s_\alpha q_\alpha} - \rho_{s_\alpha p_\alpha} \rho_{r_\alpha q_\alpha} + \kappa_{p_\alpha q_\alpha}^* \kappa_{r_\alpha s_\alpha} \right]$$

```cpp
double energy_wick(const std::vector<Monomial>& mons,
                    const Eigen::MatrixXcd& rho,
                    const Eigen::MatrixXcd& kappa) {
    std::complex<double> E = 0.0;
    for (const auto& t : mons) {
        std::complex<double> wick =
            rho(t.r, t.p) * rho(t.s, t.q)
          - rho(t.s, t.p) * rho(t.r, t.q)
          + std::conj(kappa(t.p, t.q)) * kappa(t.r, t.s);
        E += t.coeff * wick;
    }
    if (std::abs(E.imag()) > 1e-8) {
        throw std::runtime_error("Wick energy has imaginary part: "
                                  + std::to_string(E.imag()));
    }
    return E.real();
}
```

### 8.2 Cross-check: trace form

$$E = \mathrm{Tr}(h^0 \rho) + \tfrac{1}{2} \mathrm{Tr}(\Gamma \rho) + \tfrac{1}{2} \sum_{pq} \Delta_{pq} \kappa_{pq}^*$$

```cpp
double energy_trace(const Eigen::MatrixXcd& h0,
                     const Eigen::MatrixXcd& Gamma,
                     const Eigen::MatrixXcd& Delta,
                     const Eigen::MatrixXcd& rho,
                     const Eigen::MatrixXcd& kappa) {
    std::complex<double> E1 = (h0 * rho).trace();
    std::complex<double> E2 = 0.5 * (Gamma * rho).trace();
    std::complex<double> E3 = 0.0;
    const int M = static_cast<int>(rho.rows());
    for (int p = 0; p < M; ++p)
        for (int q = 0; q < M; ++q)
            E3 += 0.5 * Delta(p,q) * std::conj(kappa(p,q));

    std::complex<double> E = E1 + E2 + E3;
    if (std::abs(E.imag()) > 1e-8) {
        throw std::runtime_error("Trace energy has imaginary part");
    }
    return E.real();
}
```

The two energies must agree to $\sim 10^{-8}$ at every iteration. Disagreement indicates an inconsistency between the monomial expansion, the antisymmetrized $V$, or the field formulas.

---

## 9. Chemical potential update

### 9.1 Why no nested $\mu$ solve

A fully nested inner solve $\text{Tr}\,\rho(\mu) = N_\text{target}$ at fixed $\Gamma, \Delta$ over-converges the wrong inner problem: $\Gamma$ depends on $\rho$, so the inner solution becomes inconsistent the moment fields are rebuilt. Use **one** $\mu$ update per outer SCF iteration.

### 9.2 Default: damped secant with bracket fallback

Track the two most recent $(\mu, N)$ pairs.

```cpp
struct MuState {
    double mu_curr = 0.0;
    double mu_prev = 0.0;
    double N_curr  = 0.0;
    double N_prev  = 0.0;
    bool   have_prev = false;
};

double update_mu(MuState& s, double N_target,
                  double max_step = 0.5,
                  double damping = 0.7) {
    double mu_new;
    if (!s.have_prev || std::abs(s.mu_curr - s.mu_prev) < 1e-12) {
        // Proportional fallback. eta tuned to typical bandwidth ~ J*z.
        const double eta = 0.1;
        mu_new = s.mu_curr + eta * (N_target - s.N_curr);
    } else {
        double dN_dmu = (s.N_curr - s.N_prev) / (s.mu_curr - s.mu_prev);
        if (std::abs(dN_dmu) < 1e-6) {
            mu_new = s.mu_curr + 0.1 * (N_target - s.N_curr);
        } else {
            double step = damping * (N_target - s.N_curr) / dN_dmu;
            step = std::clamp(step, -max_step, max_step);
            mu_new = s.mu_curr + step;
        }
    }
    s.mu_prev = s.mu_curr;
    s.N_prev  = s.N_curr;
    s.mu_curr = mu_new;
    s.have_prev = true;
    return mu_new;
}
```

### 9.3 Target

For the J1-J2 spin model in Abrikosov representation at the physical filling: $N_\text{target} = L = 256$ for 16×16. We do not enforce this site-by-site (no local Lagrange multipliers in this build).

---

## 10. SCF loop with Anderson mixing

### 10.1 State packing

Pack $(\rho, \kappa)$ into a real vector for the mixer. Store only independent components: $\rho$ Hermitian → $M^2$ real numbers (real diagonal + complex upper triangle); $\kappa$ antisymmetric → $M(M-1)$ real numbers (complex strict upper triangle). Total: $2M^2 - M$ reals.

```cpp
Eigen::VectorXd pack(const Eigen::MatrixXcd& rho,
                      const Eigen::MatrixXcd& kappa);
void unpack(const Eigen::VectorXd& v,
             Eigen::MatrixXcd& rho, Eigen::MatrixXcd& kappa);
```

After every unpack, project back to the symmetry manifold:

$$\rho \leftarrow \tfrac{1}{2}(\rho + \rho^\dagger), \qquad \kappa \leftarrow \tfrac{1}{2}(\kappa - \kappa^T)$$

### 10.2 Anderson mixer (depth $m$)

Given iterates $X_k$ and SCF map $F$, residual $r_k = F(X_k) - X_k$.

```cpp
class AndersonMixer {
public:
    AndersonMixer(int depth, double alpha)
        : m_(depth), alpha_(alpha) {}

    Eigen::VectorXd step(const Eigen::VectorXd& X,
                          const Eigen::VectorXd& FX) {
        Eigen::VectorXd r = FX - X;
        Xs_.push_back(X);
        Fs_.push_back(FX);
        Rs_.push_back(r);
        while (static_cast<int>(Xs_.size()) > m_ + 1) {
            Xs_.pop_front(); Fs_.pop_front(); Rs_.pop_front();
        }

        const int k = static_cast<int>(Rs_.size()) - 1;
        if (k == 0) return X + alpha_ * r;  // linear step

        // Build dR matrix: columns are r_k - r_{k-1}, ..., r_k - r_{k-m}
        Eigen::MatrixXd dR(r.size(), k);
        Eigen::MatrixXd dF(r.size(), k);
        for (int i = 0; i < k; ++i) {
            dR.col(i) = Rs_[k] - Rs_[k - 1 - i];
            dF.col(i) = Fs_[k] - Fs_[k - 1 - i];
        }
        // Solve least squares dR * gamma = r_k
        Eigen::VectorXd gamma = dR.colPivHouseholderQr().solve(Rs_[k]);
        return Fs_[k] - dF * gamma;
    }

    void reset() { Xs_.clear(); Fs_.clear(); Rs_.clear(); }

private:
    int m_;
    double alpha_;
    std::deque<Eigen::VectorXd> Xs_, Fs_, Rs_;
};
```

### 10.3 Mixing schedule

```cpp
struct MixingPolicy {
    int    linear_warmup_iters  = 5;
    int    anderson_depth       = 6;
    double linear_alpha         = 0.3;
    double max_step_norm        = 0.5;     // safeguard
    int    fallback_after_bad   = 3;       // # bad steps -> reset to linear
};
```

Heuristics:
- Iterations 0 to `linear_warmup_iters - 1`: pure linear mixing.
- After warmup: Anderson with depth 6.
- If energy increases for 3 consecutive Anderson steps, or the residual norm grows by >10×, **reset** the mixer history and run 3 linear steps.
- Cap the Anderson step norm at `max_step_norm` (rescale if exceeded).

### 10.4 Full SCF driver

```cpp
struct HFBParams {
    int    max_iter        = 500;
    double rho_tol         = 1e-7;
    double kappa_tol       = 1e-7;
    double energy_tol      = 1e-9;
    double number_tol      = 1e-5;
    double mu_init         = 0.0;
    double beta            = 0.0;          // 0 means T=0
    bool   verbose         = true;
    MixingPolicy mix;
};

struct HFBResult {
    HFBState state;
    double   energy        = 0.0;
    double   mu            = 0.0;
    double   particle_num  = 0.0;
    int      iterations    = 0;
    bool     converged     = false;
    std::vector<double> energy_history;
    std::vector<double> residual_history;
};

HFBResult run_scf(const AbrikosovHamiltonian& ham,
                   const Eigen::MatrixXcd& h0,
                   HFBState state0,
                   const HFBParams& p) {
    const int M = ham.M;
    HFBState state = state0;
    double mu = p.mu_init;
    MuState mu_state{mu, mu, 0.0, 0.0, false};

    AndersonMixer mixer(p.mix.anderson_depth, p.mix.linear_alpha);
    Eigen::MatrixXcd Gamma(M,M), Delta(M,M);
    double E_prev = std::numeric_limits<double>::infinity();
    int bad_steps = 0;

    HFBResult res;

    for (int it = 0; it < p.max_iter; ++it) {
        // 1) Build mean fields from current densities.
        ham.build_fields(state.rho, state.kappa, Gamma, Delta);

        // 2) Effective single-particle Hamiltonian.
        Eigen::MatrixXcd h_eff = h0 + Gamma;
        h_eff.diagonal().array() -= mu;

        // 3) BdG matrix.
        Eigen::MatrixXcd H = build_bdg(h_eff, Delta);

        // 4) Diagonalize and form output density.
        HFBState candidate = (p.beta > 0.0)
            ? density_finiteT(H, p.beta)
            : density_zeroT(H);

        // 5) Update mu using current N.
        mu_state.N_curr = candidate.rho.trace().real();
        mu = update_mu(mu_state, static_cast<double>(ham.L));

        // 6) Mix.
        Eigen::VectorXd X  = pack(state.rho, state.kappa);
        Eigen::VectorXd FX = pack(candidate.rho, candidate.kappa);
        Eigen::VectorXd Xnew_v;

        if (it < p.mix.linear_warmup_iters) {
            Xnew_v = X + p.mix.linear_alpha * (FX - X);
        } else {
            Xnew_v = mixer.step(X, FX);
            // Safeguard: cap the step size.
            Eigen::VectorXd step = Xnew_v - X;
            double sn = step.norm();
            if (sn > p.mix.max_step_norm) {
                Xnew_v = X + (p.mix.max_step_norm / sn) * step;
            }
        }

        Eigen::MatrixXcd rho_new(M,M), kappa_new(M,M);
        unpack(Xnew_v, rho_new, kappa_new);
        rho_new   = 0.5 * (rho_new + rho_new.adjoint());
        kappa_new = 0.5 * (kappa_new - kappa_new.transpose());

        // 7) Diagnostics.
        double drho   = (rho_new - state.rho).norm();
        double dkappa = (kappa_new - state.kappa).norm();
        double E_wick  = ham.energy_wick(rho_new, kappa_new);
        double E_trace = energy_trace(h0, Gamma, Delta, rho_new, kappa_new);
        if (std::abs(E_wick - E_trace) > 1e-6) {
            throw std::runtime_error("Energy cross-check failed: dE = "
                + std::to_string(E_wick - E_trace));
        }
        double Ncur = rho_new.trace().real();
        double dE   = E_wick - E_prev;

        res.energy_history.push_back(E_wick);
        res.residual_history.push_back(std::sqrt(drho*drho + dkappa*dkappa));

        if (p.verbose) {
            std::printf("iter %4d  E=%.10f  dE=%+.3e  drho=%.3e  "
                        "dkap=%.3e  N=%.4f  mu=%.4f\n",
                        it, E_wick, dE, drho, dkappa, Ncur, mu);
        }

        // 8) Mixer reset on pathology.
        if (it >= p.mix.linear_warmup_iters) {
            if (dE > 1e-6) {
                ++bad_steps;
                if (bad_steps >= p.mix.fallback_after_bad) {
                    if (p.verbose) std::puts("  [mixer reset]");
                    mixer.reset();
                    bad_steps = 0;
                }
            } else {
                bad_steps = 0;
            }
        }

        state.rho = rho_new;
        state.kappa = kappa_new;

        // 9) Convergence.
        if (drho < p.rho_tol && dkappa < p.kappa_tol
            && std::abs(dE) < p.energy_tol
            && std::abs(Ncur - ham.L) < p.number_tol) {
            res.converged = true;
            res.iterations = it + 1;
            break;
        }
        E_prev = E_wick;
        res.iterations = it + 1;
    }

    res.state = state;
    res.energy = ham.energy_wick(state.rho, state.kappa);
    res.mu = mu;
    res.particle_num = state.rho.trace().real();
    return res;
}
```

---

## 11. Initialization

### 11.1 Why $\kappa \ne 0$ is required

If $\kappa_0 = 0$, then $\Delta = 0$ identically, and every iterate stays in the normal sector. To explore the paired/spin-liquid landscape, seed nonzero pairing.

### 11.2 $S_z$-preserving pairing seeds

Since the J1-J2 Hamiltonian conserves $S_z$ and we expect singlet/$S_z = 0$ pairing to dominate, seed only $\kappa_{i\uparrow, j\downarrow}$ (and the antisymmetric partner $\kappa_{j\downarrow, i\uparrow} = -\kappa_{i\uparrow, j\downarrow}$). The SCF then stays in this sector up to numerical noise.

```cpp
HFBState init_singlet_seed(int L, const Eigen::MatrixXd& J,
                            double amplitude,
                            const std::string& pattern,
                            std::mt19937& rng) {
    const int M = 2 * L;
    HFBState s;
    s.rho   = 0.5 * Eigen::MatrixXcd::Identity(M, M);  // half filling
    s.kappa = Eigen::MatrixXcd::Zero(M, M);

    auto sign = [&](int i, int j) -> double {
        if (pattern == "uniform")    return 1.0;
        if (pattern == "dwave") {
            // d_{x^2-y^2}: +1 on x-bonds, -1 on y-bonds (NN only).
            int xi = i % /*Lx=*/16, yi = i / 16;
            int xj = j % 16,        yj = j / 16;
            int dx = std::abs(xi - xj); int dy = std::abs(yi - yj);
            // periodic wrap not handled here; flesh out for production
            if (dy == 0) return  1.0;
            if (dx == 0) return -1.0;
            return 1.0;
        }
        if (pattern == "random") {
            std::uniform_real_distribution<double> u(-1.0, 1.0);
            return u(rng);
        }
        return 1.0;
    };

    for (int i = 0; i < L; ++i) {
        for (int j = 0; j < L; ++j) {
            if (i == j) continue;
            if (J(i,j) == 0.0) continue;
            int p_up = so(i, UP, L);
            int q_dn = so(j, DOWN, L);
            std::complex<double> v = amplitude * sign(i, j);
            // kappa_{p,q} = <c_q c_p>; antisymmetry: kappa_{q,p} = -kappa_{p,q}
            s.kappa(p_up, q_dn) =  v;
            s.kappa(q_dn, p_up) = -v;
        }
    }
    s.rho   = 0.5 * (s.rho + s.rho.adjoint());
    s.kappa = 0.5 * (s.kappa - s.kappa.transpose());
    return s;
}
```

### 11.3 Recommended workflow

Run the SCF from multiple seeds and pick the lowest-energy converged solution:

| Seed | Pattern | Amplitude |
|---|---|---|
| 1 | `uniform` (NN) | $10^{-2}$ |
| 2 | `dwave` (NN) | $10^{-2}$ |
| 3 | `uniform` (NN+NNN) | $10^{-2}$ |
| 4 | `random` (NN+NNN) | $10^{-3}$ × 5 different seeds |

For frustrated J1-J2, multiple metastable HFB minima are expected. This is a feature: their relative energies and order parameters tell you about the mean-field phase diagram.

---

## 12. Validation suite (mandatory before any 16×16 run)

### 12.1 Two-site exact spectrum

Build the full 16-dimensional Fock-space Hamiltonian for $L=2$ from the monomial list (no HFB), diagonalize, and verify that the singly-occupied subspace gives the spin-1/2 dimer spectrum.

For physical coupling $J = 1$ (so input matrix $J_{12} = J_{21} = 1/2$ in `SymmetricPhysicalBonds` mode), the spin Hamiltonian is $\mathbf{S}_1 \cdot \mathbf{S}_2$ with eigenvalues:

$$E_\text{singlet} = -\tfrac{3}{4}, \qquad E_\text{triplet} = +\tfrac{1}{4} \text{ (×3)}$$

```cpp
TEST(TwoSite, SpinSpectrum) {
    int L = 2;
    Eigen::MatrixXd J = Eigen::MatrixXd::Zero(2,2);
    J(0,1) = 0.5; J(1,0) = 0.5;          // physical J = 1

    auto ham = AbrikosovHamiltonian::build(J);
    auto Hfock = build_full_fock_hamiltonian(ham);  // 16x16 dense

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> es(Hfock);
    auto evals = filter_singly_occupied_eigenvalues(es, ham);
    std::sort(evals.begin(), evals.end());

    EXPECT_NEAR(evals[0], -0.75, 1e-12);
    EXPECT_NEAR(evals[1],  0.25, 1e-12);
    EXPECT_NEAR(evals[2],  0.25, 1e-12);
    EXPECT_NEAR(evals[3],  0.25, 1e-12);
}
```

This single test catches: monomial sign errors, factor-of-two errors in bond counting, index-decoding errors in $(p,q,r,s)$, and basic operator-ordering mistakes. **No HFB result on a larger lattice should be trusted until this passes.**

### 12.2 Tensor antisymmetry

```cpp
TEST(VTensor, FullAntisymmetry) {
    auto V = ham.V;
    // Build a fast lookup: (a,b,c,d) -> val
    std::map<std::tuple<int,int,int,int>, std::complex<double>> Vmap;
    for (auto& e : V) Vmap[{e.a,e.b,e.c,e.d}] = e.val;

    for (auto& [key, val] : Vmap) {
        auto [a,b,c,d] = key;
        auto check = [&](int x, int y, int z, int w, std::complex<double> expected) {
            auto it = Vmap.find({x,y,z,w});
            std::complex<double> v = (it != Vmap.end()) ? it->second : 0.0;
            EXPECT_NEAR(std::abs(v - expected), 0.0, 1e-12);
        };
        check(b,a,c,d, -val);
        check(a,b,d,c, -val);
        check(b,a,d,c,  val);
    }
}
```

### 12.3 Field symmetries (pre-cleanup)

For random valid $\rho, \kappa$ inputs, check $\Gamma = \Gamma^\dagger$ and $\Delta = -\Delta^T$ to better than $10^{-9}$ before any symmetrization is applied.

### 12.4 BdG Hermiticity

$\|\mathcal{H} - \mathcal{H}^\dagger\| < 10^{-9}$ at every iteration, before any cleanup.

### 12.5 Eigenvalue pairing

Sort eigenvalues ascending; for each $a \in [0, M)$:

$$|E_a + E_{2M - 1 - a}| < 10^{-8}$$

### 12.6 Density symmetries

After every density extraction:
- $\|\rho - \rho^\dagger\| < 10^{-10}$
- $\|\kappa + \kappa^T\| < 10^{-10}$

### 12.7 Idempotency at $T=0$ before mixing

Right after `density_zeroT`, $\|R^2 - R\|_\text{F} < 10^{-9}$.

### 12.8 Energy cross-check

At every iteration, $|E_\text{Wick} - E_\text{trace}| < 10^{-8}$. Disagreement is a hard error.

### 12.9 Number constraint

At convergence, $|\mathrm{Tr}\,\rho - L| < 10^{-5}$.

### 12.10 Free-fermion limit

Set $V = 0$ (zero out all monomials), put a small NN hopping in $h^0$, run the solver, and compare to a direct $h^0$ diagonalization. Energies and $\rho$ must match.

### 12.11 Small-lattice ED comparison

For $L \in \{4, 6, 8\}$ (a $1 \times L$ chain with PBC, projected to single occupancy), compare the HFB energy density to ED in the singly occupied sector. HFB will be variational — energy should be above ED. The gap quantifies mean-field error.

---

## 13. Complexity, memory, and performance

### 13.1 Per-iteration cost

| Component | Cost | Wall-time share (16×16) |
|---|---|---|
| `build_fields` | $O(\|V\|) \approx 10^4$ ops | <1% |
| `build_bdg` | $O(M^2)$ copies | <1% |
| Diagonalization (`zheevd`, $2M = 1024$) | $O((2M)^3) \approx 10^9$ | ~95% |
| Form $R = W W^\dagger$ ($M$ vecs) | $O(M^3)$ | ~3% |
| Mixer + diagnostics | $O(M^2)$ | <1% |

Diagonalization dominates. Use multithreaded LAPACK (MKL or OpenBLAS) and link against a complex-double-tuned build.

### 13.2 Memory

| Object | Size (16×16) |
|---|---|
| $\rho, \kappa, \Gamma, \Delta, h$ | $5 \times 512^2 \times 16$ B $\approx 21$ MB |
| BdG matrix | $1024^2 \times 16$ B $= 16$ MB |
| Eigenvectors (full) | $1024^2 \times 16$ B $= 16$ MB |
| Anderson history (depth 6, packed reals) | $\sim 7 \times 5 \times 10^5 \times 8$ B $\approx 30$ MB |
| Sparse $V$ entries | $\lesssim 10^5 \times 24$ B $= 2.4$ MB |

Total well under 200 MB. Fits comfortably in cache-conscious workstation memory.

### 13.3 Diagonalization backend

Prototype:
```cpp
Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> es(H);
```
Production (with `LAPACKE`):
```cpp
LAPACKE_zheevd(LAPACK_COL_MAJOR, 'V', 'U', n,
               reinterpret_cast<lapack_complex_double*>(H.data()),
               n, evals.data());
```
Use `'V'` to compute eigenvectors. Use `'U'` (upper triangle) to match Eigen's storage. For very large lattices later, `zheevr` with a range query on the lowest $M$ eigenvalues would save the upper-half work, but for 16×16 the savings are not worth the complexity.

---

## 14. Project layout

```
hfb_j1j2/
├── CMakeLists.txt
├── include/hfb/
│   ├── lattice.hpp          # spatial geometry, J1J2 matrix builder
│   ├── monomial.hpp         # Monomial, VEntry types
│   ├── hamiltonian.hpp      # AbrikosovHamiltonian
│   ├── fields.hpp           # build_fields, build_bdg
│   ├── density.hpp          # density_zeroT, density_finiteT
│   ├── energy.hpp           # energy_wick, energy_trace
│   ├── mu.hpp               # MuState, update_mu
│   ├── mixer.hpp            # AndersonMixer
│   ├── init.hpp             # init_singlet_seed
│   ├── solver.hpp           # HFBParams, HFBResult, run_scf
│   └── diagnostics.hpp      # validation checks
├── src/
│   ├── lattice.cpp
│   ├── hamiltonian.cpp
│   ├── fields.cpp
│   ├── density.cpp
│   ├── solver.cpp
│   ├── mixer.cpp
│   ├── diagnostics.cpp
│   └── main.cpp
├── tests/
│   ├── test_two_site.cpp        # MUST PASS FIRST
│   ├── test_v_antisymmetry.cpp
│   ├── test_field_symmetry.cpp
│   ├── test_bdg_hermiticity.cpp
│   ├── test_density_idempotent.cpp
│   ├── test_energy_crosscheck.cpp
│   ├── test_free_fermion_limit.cpp
│   └── test_small_chain_vs_ed.cpp
└── README.md
```

Build with CMake; depend on Eigen3 and (optionally) a LAPACK provider (MKL preferred).

---

## 15. Out-of-scope / deferred features

Explicitly **not** in this build (each is a future extension):

1. **Local single-occupancy Lagrange multipliers.** The user has stated this is acceptable; only the global parton number is constrained.
2. **DIIS.** Anderson mixing covers the same ground for HFB and is simpler to make robust.
3. **Explicit $S_z$-block-reduced BdG.** Use full unrestricted BdG with $S_z$-preserving initialization. The block reduction is a $\sim 4\times$ speedup that can be added once the unrestricted code is verified correct.
4. **Symmetry-projected HFB / number projection.** Possible later, but requires storing the full $W$ matrix and doing post-projection integrals.
5. **Iterative sparse eigensolvers.** HFB needs half the spectrum; dense $\texttt{zheevd}$ is the right choice for 1024×1024.
6. **Sparse supermatrix representation of $\rho \mapsto \Gamma$.** The direct sparse-tensor loop is equivalent and simpler.

---

## 16. Summary of conventions and corrections

### What changed from the original plan

| Item | Original | Final |
|---|---|---|
| Wick anomalous sign | $+\kappa_{pq}^* \kappa_{rs}$ | $+\kappa_{pq}^* \kappa_{rs}$ ✓ (confirmed correct under stated convention) |
| $\Delta$ definition | $\frac{1}{2} \sum V_{pqrs} \kappa_{rs}$ | unchanged ✓ |
| Bond-counting | implicit | explicit `JConvention` enum + halving helper |
| Eigenvalue selection | "negative-energy" with fallback | always lowest $M$ eigenvectors |
| $\mu$ update | nested inner solve | one secant step per outer SCF iteration |
| Mixing | linear, defer DIIS | linear warmup → Anderson depth 6 with safeguards |
| $V$ entries | uncompressed | sort + dedupe before SCF |
| Symmetrization | unconditional | check first, throw on violation, then clean |
| Energy | Wick only | Wick + trace cross-check at every iteration |
| Spin sector | full unrestricted | full unrestricted with $S_z$-preserving init |
| Two-site test | mentioned | **mandatory gate** before larger runs |

### Key conventions, restated

$$\rho_{pq} = \langle c_q^\dagger c_p \rangle, \qquad \kappa_{pq} = \langle c_q c_p \rangle$$

$$H_2 = \tfrac{1}{4} \sum_{pqrs} V_{pqrs} \, c_p^\dagger c_q^\dagger c_s c_r, \quad V \text{ fully antisymmetric in } (p,q) \text{ and in } (r,s)$$

$$\Gamma_{pq} = \sum_{rs} V_{prqs} \rho_{sr}, \qquad \Delta_{pq} = \tfrac{1}{2} \sum_{rs} V_{pqrs} \kappa_{rs}$$

$$\langle c_p^\dagger c_q^\dagger c_s c_r \rangle = \rho_{rp} \rho_{sq} - \rho_{sp} \rho_{rq} + \kappa_{pq}^* \kappa_{rs}$$

$$\mathcal{H} = \begin{pmatrix} h - \mu I & \Delta \\ -\Delta^* & -(h-\mu I)^* \end{pmatrix}, \qquad R_\text{T=0} = W_\text{low} W_\text{low}^\dagger$$

These four equations are the contract that every component of the code must satisfy. The two-site test verifies the contract end-to-end.

---

## 17. Implementation order (suggested)

1. **Day 1.** Lattice, `Monomial`, `J1J2_matrix` builder. Two-site Fock-space test passes.
2. **Day 2.** `VEntry` expansion, compression, antisymmetry tests.
3. **Day 3.** `build_fields`, field-symmetry tests on random densities.
4. **Day 4.** `build_bdg`, `density_zeroT`, idempotency and pairing tests.
5. **Day 5.** `energy_wick`, `energy_trace`, cross-check test.
6. **Day 6.** `MuState`, linear-mixing SCF loop. Free-fermion-limit test passes.
7. **Day 7.** Anderson mixer, mixer-reset safeguards.
8. **Day 8.** Initialization helpers, multi-seed driver.
9. **Day 9.** Small-chain vs ED comparison; tighten convergence tolerances.
10. **Day 10+.** First 16×16 production runs across $J_2/J_1 \in [0, 1]$.

Do **not** skip ahead. The validation gates exist because HFB sign and convention errors are silent: the code will run, converge to *something*, and produce plausible-looking energies that are physically wrong.
