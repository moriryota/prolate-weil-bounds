# prolate-weil-bounds

This repository contains the computations for eight papers by Ryota Mori, together with their LaTeX sources and PDFs.

**Paper I** (`paper-I/`). *A non-asymptotic bound with the sharp power of c for the eigenvalue defect 1 − λ_n(c) of time and frequency limiting.*

For 0 ≤ n ≤ 4 and every c ≥ 10π,

  1 − λ_n(c) ≤ C_n c^{n+1/2} e^{−2c},

with C_n = 18.34, 155.0, 703.4, 2340 and 6621. The power of c is the one in Fuchs' asymptotic formula, and the constants are within a factor 2.6–5.5 of Fuchs' constants.

**Paper II** (`paper-II/`). *The Weil quadratic form of the Connes–Consani–Moscovici prolate vector, and a μ⁸ upper bound for windowed Weil forms.* All three results are unconditional:
- **Theorem C:** |W(k_λ)| ≤ P(μ)ρ_out‖h_λ‖² + e^{−16μ}‖h_λ‖².
- **Theorem D:** ‖k_λ‖² ≥ 0.109‖h_λ‖² for μ ≥ 50.
- **Theorem E:** λ_min(λ) ≤ 5.046×10¹⁷ μ⁸ (log μ)³ e^{−4πμ} for μ ≥ 50 (version 2; see the corrections below).

**Paper III** (`paper-III/`). *An endpoint-smoothed prolate vector and a μ^{9/2}(log μ)⁴ upper bound for windowed Weil forms.* The vector of paper II is smoothed in a layer of width ≍ c⁻² at its truncation points. Unconditionally:
- **Theorem F:** |W(k̃_λ)| ≤ 6.048×10¹⁴ (log μ)⁴ B for 50 ≤ μ ≤ 10⁴, and |W(k̃_λ)| ≤ 5.2×10¹⁴ √μ (log μ)⁵ B for all μ ≥ 50, where B = ρ_out‖h_λ‖².
- **Theorem G:** λ_min(λ) ≤ 1.437×10²³ μ^{9/2} (log μ)⁴ e^{−4πμ} for 50 ≤ μ ≤ 10⁴, and λ_min(λ) ≤ 1.24×10²³ μ⁵ (log μ)⁵ e^{−4πμ} for all μ ≥ 50. An explicit refinement is below 0.832 times the bound of paper II, Theorem E, for all μ ≥ 50.

**Paper IV** (`paper-IV/`). *A conditional lower bound for windowed Weil forms from sampling on the zeros of ζ.* Under the Riemann Hypothesis and a local hypothesis (LP) on close pairs of zeros, 4πμ − O(log μ) ≤ −log λ_min(λ) ≤ Cμ log log μ; under RH alone, effective but weak lower bounds. The proofs are analytic; `numerics-IV/` contains numerical illustrations only.

**Paper V** (`paper-V/`). *A μ^{9/2} log μ upper bound for windowed Weil forms.* Two of the three sources of the logarithms in paper III are removed: a translation in the Sobolev step, and the exact Prüfer phase with an L² (non-stationary phase) estimate of the Poisson sum, which replaces paper II, Lemma 4.1 by ‖r_τ‖² ≤ λ^{2τ}·729B/(1−2τ). Unconditionally:
- **Theorem J:** |W(k̃_λ)| ≤ 1.47×10¹⁶ (log μ) B for 50 ≤ μ ≤ 10⁴, and ≤ 3.19×10¹⁶ √μ (log μ)³ B for all μ ≥ 50.
- **Theorem K:** λ_min(λ) ≤ 3.50×10²⁴ μ^{9/2} log μ e^{−4πμ} for 50 ≤ μ ≤ 10⁴, and ≤ 7.58×10²⁴ μ⁵ (log μ)³ e^{−4πμ} for all μ ≥ 50.
- Numerically (not certified), λ_min/(1−λ₄(2πμ)) is about 4.1–5.4 at the integers 7 ≤ μ ≤ 20 (`numerics-V/`).

**Paper VI** (`paper-VI/`, `certificate-VI/`). *Weil positivity on the window of width log 10: an interval-arithmetic certificate.* Unconditionally, Q(f) ≥ 2^−49162 ‖f‖² for every complex C² function f supported in [−½ log 10, ½ log 10]. The proof follows Liu's bounded-comparison framework and reduces positivity to a finite matrix inequality (N = 640), certified in Arb with all quadrature, rounding and inversion errors paid. It is an external interval-arithmetic certificate, not a formal proof. `certificate-VI/` contains the code and all data (about 225 MB; the files are byte-identical to the author's working tree, so the hashes pinned inside the certificate stay valid); see `certificate-VI/README.md` for replaying the check.

**Paper VII** (`paper-VII/`, `certificate-VII/`). *Weil positivity on the window of width log 12, and a lower bound for the prolate eigenvalue defect.* Unconditionally:
- **Theorem M:** Q(f) ≥ 2^−1518 ‖f‖² for every complex C² function f supported in [−½ log 12, ½ log 12] (finite matrix inequality at N = 832, certified in Arb with all errors paid; an external interval-arithmetic certificate, not a formal proof).
- **Theorem N:** 1 − λ₀(c) ≥ e^{−2(c+1)}/(4(c+1)(c+4)) for every c > 0 (an analytic proof). It replaces the out-of-band estimate of paper VI and improves its constant to 2^−876 without recomputing the finite matrices.
`certificate-VII/` contains the code and all data (about 350 MB); see `certificate-VII/README.md`.

**Paper VIII** (`paper-VIII/`, `certificate-VIII/`). *Weil positivity on the window of width log 15.* Unconditionally, Q(f) ≥ 2^−2904 ‖f‖² for every complex C² function f supported in [−½ log 15, ½ log 15] (the prime 13 enters; finite matrix inequality at N = 1664, certified in Arb with all errors paid; an external interval-arithmetic certificate, not a formal proof). `certificate-VIII/` contains the code only; the matrices (about 1.8 GB) are in the separate repository [weil-log15-certificate-data](https://github.com/moriryota/weil-log15-certificate-data) ([10.5281/zenodo.23281176](https://doi.org/10.5281/zenodo.23281176)); see `certificate-VIII/README.md`.

Paper II is a sequel to *Unconditional doubly exponential upper bounds for the bottom of windowed Weil quadratic forms* ([code](https://github.com/moriryota/weil-window-upper-bounds), [paper](https://doi.org/10.5281/zenodo.23059297)).

DOIs:
- Code (all versions): [10.5281/zenodo.23092422](https://doi.org/10.5281/zenodo.23092422)
- Paper I: [10.5281/zenodo.23092454](https://doi.org/10.5281/zenodo.23092454)
- Paper II: version 2 [10.5281/zenodo.23162895](https://doi.org/10.5281/zenodo.23162895); version 1 [10.5281/zenodo.23092542](https://doi.org/10.5281/zenodo.23092542)
- Paper III: version 2 [10.5281/zenodo.23162896](https://doi.org/10.5281/zenodo.23162896); version 1 [10.5281/zenodo.23119609](https://doi.org/10.5281/zenodo.23119609)
- Paper IV: [10.5281/zenodo.23134085](https://doi.org/10.5281/zenodo.23134085)
- Paper V: [10.5281/zenodo.23175061](https://doi.org/10.5281/zenodo.23175061)
- Paper VI: [10.5281/zenodo.23238480](https://doi.org/10.5281/zenodo.23238480) (certificate data in `certificate-VI/` of this repository)
- Paper VII: version 2 [10.5281/zenodo.23265728](https://doi.org/10.5281/zenodo.23265728); version 1 [10.5281/zenodo.23261744](https://doi.org/10.5281/zenodo.23261744) (certificate data in `certificate-VII/` of this repository)
- Paper VIII: [10.5281/zenodo.23281236](https://doi.org/10.5281/zenodo.23281236) (certificate code in `certificate-VIII/`; data [10.5281/zenodo.23281176](https://doi.org/10.5281/zenodo.23281176))

## Changes in v1.5.1 (paper VII, version 2)

`certificate-VII/SHA256SUMS` listed two Python cache files that are not shipped, so the integrity command of v1.5.0 exited with status 1 although every shipped file matched; it now lists exactly the shipped files. The mutation tests are replayed with `certificate-VII/replay_mutations.py` in a disposable copy (running the original script directly stops at the end, because it refuses to overwrite the shipped results). No data file, checker or result changed. Paper VII version 2 adds three clarifications (the functions b_j of Lemma 5.1, the sign of the transform in Lemma 3.2, the identity behind W(X) ⪰ G); its theorems and constants are unchanged.

## Corrections in v1.2.1 (papers II and III, version 2)

Versions 1 of papers II and III used the zero-counting bound of T. Trudgian (J. Number Theory 134 (2014) 280–292, Cor. 1) with the constants 0.111, 0.275, 2.450 of its arXiv version (arXiv:1208.5846v2), while citing the published version, whose constants are 0.112, 0.278, 2.510. In v1.2.1 every script uses the published constants (`proofs-II/theoremC_constants.py`, `proofs-III/*.py`), and the outputs in `outputs/` were regenerated. Changed constants:
- paper II: P(5) ≤ 1.872×10¹³ (was 1.854×10¹³), P(15) ≤ 1.159×10¹⁵, P(100) ≤ 1.568×10¹⁸, Theorem E 5.046×10¹⁷ (was 5.011×10¹⁷);
- paper III: Theorem F(a) 6.048×10¹⁴ (was 6.034×10¹⁴), Theorem G(a) 1.437×10²³ (was 1.434×10²³), C_M and C_T, Table 1. Theorems F(b), G(b) and G(c) keep their printed constants.

The powers of μ and log μ, the statements and the proofs are unchanged. Papers I and IV are not affected.

Papers I–III and V give upper bounds only and say nothing about positivity; papers VI–VIII prove positivity on fixed windows only. None of them addresses the Riemann Hypothesis. Paper IV gives lower bounds under RH (and, for the main theorem, an additional hypothesis); it does not prove RH.

## Layout

| Directory | Contents |
|---|---|
| `paper-I/` … `paper-VIII/` | LaTeX sources and PDFs. |
| `certificate-VI/`, `certificate-VII/`, `certificate-VIII/` | Code and data of the positivity certificates of papers VI and VII; code of paper VIII (its data are in a separate repository). |
| `proofs-I/`, `proofs-II/`, `proofs-III/`, `proofs-V/` | Computations that **prove** the numerical constants of the papers. They use ball arithmetic (Arb, through python-flint) or exact rationals. Every printed upper bound is rounded upward, every lower bound downward, and each is asserted against the ball. |
| `numerics-I/` … `numerics-V/` | Numerical illustrations. These scripts are **not** part of any proof. |
| `common/` | Shared code for the numerics: prolate spheroidal wave functions from Legendre expansions, and the Weil form on a Legendre basis in Arb. |
| `outputs/` | The raw outputs of all scripts, with the same file names that `run_all.sh` produces. |

## Running

- Use Python 3.12 and run `pip install -r requirements.txt`. The versions listed there are the ones that produced `outputs/`.
- Run everything from the repository root.
  - `./run_all.sh` runs the proof scripts (about 30 minutes) and compares each output with `outputs/`.
  - `./run_all.sh all` also runs the numerical illustrations. This takes several hours; the largest Galerkin computation alone takes about 50 minutes.
  - Set `PYTHON=/path/to/python` to choose the interpreter.
- Some proof scripts also write a JSON file next to themselves. In those files, timing fields change from run to run.

**Trial polynomials.** `proofs-I/trials.json` contains the trial polynomials f₀,…,f₄ of paper I, Lemma 5.1, as exact decimal rationals.
- `proofs-I/certify.py` reads `trials.json` and writes the certified matrices and bounds to `proofs-I/certificate.json`.
- `proofs-I/make_trials.py` documents how the trials were chosen. It is not part of the proof, and on another platform it may produce slightly different (equally valid) trials.

**Galerkin values of λ_min (paper V).** `numerics-V/galerkin_lambda_min.py mu N prec nq` computes the bottom of the even part. Convergence in N needs N/(a c) ≳ 2.3 (a = ½ log μ, c = 2πμ), and the precision must grow with N (μ = 20, N = 570 gives a negative value with 840 bits and 1.09899×10⁻⁹⁶ with 1300 bits; `numerics-V/precision_check.py`). The raw outputs of all runs (hours each, on a 4-vCPU machine) are in `outputs/numerics-V/runs/`, and `outputs/numerics-V/ratio_table.txt` collects the ratios. The λ_min(even) values printed by `numerics-II/galerkin_W.py` are not converged at the N listed there; only its W(k)/‖k‖² values were used in paper II, and even these change by a few percent at larger N.

**Galerkin values.** The values of W(k_λ)/‖k_λ‖² in the table of paper II come from `numerics-II/galerkin_W.py`. They need the parameters recorded there, in particular the quadrature size nq = N + 40 for μ ≥ 13. With smaller N or nq the digits are wrong.

## Where each constant is proved

| Paper | Statement | Script |
|---|---|---|
| I | Lemma 2.2: χ_n/c² < 0.3 | `proofs-I/chi_bound.py` |
| I | Lemma 3.1: D₀, D₁, D₂ | `proofs-I/constants_proof.py` |
| I | Lemma 4.1 (BL1) | `proofs-I/bl1_independent.py` |
| I | Lemma 4.1 (BL2) | `proofs-I/bessel_compact.py` |
| I | Lemma 4.1 (BL3); Theorem 2 (factor 1 + 15/c²) | `proofs-I/constants_proof.py`, `proofs-I/assembly_constant.py`, `proofs-I/algebra_checks.py` |
| I | Lemma 5.1: starting values U_n | `proofs-I/certify.py` (with `trials.json`) |
| I | Lemma 2.3 and Theorem 1: τ_n, β_n, C_n | `proofs-I/assemble.py` |
| II | Lemma 3.1 (eigenvalues, variation) | analytic; arithmetic checks in `proofs-II/S1e_selfcheck.py` |
| II | Lemma 3.2 (a), (c) | `proofs-II/chi_bound.py`, `proofs-II/lg_delta_bound.py` |
| II | Lemma 3.3 (endpoint value) | `proofs-II/L5b_constants.py` |
| II | Lemmas 3.4–3.6, 4.1 and P(μ) (Theorem C) | `proofs-II/theoremC_constants.py` |
| II | κ_∞ | `proofs-II/kappa_inf_enclosure.py` |
| II | Theorem D | `proofs-II/L6_effective_constants.py` |
| II | Theorem E | `proofs-II/theoremE_constants.py` |
| III | Lemmas 3.1, 3.2, Proposition 3.3, m₀(c) in (3.3) | `proofs-III/layer_constants.py` |
| III | Lemma 5.1 (trace J); zero-count constants | `proofs-III/trace_constants.py` |
| III | Lemma 4.1 (symbolic identities) | `proofs-III/cancellation_identities.py` |
| III | Lemma 4.1, Propositions 4.2 and 5.3; Theorems F(a), G(a) | `proofs-III/theoremF_a_constants.py` |
| III | P″(μ); Theorems F(b), G(b), G(c) | `proofs-III/theoremF_b_constants.py` |
| III | Sections 6–7 (symbolic identities) | `proofs-III/all_mu_identities.py` |

## Licence

- Code: Apache License 2.0 (`LICENSE`, `NOTICE`).
- Papers: Creative Commons Attribution 4.0 International (CC BY 4.0).

## Use of generative AI

This work was carried out with extensive use of generative AI systems, namely Anthropic Claude, Google Gemini and OpenAI GPT. Under the direction of the author, they drafted the proofs and the code, ran the computations, reviewed one another's work adversarially, and drafted the manuscripts. The author chose the problems and the strategy, assigned and reviewed the tasks, and takes full responsibility for the content. The AI systems are not authors. The papers have not yet been reviewed by independent human experts.
