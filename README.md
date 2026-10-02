# prolate-weil-bounds

This repository contains the computations for two papers by Ryota Mori, together with their LaTeX sources and PDFs.

**Paper I** (`paper-I/`). *A non-asymptotic bound with the sharp power of c for the eigenvalue defect 1 − λ_n(c) of time and frequency limiting.*

For 0 ≤ n ≤ 4 and every c ≥ 10π,

  1 − λ_n(c) ≤ C_n c^{n+1/2} e^{−2c},

with C_n = 18.34, 155.0, 703.4, 2340 and 6621. The power of c is the one in Fuchs' asymptotic formula, and the constants are within a factor 2.6–5.5 of Fuchs' constants.

**Paper II** (`paper-II/`). *The Weil quadratic form of the Connes–Consani–Moscovici prolate vector, and a μ⁸ upper bound for windowed Weil forms.* All three results are unconditional:
- **Theorem C:** |W(k_λ)| ≤ P(μ)ρ_out‖h_λ‖² + e^{−16μ}‖h_λ‖².
- **Theorem D:** ‖k_λ‖² ≥ 0.109‖h_λ‖² for μ ≥ 50.
- **Theorem E:** λ_min(λ) ≤ 5.011×10¹⁷ μ⁸ (log μ)³ e^{−4πμ} for μ ≥ 50.

Paper II is a sequel to *Unconditional doubly exponential upper bounds for the bottom of windowed Weil quadratic forms* ([code](https://github.com/moriryota/weil-window-upper-bounds), [paper](https://doi.org/10.5281/zenodo.23059297)).

DOIs:
- Code (all versions): [10.5281/zenodo.23092422](https://doi.org/10.5281/zenodo.23092422)
- Paper I: [10.5281/zenodo.23092454](https://doi.org/10.5281/zenodo.23092454)
- Paper II: [10.5281/zenodo.23092542](https://doi.org/10.5281/zenodo.23092542)

These are upper bounds only. They say nothing about positivity and do not address the Riemann Hypothesis.

## Layout

| Directory | Contents |
|---|---|
| `paper-I/`, `paper-II/` | LaTeX sources and PDFs. |
| `proofs-I/`, `proofs-II/` | Computations that **prove** the numerical constants of the papers. They use ball arithmetic (Arb, through python-flint) or exact rationals. Every printed upper bound is rounded upward, every lower bound downward, and each is asserted against the ball. |
| `numerics-I/`, `numerics-II/` | Numerical illustrations. These scripts are **not** part of any proof. |
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

## Licence

- Code: Apache License 2.0 (`LICENSE`, `NOTICE`).
- Papers: Creative Commons Attribution 4.0 International (CC BY 4.0).

## Use of generative AI

This work was carried out with extensive use of generative AI systems, namely Anthropic Claude, Google Gemini and OpenAI GPT. Under the direction of the author, they drafted the proofs and the code, ran the computations, reviewed one another's work adversarially, and drafted the manuscripts. The author chose the problems and the strategy, assigned and reviewed the tasks, and takes full responsibility for the content. The AI systems are not authors. The papers have not yet been reviewed by independent human experts.
