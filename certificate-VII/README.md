# Weil positivity on the window of width log 12 — certificate data and code

Certificate data and code for: R. Mori, *Weil positivity on the window of width log 12, and a lower bound for the prolate eigenvalue defect*, Zenodo, 2026, doi:10.5281/zenodo.23261744 (paper VII).
This folder is part of https://github.com/moriryota/prolate-weil-bounds (release v1.5.1; first released in v1.5.0), archived on Zenodo, doi:10.5281/zenodo.23092422 (all versions). The data files are included.

Theorem M: for every complex C² function f supported in [−½ log 12, ½ log 12], Weil's quadratic form satisfies
Q(f) ≥ 2^−1518 ‖f‖². The proof is analytic (paper VII, Sections 2–4) plus this external interval-arithmetic certificate
(Arb/FLINT). It is not a formal (Lean) proof.

The directory layout and all files are byte-identical to the author's working tree, so the hashes pinned inside the
certificate remain valid. (The folder names record the internal task numbers of the work; some documentation inside
the work notes is not included, the paper is the reference.)

## Contents
| path | role |
|---|---|
| `python/weil_spectral.py` | Legendre basis, Gauss–Legendre nodes (used by the generator) |
| `notes/0784_gpt_A3_mu12/generate.py` (+ `matrix_utils.py`, `diagnostic_support.py`, `prepare.py`, `run_guarded.py`) | generate the source matrices A, B, J (N = 832, 1920 bits, 256 nodes per frequency panel) as Arb balls |
| `notes/0784_gpt_A3_mu12/n832p1920q256_source_manifest.json`, `n832p1920q256_result.json` | manifest of the source matrices (SHA-256) and the search run |
| `notes/0784_gpt_A3_mu12/budget.py`, `rational_budget.py` (+ `.json`) | the projection-error budget e, n, h (ρ = 3) and α + βb² < 2^−303 in exact rationals |
| `notes/0785_gpt_A4_mu12_certificate/sources/` | the six source matrices used by the certificate (copies of the generated ones, pinned in `SOURCE_PINS.json`) |
| `notes/0785_gpt_A4_mu12_certificate/propose.py` | the proposer: rational matrices and congruence transforms (not trusted) |
| `notes/0785_gpt_A4_mu12_certificate/check_certificate.py`, `codec.py`, `bounds.py` | the checker (does not import the proposer, computes no eigenvalues; re-evaluates the prime Schur bound, the tail scalars, the projection and quadrature budgets) |
| `notes/0785_gpt_A4_mu12_certificate/packet/` | the certificate: 18 dyadic rational matrices and `manifest.json` |
| `notes/0785_gpt_A4_mu12_certificate/verification.json`, `verification_1536.json`, `mutation_results.json` | recorded results (1280 and 1536 bits; seven mutations rejected) |
| `notes/0787_gpt_review_mu12_certificate/independent_certificate.py` (+ `.json`, `.out`) | an independent reimplementation of the final check at 1664 bits (shares no code with the checker) |
| `notes/0787_gpt_review_mu12_certificate/source_spotchecks.py`, `scalar_review.py` | independent recomputation of nine source entries; independent check of the scalar bounds |
| `notes/0786_claude_review_0785/spot_check_*.py` | independent mpmath check of A[0][0], B[0][0] (both parities) and J[0][0] (even) |
| `notes/0782_claude_review_0778/check_defect.py`, `check_general_tail.py` | numerical illustration of Theorem N (c ≤ 20) and of the bounds for e₀, e₁ (not part of the proof) |
| `replay_mutations.py` | runs the mutation tests in a disposable copy (added in v1.5.1; see below) |

## Integrity
    shasum -a 256 -c SHA256SUMS      # from this folder; lists exactly the files shipped in this folder (no Python caches)

Run this first: the independent check below rewrites `independent_certificate.json` (restore it with
`git checkout -- notes/0787_gpt_review_mu12_certificate/independent_certificate.json`).

## Reproducing (from this directory; tested with Python 3.12, python-flint 0.9.0, mpmath 1.3.0, numpy)
    python notes/0785_gpt_A4_mu12_certificate/check_certificate.py --bits 1280 --output replay.json   # ~32 s, ~1.3 GB → PASS_EXTERNAL_CERTIFICATE
    python notes/0787_gpt_review_mu12_certificate/independent_certificate.py                         # ~60 s; rewrites independent_certificate.json → PASS_INDEPENDENT_CERTIFICATE
    python replay_mutations.py                                                                       # ~60 s; seven mutations rejected for the intended reasons

`replay_mutations.py` runs the original `test_mutations_v2.py` in a disposable hard-link copy under `replay/` and
compares the rejection reasons with the shipped `mutation_results.json`. (Running `test_mutations_v2.py` directly in
this folder performs the seven tests but then stops, because it refuses to overwrite the shipped results; the script
itself is unchanged, since its hash is recorded in the provenance files.)

Changes in v1.5.1: `SHA256SUMS` no longer lists two Python cache files that are not shipped (in v1.5.0 the
integrity command therefore exited with status 1, while all shipped files matched), and `replay_mutations.py` was
added. No data file, checker or result changed.

Regenerating the source matrices takes about 16 minutes and 2.5 GB:

    python notes/0784_gpt_A3_mu12/generate.py 832 1920 256 regen      # writes regen_* files; refuses to overwrite

and compare the SHA-256 of the six `regen_p*_{A,B,J}.json.gz` with `n832p1920q256_source_manifest.json` (the files
in `notes/0785_gpt_A4_mu12_certificate/sources/` are byte-identical copies). `prepare.py` only documents how
`generate.py` was derived from an earlier script and is not needed.

Licence: code Apache-2.0, data and documentation CC BY 4.0.
