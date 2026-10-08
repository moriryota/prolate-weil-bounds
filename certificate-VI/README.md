# Weil positivity on the window of width log 10 — certificate data and code

Certificate data and code for: R. Mori, *Weil positivity on the window of width log 10: an interval-arithmetic certificate*, Zenodo, 2026, doi:10.5281/zenodo.23238480 (paper VI).
This folder is part of https://github.com/moriryota/prolate-weil-bounds (release v1.4.0), archived on Zenodo, doi:10.5281/zenodo.23092422 (all versions). The data files are included.

Theorem: for every complex C² function f supported in [−½ log 10, ½ log 10], Weil's quadratic form satisfies
Q(f) ≥ 2^−49162 ‖f‖². The proof is analytic plus this external interval-arithmetic certificate (Arb/FLINT).
It is not a formal (Lean) proof.

The directory layout and all files are byte-identical to the author's working tree, so the hashes pinned inside the
certificate remain valid. (The folder names record the internal task numbers of the work.)

## Contents
| path | role |
|---|---|
| `python/weil_spectral.py` | Legendre basis, Gauss–Legendre nodes (used by the generator) |
| `notes/0775_gpt_A0_N640/generate640.py`, `diagnostic640.py` | generate the source matrices A, B, J (N = 640, 1280 bits, 192 nodes per frequency panel) as Arb balls |
| `notes/0775_gpt_A0_N640/n640p1280q192_*` | the saved source matrices and their manifest (pinned by SHA-256) |
| `notes/0775_gpt_A0_N640/verify_exact_parameters.py` | exact rational check of the projection errors e, n, h, α, β and n < m |
| `notes/0770_gpt_A0_mu10/liu_route_scalars_v2.py` (+ `.json`) | interval check of ‖C_p‖ ≤ 3 (1792 cells) and the scalar inputs of the tail lemma |
| `notes/0776_gpt_A1_certificate/propose.py` | the proposer: rational matrices and congruence transforms (not trusted) |
| `notes/0776_gpt_A1_certificate/check_certificate.py`, `codec.py`, `bounds.py` | the checker (does not import the proposer, computes no eigenvalues) |
| `notes/0776_gpt_A1_certificate/packet/` | the certificate: 18 dyadic rational matrices and `manifest.json` |
| `notes/0776_gpt_A1_certificate/SOURCE_PINS.json`, `MATRIX_SHA256.json`, `verification.json` | pins and the recorded verification result |
| `notes/0776_gpt_A1_certificate/test_mutations_v2.py`, `scalar_checks.py` | five mutation tests; scalar checks |
| `notes/0777_gpt_review_A1/independent_certificate.py` (+ `SOURCE_PINS.json`, `independent_results.json`) | an independent reimplementation of the final check at 1536 bits; it reads its own copy of the packet (see below) |
| `notes/0777_gpt_review_A1/analytic_checks.py` | rational checks of the analytic lemmas (δ, τ, e₀, e₁, digamma) |

## Integrity
    shasum -a 256 -c SHA256SUMS      # from this folder

## Reproducing (from this directory; tested with Python 3.12.11, python-flint 0.9.0, mpmath 1.3.0, numpy 2.5.3)
    python notes/0776_gpt_A1_certificate/check_certificate.py --output replay.json      # ~20 s, ~1 GB → PASS_EXTERNAL_CERTIFICATE
    cp -R notes/0776_gpt_A1_certificate/packet notes/0777_gpt_review_A1/packet          # the reviewer's copy is byte-identical
    python notes/0777_gpt_review_A1/independent_certificate.py                          # overwrites independent_results.json → PASS
    python notes/0776_gpt_A1_certificate/test_mutations_v2.py                           # 5 mutations rejected for the intended reasons
Regenerating the source matrices (`generate640.py`) takes about 6 minutes and writes new files; compare their SHA-256
with `n640p1280q192_source_manifest.json`.

Licence: code Apache-2.0, data and documentation CC BY 4.0.
