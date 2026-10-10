# Weil positivity on the window of width log 15 — certificate code

Certificate code for: R. Mori, *Weil positivity on the window of width log 15*, Zenodo, 2026, doi:10.5281/zenodo.23281236 (paper VIII).
This folder is part of https://github.com/moriryota/prolate-weil-bounds (release v1.6.0), archived on Zenodo, doi:10.5281/zenodo.23092422 (all versions).
**The matrices (about 1.8 GB) are not in this repository.** They are in https://github.com/moriryota/weil-log15-certificate-data (archived on Zenodo, doi:10.5281/zenodo.23281176); see "Data" below.

Theorem O: for every complex C² function f supported in [−½ log 15, ½ log 15], Weil's quadratic form satisfies
Q(f) ≥ 2^−2904 ‖f‖². The proof is analytic (paper VIII and paper VII, Sections 2–3) plus this external
interval-arithmetic certificate (Arb/FLINT). It is not a formal (Lean) proof.

The directory layout and the files are byte-identical to the author's working tree, so the hashes pinned inside the
certificate remain valid, with one exception: `notes/0803_gpt_A5_mu15_certificate/REVIEW_0802_PINS.json`. In the working
tree it also pins four internal review documents of the analytic budget (not published; their content is in paper VIII).
The shipped copy lists only the five files of `notes/0802_gpt_review_mu15_budget/` that are included here. The checker
reads this list and verifies each listed file, so the check is otherwise unchanged. The original list has SHA-256
8a0da02ab4444e398aa0b65c86f4d5efd7811b49a0d17dc95198d108fab5640e, which is the value recorded in
`FINAL_CHECKER_CODE_PINS.json`. (The folder names record the internal task numbers of the work; most documentation
inside the work notes is not included, the paper is the reference.)

## Contents
| path | role |
|---|---|
| `python/weil_spectral.py` | Legendre basis, Gauss–Legendre nodes (used by the generator) |
| `notes/0792_claude_lowmem_generator/generate_lowmem.py` (+ `matrix_utils.py`, `run_guarded.py`) | generate the source matrices A, B, J one parity at a time (N = 1664, 3328 bits, 224 nodes per frequency panel) as Arb balls |
| `notes/0792_claude_lowmem_generator/runs/*_source_manifest.json`, `*_diagnostic.json` | SHA-256 manifests of the source matrices; midpoint LDL diagnostic (not part of the proof) |
| `notes/0802_gpt_review_mu15_budget/independent_scalars.py` (+ `.json`), `source_audit.json` | independent check of the analytic scalars (Schur bound with the prime 13, tail scalars, projection budget α + βb² < 2^−697, quadrature < 2^−319) and of the source provenance |
| `notes/0803_gpt_A5_mu15_certificate/propose.py` | the proposer: rational matrices and congruence transforms (not trusted) |
| `notes/0803_gpt_A5_mu15_certificate/check_certificate.py`, `codec.py`, `bounds.py` | the checker (does not import the proposer, computes no eigenvalues; re-evaluates the Schur bound, the tail scalars, the projection and quadrature budgets) |
| `notes/0803_gpt_A5_mu15_certificate/packet/manifest.json` | SHA-256 sums of the 18 certificate matrices (data repository) |
| `notes/0803_gpt_A5_mu15_certificate/check1280_final.json`, `check1536_final.json`, `mutation_results.json`, `mutations/*/result.json` | recorded results (1280 and 1536 bits; seven mutations rejected) |
| `notes/0805_sol_review_mu15_certificate/independent_certificate.py` (+ `.json`, `.out`) | an independent reimplementation of the final check at 1664 bits (shares no code with the checker) |
| `notes/0805_sol_review_mu15_certificate/source_spotchecks.py`, `connection_scalars.py`, `compare_results.py` | independent recomputation of ten source entries (including the prime 13), of the connection scalars, and comparison with the recorded results |
| `notes/0804_claude_review_0803/spot_check.py`, `scalars.py` (+ outputs) | independent mpmath recomputation of 18 source entries; scalar checks |

## Data
Clone or download the data repository and copy the matrices into this folder:

    python restore.py /path/to/prolate-weil-bounds/certificate-VIII     # run inside the data repository

`restore.py` copies the 18 certificate matrices, joins the six source matrices from their 48 MiB parts, and checks
every file against the SHA-256 sums in `packet/manifest.json` and `SOURCE_PINS.json` of this folder.

## Integrity
    shasum -a 256 -c SHA256SUMS      # from this folder; lists exactly the code and result files shipped here

Run this first: the independent check below rewrites `independent_certificate.json` (restore it with
`git checkout -- notes/0805_sol_review_mu15_certificate/independent_certificate.json`).

## Reproducing (from this directory, after restoring the data; tested with Python 3.12, python-flint 0.9.0, mpmath 1.3.0, numpy)
    python notes/0803_gpt_A5_mu15_certificate/check_certificate.py --bits 1280 --output replay.json   # ~200 s, ~3.2 GB → PASS_EXTERNAL_CERTIFICATE
    python notes/0805_sol_review_mu15_certificate/independent_certificate.py                         # ~210 s, ~2.9 GB → both parities PASS
    python notes/0805_sol_review_mu15_certificate/compare_results.py                                 # compares with the recorded results

The certified lower bounds of the finite condition are 2.3979369226…e−70 (even) and 3.9458005898…e−66 (odd).

Regenerating the source matrices takes about 7800 s and 5 GB per parity (one parity at a time):

    python notes/0792_claude_lowmem_generator/generate_lowmem.py mu15 1664 3328 224 0 regen --out regen
    python notes/0792_claude_lowmem_generator/generate_lowmem.py mu15 1664 3328 224 1 regen --out regen

This writes `regen/regen_p{0,1}_{A,B,J,...}.json.gz` (it refuses to overwrite an earlier run). Compare them with
`m15n1664p3328q224_p{0,1}_{A,B,J}.json.gz`: entrywise the balls must overlap (for the window log 12 the same generator
reproduced the matrices of paper VII in this sense). Byte-identical output is expected with the same software versions,
but we have not tested it.

This work used Anthropic Claude and OpenAI GPT to derive bounds, implement the computations and review the results,
under the author's direction. The author is responsible for the content.

Licence: code Apache-2.0, data and documentation CC BY 4.0.
