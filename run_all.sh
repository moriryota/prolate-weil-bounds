#!/usr/bin/env bash
# Run every script from the repository root and compare its output with outputs/.
# Usage: ./run_all.sh [proofs|numerics|all]   (default: proofs; "all" includes numerics that take hours)
set -u
PY=${PYTHON:-python3}
MODE=${1:-proofs}
NEW=run_all_outputs
mkdir -p $NEW/proofs-I $NEW/proofs-II $NEW/numerics-I $NEW/numerics-II
status=0
run() {  # run DIR SCRIPT [ARGS...]  -> NEW/DIR/SCRIPT[_ARGS].out, compared with outputs/DIR/same name
  d=$1; s=$2; shift 2
  tag=$s; [ $# -gt 0 ] && tag="${s}_$(echo "$@" | tr ' ' '_')"
  echo "== $d/$s.py $*"
  "$PY" $d/$s.py "$@" > $NEW/$d/$tag.out 2>&1 || { echo "   FAILED (exit $?)"; status=1; }
  if [ -f outputs/$d/$tag.out ]; then
    if diff -q $NEW/$d/$tag.out outputs/$d/$tag.out > /dev/null; then echo "   identical to outputs/$d/$tag.out"
    else echo "   DIFFERS from outputs/$d/$tag.out"; status=1; fi
  fi
}
# ---- proofs (paper I)
for s in chi_bound constants_proof bessel_compact bl1_independent assembly_constant algebra_checks certify assemble; do run proofs-I $s; done
# ---- proofs (paper II)
for s in chi_bound lg_delta_bound L5b_constants S1e_selfcheck theoremC_constants kappa_inf_enclosure L6_effective_constants theoremE_constants; do run proofs-II $s; done
if [ "$MODE" = all ] || [ "$MODE" = numerics ]; then
  for s in theorem2_ratio_check slepian_identity_check liouville_delta_scan bessel_constants_scan BL3_spot_check; do run numerics-I $s; done
  for s in kappa_limit hermite_residual theoremC_sanity kvector_endpoint kvector_outofband L6_checks weighted_error_check; do run numerics-II $s; done
  run numerics-II lemma6_check 5
  run numerics-II galerkin_W 5 7 11 13 15 17
fi
echo; [ $status = 0 ] && echo "All runs completed and all compared outputs are identical." || echo "Some runs failed or differ (see above)."
exit $status
