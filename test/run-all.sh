#!/bin/bash
bin=../sd_v2.x
ok=0
notok=0
dunno=0
notoklist=""
dunnolist=""
#
function test () {
  for inp in "$@" ; do
    out="$(echo "$inp" | sed -e 's/\.inp$/.out/')"
    chk="$(echo "$inp" | sed -e 's/\.inp$/.chk/')"
    ref="$(echo "$inp" | sed -e 's/\.inp$/.out_ref/')"
    echo "============================================================"
    if [ -r "$out" ] ; then
      echo "Output file $out already exists. Skipping the run"
    else
      echo "Executing $inp"
      $bin < "$inp" > "$out" 2>&1
    fi
    if [ -x "./$chk" -a -r "$ref" ] ; then
      ./"$chk" "$out" "$ref"
      if [ $? -eq 0 ] ; then
        ok=$(($ok+1))
      else
        notok=$(($notok+1))
        notoklist="$notoklist $out"
      fi
    else
      dunno=$(($dunno+1))
      dunnolist="$dunnolist $out"
      [ ! -x "./$chk" ] && echo "WARNING: No check script $chk, skipping verification of the results"
      [ ! -r   "$ref" ] && echo "WARNING: No reference output $ref, skipping verification of the results"
    fi
  done
  echo "============================================================"
  echo "$ok tests gave expected results"
  echo "$notok tests failed"
  echo "$dunno runs can't be checked"
  #
  [ -n "$notoklist" ] && echo "Failing tests are: $notoklist"
  [ -n "$dunnolist" ] && echo "Undecided tests are: $dunnolist"
  [ -n "$notoklist" ] && status=1
  }
#
cheap="h2a-slater_dysonx_sym.inp h2b-slater_dysonx_sym.inp h2p-slater_srdm.inp \
       c6h5f-lowdin_rdm_nosym.inp c6h5f-lowdin_rdm_sym.inp co2-slater_dyson_nosym.inp co2-slater_dyson_sym.inp \
       co-rhf_dipole_sym.inp co-rhf_en_C_sym.inp co-rhf_en_O_sym.inp co2-slater_dysonx_sym.inp \
       sf6c-slater_rdm_sym.inp so2-slater_dipole_nosym.inp so2-slater_dipole_sym.inp \
       so2-slater_overlap_nosym.inp so2-slater_overlap_sym.inp ar2p-slater_srdm_nosym.inp"
medium="co2-lowdin_dysonx_nosym.inp co2-lowdin_dysonx_sym.inp co2-slater_dysonx_nosym.inp n2-lowdin_prop_sym.inp \
       n2-lowdin_rdm_sym.inp sf6c-slater_rdm_nosym.inp so2-lowdin_dipole_sym.inp"
long="n2-lowdin_prop_nosym.inp n2-lowdin_rdm_nosym.inp so2-lowdin_dipole_nosym.inp"
#
status=0
echo "Running cheap tests"
test $cheap
[ "$1" == "cheap" ] && exit $status
#
echo "Running moderately-expensive tests"
test $medium
[ "$1" == "medium" ] && exit $status
#
echo "Running expensive tests"
test $long
exit $status
