#!/bin/bash
#
# The calling sequence is: 
#
#  overlap.sh left_out left_state right_out right_state [task] [operator]
#
#  {left,right}_out   should specify GAMESS .dat file containing the orbitals
#                     and optionally the matching .out and/or .inp files
#  {left,right}_state may be an integer state index, in which case GAMESS .out
#                     specofoed by {left,right}_out must exist. Alternatively,
#                     it may also give a determinant list in the internal
#                     superdyson format.
#
#  Optional task and operator are as defined in superdyson; see the source.
#
#  All files may be compressed with bzip2.
#
ls="$(echo "$1"|sed -e 's/\.inp.*$//' -e 's/\.out.*$//' -e 's/\.dat.*$//')"
ln="$(echo "$2"|sed -e 's/.bz2$//')"
rs="$(echo "$3"|sed -e 's/\.inp.*$//' -e 's/\.out.*$//' -e 's/\.dat.*$//')"
rn="$(echo "$4"|sed -e 's/.bz2$//')"
ta="${5:-"braket"}"
op="${6:-"none"}"
#
sd_loc="$(dirname "$0")/"
debug=false
#
if $debug ; then
  tmpdir="ps.$$"
  mkdir -p "${tmpdir}"
  echo "Temporary files are kept in ${tmpdir}"
else
  tmpdir="/tmp/ps.$$"
  mkdir -p "${tmpdir}"
  trap "rm -rf ${tmpdir}" 0
fi
#
#  Prepare header for superdyson.
#  Because GAMESS does not copy the $ECP group to the .dat file, we'll try to
#  extract it from the input.
#
function header () {
  local file="$1" ;
  # Extract C1 geometry section from the .dat file
  bzcat --force ${file}.dat* | awk -f ${sd_loc}/extract_dat.awk -
  bzcat --force ${file}.inp* | awk '/^ \$ECP/,/^ \$END/'
  }
#
function orbitals () {
  local file="$1" ;

  header "${file}" ;
  if bzcat --force ${file}.dat* | grep -sq 'ALPHA MOS: Orbitals at the entry to CI calculation' > /dev/null ; then
    bzcat --force ${file}.dat* | awk '/ALPHA MOS: Orbitals at the entry to CI calculation/,/^ \$END/'
  else
    echo "WARNING: CI MOs not found, extracting \$VEC section" > /dev/stderr
    bzcat --force ${file}.dat* | awk '/^ \$VEC/,/^ \$END/'
  fi
  }
#
function determinants () {
  local file="$1" ;
  local state="$2" ;
  # If state is longer than 3 characters, treat it as a determinant list
  if [ "$(echo "$state" | wc -c)" -gt 4 ] ; then
    bzcat --force ${state}*
  else
    bzcat --force ${file}* | ${sd_loc}/det_conv.awk state="${state}"
  fi
# Delete the 1s core orbitals from the comparison - we do not care about those
# bzcat --force ${file}* | ~/superdyson/det_conv.awk state="${state}" | sed -e 's/  2 2 /  0 0 /'
  }
#
orbitals "${ls}" > "${tmpdir}/ls_orbs.dat"
orbitals "${rs}" > "${tmpdir}/rs_orbs.dat"
#
determinants "${ls}.out" "${ln}" > "${tmpdir}/ls_dets.dets" &
determinants "${rs}.out" "${rn}" > "${tmpdir}/rs_dets.dets" &
wait
#
lmos=$(awk '(NR==1){print NF-1;exit}' "${tmpdir}/ls_dets.dets")
rmos=$(awk '(NR==1){print NF-1;exit}' "${tmpdir}/rs_dets.dets")
ldets=$(wc -l < "${tmpdir}/ls_dets.dets")
rdets=$(wc -l < "${tmpdir}/rs_dets.dets")
#
#
cat > "${tmpdir}/sd.inp" <<eoi
 &sd_v2
   verbose         = 0,
   task            = '${ta}',
   braket_operator = '${op}',
   comment         = '${ta}: ${ls} state ${ln} with ${rs} state ${rn}',
   nmobra          = ${lmos},
   nmoket          = ${rmos},
   ndetbra         = ${ldets},
   ndetket         = ${rdets},
   eps_cdet        = 1e-10,
   dont_check_s2   = .true.
 ! dont_check_s2   = .false.
   file_cbra       = '${tmpdir}/ls_orbs.dat'
   file_cket       = '${tmpdir}/rs_orbs.dat'
   file_detbra     = '${tmpdir}/ls_dets.dets'
   file_detket     = '${tmpdir}/rs_dets.dets'
 /
eoi
#
echo "WARNING: Header corresponds to the bra wavefunction"
header "${ls}"
${sd_loc}/sd_v2.x < "${tmpdir}/sd.inp"
