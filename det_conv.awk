#!/usr/bin/awk -f
#
#  Version 18 June 2011, with checks for common input mistakes in GAMESS runs
#          04 Feb  2015, supressing zero-amplitude determinants altogether
#          27 Jan  2026, add determinantal-CIS support
#
#  Extract determinant list from GAMESS CI run.
#
#  We support the following CITYP:
#    ALDET - Determinantal full CI
#    GEN   - Determinantal general CI
#    ORMAS - Determinantal occupation-restricted CI
#    CIS   - Requires HAMTYP=DETS
#
#  There is only one input parameter:
#
#    state=  The ordinal number of the state in the output.
#
BEGIN {
  infty    = 1e20 ;
  sl       = infty ;
  scis     = infty ;
  cango    = 0 ;
  type     = "bad" ;
  notmcscf = 1 ;  # Orbital indices in ORMAS output are printed differently
  verbose  = 0 ;
  }
#
#  Check for the most common mistakes in GAMESS CI preparation
#  These warnings should no longer be necessary. now that we dump the orbitals
#  right before the CI calculation starts.
#
#/^ SCFTYP=/ && !/^ SCFTYP=NONE / {
#  printf "WARNING: ORBITALS WILL BE RECOMPUTED (%s). THE TOTAL WAVEFUNCTION IS LIKELY WRONG!\n", $1 > "/dev/stderr" ;
#  }
#/^  *TOLZ  *= .* TOLE  *= / && ( ($3>0.0) || ($6>0.0) ) {
#  printf "WARNING: ORBITALS WILL BE MODIFIED BY $GUESS (TOLZ OR TOLE ARE NOT ZERO). THE TOTAL WAVEFUNCTION IS LIKELY WRONG!\n" > "/dev/stderr" ;
#  }
#/^  *SYMDEN= .* PURIFY= / && ( ($2!="F") || ($4!="F") ) {
#  printf "WARNING: ORBITALS WILL BE SYMMETRIZED BY $GUESS (SYMDEN OR PURIFY ARE TRUE). THE TOTAL WAVEFUNCTION IS LIKELY WRONG!\n" > "/dev/stderr" ;
#  }
#/^ INPUT CARD>[^!]* IPURFY=[12]/ {
#  printf "WARNING: IPURFY MAY BE ACTIVE IN $TRANS. THE TOTAL WAVEFUNCTION MAY BE WRONG!\n" > "/dev/stderr" ;
#  }
#
function vp() {
  if (verbose>0) printf "%08d: %s\n", NR, $0 > "/dev/stderr" ;
  }
function abs(x) { return (x<0)?(-x):(x) ; }
#
#  Choose which set of determinats to use.
#  For CI, we want the first (and only) set.
#  For MCSCF we want the second set.
#
/^ *MCSCF CALCULATION *$/                        { vp() ; cango-- ; notmcscf = 0 ; }
/^ *FINAL MCSCF ENERGY IS/                       { vp() ; cango++ ; }
/^ *DIRECT DETERMINANT ORMAS-CI *$/              { vp() ; cango++ ; type = "ormas" ; }
/^ *AMES LABORATORY DETERMINANTAL GENERAL CI *$/ { vp() ; cango++ ; type = "aldet" ; }
/^ *AMES LABORATORY DETERMINANTAL FULL CI *$/    { vp() ; cango++ ; type = "aldet" ; }
/^ *RESULTS FROM DETERMINANT BASED ATOMIC ORBITAL CI-SINGLES *$/ {
                                                   vp() ; cango++ ; type = "detcis" ; }
#
#  Sizes of active spaces.
#
/^ *NUMBER OF CORE ORBITALS  * = /   { vp() ; ncore = $6 ; nmax = ncore + nact ; }
/^ *NUMBER OF ACTIVE ORBITALS  * = / { vp() ; nact  = $6 ; nmax = ncore + nact ; }
#  CIS: Special case
/^ # CORE ORBITALS      = / { vp() ; ncore = $5 ; }
/^ # OCCUPIED ORBITALS  = / { vp() ; nocc  = $5 ; }
/^ # MOLECULAR ORBITALS = / { vp() ; nmo   = $5 ; }
/^ # BASIS FUNCTIONS    = / { vp() ; nbas  = $5 ; }
#
function count_orbs(ind,start,end,  i,count) {
  count = 0 ;
  for (i=start;i<=end;i++) { 
    if (fields[i]==ind) count++ ;
    }
  if (count>1) {
    print "det_conv.awk: GAGA in line ", NR > "/dev/stderr" ;
    exit 1 ;
    }
  return count ;
  }
#
#  Locating the right state - All but CIS
#
(cango>0)&&/^ *STATE .* ENERGY= .* S= .* SZ=/&&($2==state) {
  printf "%08d: %s\n", NR, $0 > "/dev/stderr" ;
  sl = NR + 4 ;
  if (type=="ormas") sl ++ ;
  }
(NR>=sl) && /^ *$/ { vp() ; exit ; }
(NR>=sl) && /DONE WITH DETERMINANT CI COMPUTATION/ { vp() ; exit ; }
#
#  Locating the right state - CIS. Note that CIS numbers excited states,
#  so that the ground state has to be treated separately.
#
(cango>0)&&/^ *EXCITED STATE .* ENERGY= .* S = .* SPACE SYM/&&($3==(state-1)) {
  printf "%08d: %s\n", NR, $0 > "/dev/stderr" ;
  scis = NR + 6 ;
  }
(NR>=scis) && /^ *---+ *$/ { vp() ; exit ; }
#
#  ORMAS dump. 
#  The ORMAS-formatted line requires special parsing; the Fortran
#  format line is almost (but not quite) like this:
#
#  n(a3),' |',m(a3),' |',fx.y
#
#  As the result, the numerical fields can be missing or run into
#  each other.
#
(type=="ormas")&&(NR>=sl){
  #
  # Old parsing code. Only reliable if there are at most 99 active orbitals!
  #
  # for (ib1=1    ;(ib1<=NF)&&($ib1!="|");ib1++) ;
  # for (ib2=ib1+1;(ib2<=NF)&&($ib2!="|");ib2++) ;
  # if (($ib1!="|")||($ib2!="|")) { go = infty ; next ; }
  #
  # Parse the line, placing logical fields into "fields" array.
  # This should work for up to 999 active orbitals
  #
  nfields = 0 ;  # Number of fields so far
  nsep    = 0 ;  # Number of separator symbols so far
  for (p=1;p<=length();) {
    # Is this a separator?
    token = substr($0,p,2) ;
    if (token==" |") { 
      fields[++nfields] = "|" ; ++nsep ; p += 2 ;
           if (nsep == 1) { ib1 = nfields ; }
      else if (nsep == 2) { ib2 = nfields ; break ; } # Terminate field processing after the second separator
      }
    else {
      # Is this an empty field or orbital index?
      token = substr($0,p,3) ;
           if ( token ~ /^ *$/)        { }                                 # Empty; just skip it
      else if ( token ~ /^ *[0-9]+$/ ) { fields[++nfields] = token + 0 ; } # Integer; must be orbital index
      else {
        printf "det_conv.awk: Error parsing line '%s' at column %d\n", $0, p > "/dev/stderr" ;
        exit(1) ;
        }
      p += 3 ;
      }
    }
  # The rest of the line is one field, containing a real number
  fields[++nfields] = substr($0,p) + 0.0 ;
  if (nsep!=2) { go = infty ; next ; }
  #
  for (orb=1;orb<ib2;orb++) {
    if (orb!=ib1) {
      if (fields[orb]+notmcscf*ncore>nmax) {
        printf "det_conv.awk: Increase nmax to at least %d\n", $orb+ncore > "/dev/stderr" ;
        exit 1 ;
        }
      }
    }
  #
  if (abs(fields[ib2+1])<1e-12) next ;
  printf " %16.12f ", fields[ib2+1] ;
  for (orb=1;orb<=ncore;orb++) { printf " 2" ; }
  #
  for (orb=ncore+1;orb<=nmax;orb++) {
    nalpha = count_orbs(orb-notmcscf*ncore,    1,ib1-1) ;
    nbeta  = count_orbs(orb-notmcscf*ncore,ib1+1,ib2-1) ;
    if ((nalpha==0)&&(nbeta==0)) printf "  0" ;
    if ((nalpha==1)&&(nbeta==1)) printf "  2" ;
    if ((nalpha==0)&&(nbeta==1)) printf " -1" ;
    if ((nalpha==1)&&(nbeta==0)) printf " +1" ;
    }
  printf "\n" ;
  }
#
#  ALDET and GENCI dump
#
(type=="aldet")&&(NR>=sl) {
  sa=$1 ; sb=$3 ; wgt=$5 ;
  if (abs(wgt)<1e-12) next ;
  printf " %16.12f ", wgt ;
  for (i=1;i<=ncore;i++) printf " 2" ;
  for (i=1;i<=nact;i++) {
    oa = substr(sa,i,1) + 0 ; ob = substr(sb,i,1) + 0 ;
         if (oa==1 && ob==1) printf "  2" ;
    else if (oa==1 && ob==0) printf " +1" ;
    else if (oa==0 && ob==1) printf " -1" ;
    else if (oa==0 && ob==0) printf "  0" ;
    else {
      printf "Bad spin strings for index %d: %s and %s\n", i, oa, ob > "/dev/stderr" ;
      exit ;
      }
    }
  printf "\n" ;
  }
#
#  CIS dump
#
#  Special case: The reference state, CIS
(type=="detcis")&&(state==1) {
  printf " 1.0 " ;
  for (i=1;i<=ncore;i++) {
    printf " 2" ;
    }
  for (;i<=nocc+ncore;i++) {
    printf "  2" ;
    }
  for (;i<=nmo;i++) {
    printf "  0" ;
    }
  printf "\n" ;
  exit ;
  }
(type=="detcis")&&(NR>=scis) {
  vp() ;
  src = $2 ; dst = $3 ; wgt = $4 ;
  if (abs(wgt)<1e-12) next ;
  sign = (-1)^(nocc+ncore-src) ; # Phase correction - from CIS orbital substitution to alpha//beta
  printf " %16.12f ", sign * wgt ;
  for (i=1;i<=ncore;i++) {
    printf " 2" ;
    }
  for (;i<=nocc+ncore;i++) {
    if (i!=src) { printf "  2" ; }
    else {
      if ($1=="ALPHA") { printf " -1" ; }
      else             { printf " +1" ; }
      }
    }
  for (i=nocc+ncore+1;i<=nmo;i++) {
    if (i!=dst) { printf "  0" ; }
    else {
      if ($1=="ALPHA") { printf " +1" ; }
      else             { printf " -1" ; }
      }
    }
  printf "\n" ;
  }
#{ print $0 > "/dev/stderr" }
/PUNCHED ALPHA MOS: Orbitals at the entry to CI calculation/ { nopunchwarning = 1 ; }
/PUNCHED ALPHA CIS MOS: Orbitals at the entry to CI calculation/ { nocispunchwarning = 1 ; }
END {
  if ( ! nopunchwarning && type!="detcis" ) {
    printf "WARNING: ALL GAMESS VERSIONS AFTER 2009 MAY MESS UP ORBITAL PHASE, REGARDLESS OF THE INPUT OPTIONS\n" > "/dev/stderr" ;
    printf "WARNING: PLEASE USE A HACKED VERSION, WHICH DUMPS THE ACTUAL ORBITALS IN CI TO THE PUNCH FILE\n" > "/dev/stderr" ;
    }
  if ( ! nocispunchwarning && type=="detcis" ) {
    printf "WARNING: STOCK GAMESS EXECUTABLE MESSES UP CIS ORBITAL PHASES. THE RESULTS WILL BE CHAOTIC.\n" > "/dev/stderr" ;
    printf "WARNING: PLEASE USE A HACKED VERSION, WHICH DUMPS THE ACTUAL ORBITALS IN CIS TO THE PUNCH FILE\n" > "/dev/stderr" ;
    }
  }
