#!/usr/bin/awk -f
#
#  A limited parser for the GUGA output from MCQDPT2 code. We are following the
#  CSF coupling scheme in "2026_11_12 - CSF generatort - Correct.nb", which is
#  hopefully correct - unlike our naive projection-operator implementation in 
#  "2026_09_11 - CSF generator.nb"
#
#  The code is expected to run in SCFTYP=NONE RUNTYP=TRANSITN regime; we didn't
#  bother looking at anything else.
#
#  Due to the wild range of possibilities, we have plenty of input parameters,
#  namely:
#
#    mc=i            Is the sequential number of the MCQDPT2 calculation, 
#                    corresponding to the $MCQDm group in the input
#    ty="cas"/"mp2"  Extract CASSCF or MCQDPT2-corrected wavefuntion (which is 
#                    not quite the hypothetical MCQDPT2 wavefunction!)
#    is=n            State index within the calculation
#    ms=ms           Twice the desired M_S value (twice to keep it an integer!)
#                    The MCQDPT2 calculation itself is only done for the ms=s;
#                    for all other cases we'll apply the S+/S- operator to the 
#                    extracted w.f. as needed
#
#  This is a rather complicated affair - rather than simply parsing the output,
#  we need to convert the CSFs to determinants, and then accumulate the wavefunction
#  in the determinantal representation. We therefore keep the code separate from
#  the main parser in det_conv.awk
#
#  Each determinant is represented by an array of numerical amplitude and two 
#  strings of equal length, one for spin alpha, one for spin beta. The strings 
#  contain zeros (no electrons) and ones (one electron).
#
#  A wavefunction is represented by an array of amplitudes, indexed by the alpha
#  and beta strings.
#
#  The representation here will work for relatively small CI expansions, typically
#  seen in MCQDPT2. It will fail miserably for large GUGA-CI calculations; that
#  requires a real implementation in a compiled language.
#
BEGIN {
  LINT = 1 ;              # Helps us find bugs
  eps  = 1e-13 ;          # Equivalence criterion for floats
  dbg  = 0 ;              # Debug output on
  go   = 0 ;              # Current parser state
  abrt = 0 ;              # Do not try to report the wavefunction at the end; we have failed!
  sort = "sort -k 1.2bnr" # Sneaky way of sorting by the absolute value
  }
(NR==1) {
  if (mc=="" || ty=="" || is=="" || ms=="" ) {
    printf "Required parameters are missing. mc=%s, ty=%s, is=%s, ms=%s\n", mc, ty, is, ms > "/dev/stderr" ;
    abrt = 1 ;
    exit(1) ;
    }
  }
function g_el(str,pos) { # Numerical value at position pos of a str
  return substr(str,pos,1) + 0 ;
  }
# One-electron operators
function g_sz_a(inp,pos,out) {  # sz(pos), alpha component
  out[1] = inp[1] *  0.5 * g_el(inp[2],pos) ;
  out[2] = inp[2] ;
  out[3] = inp[3] ;
  }
function g_sz_b(inp,pos,out) {  # sz(pos), beta component
  out[1] = inp[1] * -0.5 * g_el(inp[3],pos) ;
  out[2] = inp[2] ;
  out[3] = inp[3] ;
  }
function g_substr2(s,pos) {  # a more permissive version of substr, to simplify the latter code
  if (pos>length(s)) return "" ;
  return substr(s,pos) ;
  }
function g_substr3(s,pos,len) {  # a more permissive version of substr, to simplify the latter code
  if (len<=0 || pos>length(s)) return "" ;
  return substr(s,pos,len) ;
  }
function g_sp_b(inp,pos,out) {  # s+(pos), beta component (alpha component is zero)
  out[1] = inp[1] * (1-g_el(inp[2],pos)) * g_el(inp[3],pos) ;
  out[2] = g_substr3(inp[2],1,pos-1) "1" g_substr2(inp[2],pos+1) ;
  out[3] = g_substr3(inp[3],1,pos-1) "0" g_substr2(inp[3],pos+1) ;
  }
function g_sm_a(inp,pos,out) {  # s-(pos), alpha component (beta component is zero)
  out[1] = inp[1] * (1-g_el(inp[3],pos)) * g_el(inp[2],pos) ;
  out[2] = g_substr3(inp[2],1,pos-1) "0" g_substr2(inp[2],pos+1) ;
  out[3] = g_substr3(inp[3],1,pos-1) "1" g_substr2(inp[3],pos+1) ;
  }
# Many-electron operators
function g_norb(key) {  # number of orbitals from the combined key (alpha SUBSEP beta)
  return int((length(key)-length(SUBSEP))/2) ;
  }
function g_key_a(key, norb) { # spin-alpha subkey
  return substr(key,1,norb) ;
  }
function g_key_b(key, norb) { # spin-beta subkey
  return substr(key,norb+1+length(SUBSEP),norb) ;
  }
function g_fill_det(wf,key,det, norb) { # Extract one determinant from a wf, given a key
  norb   = g_norb(key) ;
  det[1] = wf[key] ;
  det[2] = g_key_a(key, norb) ;
  det[3] = g_key_b(key, norb) ;
  }
function g_small(x) {
  return (x>-eps) && (x<eps) ;
  }
function g_accu(wf,det) {
  if (! ((det[2],det[3]) in wf)) wf[det[2],det[3]] = 0.0 ; # Avoid lint warning
  wf[det[2],det[3]] += det[1] ;
  if (g_small(wf[det[2],det[3]])) delete wf[det[2],det[3]] ;
  }
function g_ssz(iwf,owf,  key,norb,iorb,inp,out) {  # Sz
  delete owf ;
  for (key in iwf) {
    norb   = g_norb(key) ;
    g_fill_det(iwf,key,inp) ;
    for (iorb=1;iorb<=norb;iorb++) {
      g_sz_a(inp,iorb,out) ; g_accu(owf,out) ;
      g_sz_b(inp,iorb,out) ; g_accu(owf,out) ;
      }
    }
  }
function g_ssp(iwf,owf,  key,norb,iorb,inp,out) {  # S+
  delete owf ;
  for (key in iwf) {
    norb   = g_norb(key) ;
    g_fill_det(iwf,key,inp) ;
    for (iorb=1;iorb<=norb;iorb++) {
      g_sp_b(inp,iorb,out) ; g_accu(owf,out) ;
      }
    }
  }
function g_ssm(iwf,owf,  key,norb,iorb,inp,out) {  # S-
  delete owf ;
  for (key in iwf) {
    norb   = g_norb(key) ;
    g_fill_det(iwf,key,inp) ;
    for (iorb=1;iorb<=norb;iorb++) {
      g_sm_a(inp,iorb,out) ; g_accu(owf,out) ;
      }
    }
  }
function g_copy(iwf,owf, key) { # owf = iwf
  delete owf ;
  for (key in iwf) {
    owf[key] = iwf[key] ;
    }
  }
function g_add(c,iwf,owf,  key,det) {  # owf += c*iwf2
  for (key in iwf) {
    g_fill_det(iwf,key,det) ;
    det[1] *= c ;
    g_accu(owf,det) ;
    }
  }
function g_scale(c,owf,  key) {  # owf *= c
  for (key in owf) {
    owf[key] *= c ;
    }
  }
# And now for the big things
function g_ss2(iwf,owf, tmp,tm2) { # owf <- S^2 iwf
  g_ssm(iwf,tmp) ; g_ssp(tmp,owf) ;      # S_+ S_- iwf
  g_ssz(iwf,tmp) ; g_add(-1.0,tmp,owf) ; # - S_z iwf
  g_ssz(tmp,tm2) ; g_add( 1.0,tm2,owf) ; # + S_z S_z iwf
  }
function g_zs(ts,s,iwf,owf) { # Spin-annihilation operator. ts is the desired spin, s is the spin to cancel
  g_ss2(iwf,owf) ;            # S^2 iwf
  g_add(-s*(s+1.0),iwf,owf) ; # - s(s+1) iwf
  g_scale(1.0/(ts*(ts+1.0)-s*(s+1.0)),owf) ; # The prefactor
  }
function g_norm(wf, key,norm) { # Wavefunction 2-norm
  norm = 0 ;
  for (key in wf) {
    norm += wf[key]^2 ;
    }
  return norm ;
  }
function g_normalize(wf,  norm) {  # Normalize wavefunction to 1.
  norm = g_norm(wf) ;
  if (!g_small(norm)) { g_scale(1.0/sqrt(norm),wf) ; }
  }
function g_maxbadspin(wf,  key,norb,inp,unc,max_unc) { # Figure out the largest possible S value for the wavefunction
  max_unc = 0 ;
  for (key in wf) {
    norb   = g_norb(key) ;
    g_fill_det(wf,key,inp) ;
    unc = 0 ;
    for (pos=1;pos<=norb;pos++) {
      unc += (g_el(inp[2],pos)!=g_el(inp[3],pos))
      }
    if (unc>max_unc) max_unc = unc ;
    }
  return max_unc ; # Keep it an integer; divide by 2 to get the actual spin
  }
function g_xs(ts,iwf,owf,  maxs,is,tmp,tm2) { # Annihilate all spins other than ts from iwf
  g_copy(iwf,tmp) ;
  maxis = g_maxbadspin(iwf) ;        # Twice the spin!
  for (is=maxis%2;is<=maxis;is+=2) {  # Index is twice the spin, to keep it an integer
    if (g_small(is/2-ts)) continue ;
    # printf "Annihilating %f of %f\n", is/2, maxis/2 ;
    g_zs(ts,is/2,tmp,tm2) ;
    g_copy(tm2,tmp) ;
    }
  g_normalize(tmp) ;
  g_copy(tmp,owf) ;
  }
# Proper implementation of the CSFs (which is is fact probably simpler than the projection version :L))
function g_upc(s,m,updn) { # Upcoupled: (S,M)+(1/2,+-1/2) -> (S+1,M+-1/2); the sign of the new electron's M is updn
  return sqrt((s+1.0+updn*m)/(2*s+1.0)) ;
  }
function g_abs(x) { return x>=0?x:-x }
function g_dnc(s,m,updn) { # Downcoupled: (S,M)+(1/2,+-1/2) -> (S-1,M+-1/2)
  if (g_abs(2*m+updn)>(2*s-1)) { return 0 ; }
  return -updn * sqrt((s-updn*m)/(2*s+1.0)) ;
  }
# Process a single position in a single partially coupled determinant
# This routine is a little complex, mostly because awk does not properly support lists.
# p is the number of spins already coupled
# (s,m) are the spin and the projection
# c is the amplitude
# a and b are the alpha and beta strings, partially coupled
# pl, sl, ml, and cl are the corresponding output arrays, indexed by the coupled a,b
function g_oneD(p,s,m,c,a,b,sl,ml,cl, sa,sb,ca,cb) {
  sa = g_el(a,p+1) ; sb = g_el(b,p+1) ;
  if (sa==sb) { # Closed-shell pair, copy the result to the output
    sl[a,b] = s ; ml[a,b] = m ; cl[a,b] = c ;
    }
  else if (sa==1) { # We are upcoupling
    # Add alpha spin
    sa = g_substr3(a,1,p) "1" g_substr2(a,p+2) ; sb = g_substr3(b,1,p) "0" g_substr2(b,p+2) ;
    sl[sa,sb] = s + 0.5 ; ml[sa,sb] = m + 0.5 ; cl[sa,sb] = c * g_upc(s,m,+1) ;
    # Add beta spin
    sa = g_substr3(a,1,p) "0" g_substr2(a,p+2) ; sb = g_substr3(b,1,p) "1" g_substr2(b,p+2) ;
    sl[sa,sb] = s + 0.5 ; ml[sa,sb] = m - 0.5 ; cl[sa,sb] = c * g_upc(s,m,-1) ;
    }
  else { # sb==1, We are downcoupling
    # Add alpha spin
    sa = g_substr3(a,1,p) "1" g_substr2(a,p+2) ; sb = g_substr3(b,1,p) "0" g_substr2(b,p+2) ;
    sl[sa,sb] = s - 0.5 ; ml[sa,sb] = m + 0.5 ; cl[sa,sb] = c * g_dnc(s,m,+1) ;
    # Add beta spin
    sa = g_substr3(a,1,p) "0" g_substr2(a,p+2) ; sb = g_substr3(b,1,p) "1" g_substr2(b,p+2) ;
    sl[sa,sb] = s - 0.5 ; ml[sa,sb] = m - 0.5 ; cl[sa,sb] = c * g_dnc(s,m,-1) ;
    }
  }
# Process a single position in a list of determinants
function g_one(p,si,mi,ci,so,mo,co, key,norb,a,b) {
  delete so ; delete mo ; delete co ;
  for (key in ci) {
    norb = g_norb(key) ;
    a = g_key_a(key,norb) ;
    b = g_key_b(key,norb) ;
    g_oneD(p,si[key],mi[key],ci[key],a,b,so,mo,co) ;
    }
  }
# Debug printout for CSF conversion
function g_csf_report(s,m,c, key,norb,a,b) {
  for (key in c) {
    norb = g_norb(key) ;
    a = g_key_a(key,norb) ;
    b = g_key_b(key,norb) ;
    printf " s= %3.1f m= %4.1f c= %18.14f %s//%s\n", s[key], m[key], c[key], a, b > "/dev/stderr" ;
    }
  }
# Convert a CSF specified by the a/b strings to a list of determinants
function g_g2d(a,b,m,wf, norb,p,si,mi,ci,so,mo,co) {
  norb = length(a) ;
  delete si ; delete mi ; delete ci ;
  si[a,b] = 0 ; mi[a,b] = 0 ; ci[a,b] = 1.0 ;
  if (dbg) { printf "\ninput determinants are:\n" > "/dev/stderr" ; g_csf_report(si,mi,ci) ; }
  for (p=0;p<norb;p++) {
    g_one(p,si,mi,ci,so,mo,co) ;
    if (dbg) { printf "\np = %d determinants are:\n", p + 1 > "/dev/stderr" ; g_csf_report(so,mo,co) ; }
    delete si ; delete mi ; delete ci ;
    for (key in co) {
      if (g_small(co[key])) continue ;
      si[key] = so[key] ; mi[key] = mo[key] ; ci[key] = co[key] ;
      }
    }
  delete wf ;
  for (key in ci) {
    if (mi[key]!=m) continue ;
    wf[key] = ci[key] ;
    }
  }
# Parity factor for going from "textbook" to the alpha//beta notation
function g_parity(a,b,  norb,io,ind,nel,sw,par) {  # alpha and beta occupations, as used in 1-electron operators
  norb = length(a) ;
  nel = 0 ;
  for (io=1;io<=norb;io++) { # Record the textbook order of the occupied orbitals; there's a gap
                             # between alpha and beta indices, but it does not matter.
    if (g_el(a,io)==1) { ind[++nel] = io ; }
    if (g_el(b,io)==1) { ind[++nel] = io + norb ; }
    }
  # Count permulation parity. We are lazy, and Posix awk lacks sort() - so we will use bubble sort.
  par = 1 ; # Parity
  do {
    sw = 0 ; # Will stay zero if no switches
    for (io=2;io<=nel;io++) {
      if (ind[io-1]<ind[io]) continue ;
      sw = ind[io-1] ; ind[io-1] = ind[io] ; ind[io] = sw ;
      par = -par ;
      }
    } while (sw) ;
  return par ;
  }
# Apply parity factors to the wavefunction; only the amplitudes change
function g_wf_parity(wf, key,det,par) {
  for (key in wf) {
    g_fill_det(wf,key,det) ;
    par = g_parity(det[2],det[3]) ;
    if (dbg) {
      printf "%s/%s: C = %.10f Parity = %d\n", det[2], det[3], det[1]*par, par > "/dev/stderr" ; 
      }
    wf[key] *= par ;
    }
  }
# Reporting for debugging
function g_report_wf(wf, key,norb) {
  for (key in wf) {
    norb   = g_norb(key) ;
    printf "%14.10f %s|%s\n", wf[key], g_key_a(key,norb), g_key_b(key,norb) > "/dev/stderr" ;
    }
  }
# Parsing section for the GAMESS MCQDPT2 output
function g_dp() {
  if (dbg) printf "%s: %s\n", NR, $0 > "/dev/stderr" ; 
  }
(go==0)&&/^ NOW ACTUAL MCQDPT COMPUTATIONS WILL START/ { g_dp() ; go++ ; }
(go==1)&&/^ *[0-9]+-[A-Z][A-Z] MCQDPT CALCULATION/ {
  g_dp() ;
  seq = substr($1,1,length($1)-2) + 0 ; # Drop the fancy "-ST", "-ND", etc suffix
  if (seq == mc) { go++ ; }
  }
(go==2)&&/^ SPIN MULTIPLICITY  *= / { g_dp() ; spinX2 = $4 - 1 } # Double it to keep it an integer!
(go==2)&&/^ SPATIAL STATE SYMMETRY  *= / { g_dp() ; spacesym = $5 } # We  don't really need it though
(go==2)&&/^ # OF STATES  *= / { g_dp() ; nstates = $5 }
(go==2)&&/^ # OF FROZEN CORE ORBITALS  *= / { g_dp() ; ncore = $7 }
(go==2)&&/^ # OF DOUBLY OCCUPIED ORBITALS  *= / { g_dp() ; ndocc = $7 }
(go==2)&&/^ # OF ACTIVE ORBITALS  *= / { g_dp() ; nact = $6 }
(go==2)&&(ty=="cas")&&/^ ###  *CAS-CI RESULTS  *###/ { g_dp() ; go++ ; } # CAS-SCF output
(go==2)&&(ty=="mp2")&&/^ ###  *MCQDPT2 RESULTS  *###/ { g_dp() ; go++ ; } # MCQDPT2 output
(go==3)&&/^ STATE # .* ENERGY = /&&($3==is) { # Ready to start parsing the WF!
  g_dp() ; go++ ; csfrec = NR + 4 ; energy = $6 ;
  delete wf ; # Clear the total wavefunction, just in case
  } 
(go==4)&&(NR>=csfrec)&&(NF!=3) { g_dp() ; exit(0) ; } # We are done, do not parse the rest of the input
function g_process_csf(coeff,str, csf,dets,tmp,i,c,a,b,msX2,par) {  # Convert 02+- style string to alpha/beta string
  delete csf ;
  if (dbg) { 
    printf "\nadding CSF %s with weight %f\n", str, coeff > "/dev/stderr" ;
    }
  a = "" ; b = "" ; msX2 = 0 ;
  for(i=1;i<=nact;i++) {
    c = substr(str,i,1) ;
         if (c=="0") { a = a "0" ; b = b "0" }
    else if (c=="2") { a = a "1" ; b = b "1" }
    else if (c=="+") { a = a "1" ; b = b "0" ; msX2++ }
    else if (c=="-") { a = a "0" ; b = b "1" ; msX2-- }
    else {
      printf "Bad occupation symbol '%s' in CSF, aborting!\n", c > "/dev/stderr" ;
      printf "%d: %s\n", NR, $0 > "/dev/stderr" ;
      abrt = 1 ;
      exit(1) ;
      }
    }
  # The projection method gives wrong results!
  # csf[a,b] = 1.0 ;
  # g_xs(spinX2/2.0,csf,dets) ; # This gives us the determinants for S=spinX2*0.5 and Ms=msX2*0.5
  g_g2d(a,b,spinX2/2.0,dets) ; # Convert CSF to a determinant list
  if (dbg) { 
    printf "determinants for Ms = %.1f are:\n", msX2*0.5 > "/dev/stderr" ; 
    g_report_wf(dets) ;
    }
  while(msX2<ms) { g_ssp(dets,tmp) ; g_copy(tmp,dets) ; msX2+=2 }
  while(msX2>ms) { g_ssm(dets,tmp) ; g_copy(tmp,dets) ; msX2-=2 }
  g_normalize(dets) ; # S+/S- change the norm
  if (dbg) { 
    printf "after shifting to Ms = %.1f, determinants are:\n", msX2*0.5 > "/dev/stderr" ; 
    g_report_wf(dets) ;
    printf "\n" > "/dev/stderr" ;
    }
  par = g_parity(a,b) ;      # THIS IS A HORRIBLE HACK - BUT IT SEEMS TO WORK FOR THE
  g_add(par*coeff,dets,wf) ; # MCQDPT2 CODE IN GAMESS. IT MAY NOT WORK FOR ANYTHING ELSE
  if (dbg) {
    printf "norm is currently %.8f\n", g_norm(wf) > "/dev/stderr" ; 
    }
  }
(go==4)&&(NR>=csfrec) {
  g_dp() ;
  g_process_csf($2,$3) ;
  }
# Report wavefunction in the superdyson format.
function g_report_superdyson(wf, key,det,io,ca,cb) {
  for (key in wf) {
    g_fill_det(wf,key,det) ;
    printf " %+20.14f ", det[1] | sort ;
    for (io=1;io<=ncore+ndocc;io++) { printf " 2" | sort }
    for (io=1;io<=nact;io++) {
      ca = substr(det[2],io,1) ; cb = substr(det[3],io,1) ;
           if (ca==1&&cb==1) { printf "  2" | sort ; }
      else if (ca==0&&cb==0) { printf "  0" | sort ; }
      else if (ca==1&&cb==0) { printf " +1" | sort ; }
      else if (ca==0&&cb==1) { printf " -1" | sort ; }
      else {                   printf "  Q" | sort ; }
      }
    printf "\n" | sort ;
    }
  }
END {
  if (abrt) exit ;
  }
END {
  g_wf_parity(wf) ;
  norm = g_norm(wf) ;
  printf "Extracting %s state %d of symmetry %d, multiplicity %d, and Ms = %.1f (block %d, line %d). E = %.8f Norm = %.8f\n", \
         ty, is, spacesym, spinX2 + 1, 0.5*ms, mc, csfrec, energy, norm > "/dev/stderr" ;
  g_report_superdyson(wf) ;
  close(sort) ;
  }
END {
  close("/dev/stderr") ;
  }
