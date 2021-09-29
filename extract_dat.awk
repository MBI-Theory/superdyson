#!/usr/bin/awk -f
#  This script attempts to extract $DATA section from GAMESS .dat file
#  There is (very) limited symmetry handling, to produce a C1 geometry.
#
BEGIN {
  infty     = 1e20 ; 
  datline   = infty ; 
  geostart  = infty ;
  geo_atom  = 0 ;
  geo_basis = 1 ;
  geostate  = geo_atom ;
  symcode   = "" ;
  eps       = 1e-3 ; # Tolerance for geometry comparison (square of the distance)
  # Symmetry replication tables. The syntax of s_list[] yis a bit weird:
  # First index is the symmetry label; second index is the operation index;
  # The remaining two indices are just a Cartesian variable index of the transformation matrix.
  #
  # C1
  symop_count["c1"] = 0 ;
  # D2
  symop_count["d2"] = 3 ;
  s_list["d2",1,1,1] =-1.0 ; s_list["d2",1,1,2] = 0.0 ; s_list["d2",1,1,3] = 0.0 ;
  s_list["d2",1,2,1] = 0.0 ; s_list["d2",1,2,2] =-1.0 ; s_list["d2",1,2,3] = 0.0 ;
  s_list["d2",1,3,1] = 0.0 ; s_list["d2",1,3,2] = 0.0 ; s_list["d2",1,3,3] = 1.0 ;
  #                                                                  
  s_list["d2",2,1,1] = 1.0 ; s_list["d2",2,1,2] = 0.0 ; s_list["d2",2,1,3] = 0.0 ;
  s_list["d2",2,2,1] = 0.0 ; s_list["d2",2,2,2] =-1.0 ; s_list["d2",2,2,3] = 0.0 ;
  s_list["d2",2,3,1] = 0.0 ; s_list["d2",2,3,2] = 0.0 ; s_list["d2",2,3,3] =-1.0 ;
  #                                                                  
  s_list["d2",3,1,1] =-1.0 ; s_list["d2",3,1,2] = 0.0 ; s_list["d2",3,1,3] = 0.0 ;
  s_list["d2",3,2,1] = 0.0 ; s_list["d2",3,2,2] = 1.0 ; s_list["d2",3,2,3] = 0.0 ;
  s_list["d2",3,3,1] = 0.0 ; s_list["d2",3,3,2] = 0.0 ; s_list["d2",3,3,3] =-1.0 ;
  # Cs
  symop_count["cs"] = 1 ;
  s_list["cs",1,1,1] = 1.0 ; s_list["cs",1,1,2] = 0.0 ; s_list["cs",1,1,3] = 0.0 ;
  s_list["cs",1,2,1] = 0.0 ; s_list["cs",1,2,2] = 1.0 ; s_list["cs",1,2,3] = 0.0 ;
  s_list["cs",1,3,1] = 0.0 ; s_list["cs",1,3,2] = 0.0 ; s_list["cs",1,3,3] =-1.0 ;
  # C2v
  symop_count["c2v"] = 3 ;
  s_list["c2v",1,1,1] =-1.0 ; s_list["c2v",1,1,2] = 0.0 ; s_list["c2v",1,1,3] = 0.0 ;
  s_list["c2v",1,2,1] = 0.0 ; s_list["c2v",1,2,2] =-1.0 ; s_list["c2v",1,2,3] = 0.0 ;
  s_list["c2v",1,3,1] = 0.0 ; s_list["c2v",1,3,2] = 0.0 ; s_list["c2v",1,3,3] = 1.0 ;
  #                                                                  
  s_list["c2v",2,1,1] = 1.0 ; s_list["c2v",2,1,2] = 0.0 ; s_list["c2v",2,1,3] = 0.0 ;
  s_list["c2v",2,2,1] = 0.0 ; s_list["c2v",2,2,2] =-1.0 ; s_list["c2v",2,2,3] = 0.0 ;
  s_list["c2v",2,3,1] = 0.0 ; s_list["c2v",2,3,2] = 0.0 ; s_list["c2v",2,3,3] = 1.0 ;
  #                                                                  
  s_list["c2v",3,1,1] =-1.0 ; s_list["c2v",3,1,2] = 0.0 ; s_list["c2v",3,1,3] = 0.0 ;
  s_list["c2v",3,2,1] = 0.0 ; s_list["c2v",3,2,2] = 1.0 ; s_list["c2v",3,2,3] = 0.0 ;
  s_list["c2v",3,3,1] = 0.0 ; s_list["c2v",3,3,2] = 0.0 ; s_list["c2v",3,3,3] = 1.0 ;
  # C2
  symop_count["c2"] = 1 ;
  s_list["c2",1,1,1] =-1.0 ; s_list["c2",1,1,2] = 0.0 ; s_list["c2",1,1,3] = 0.0 ;
  s_list["c2",1,2,1] = 0.0 ; s_list["c2",1,2,2] =-1.0 ; s_list["c2",1,2,3] = 0.0 ;
  s_list["c2",1,3,1] = 0.0 ; s_list["c2",1,3,2] = 0.0 ; s_list["c2",1,3,3] = 1.0 ;
  # C2h
  symop_count["c2h"] = 3 ;
  s_list["c2h",1,1,1] =-1.0 ; s_list["c2h",1,1,2] = 0.0 ; s_list["c2h",1,1,3] = 0.0 ;
  s_list["c2h",1,2,1] = 0.0 ; s_list["c2h",1,2,2] =-1.0 ; s_list["c2h",1,2,3] = 0.0 ;
  s_list["c2h",1,3,1] = 0.0 ; s_list["c2h",1,3,2] = 0.0 ; s_list["c2h",1,3,3] = 1.0 ;
  #                                                                  
  s_list["c2h",2,1,1] = 1.0 ; s_list["c2h",2,1,2] = 0.0 ; s_list["c2h",2,1,3] = 0.0 ;
  s_list["c2h",2,2,1] = 0.0 ; s_list["c2h",2,2,2] = 1.0 ; s_list["c2h",2,2,3] = 0.0 ;
  s_list["c2h",2,3,1] = 0.0 ; s_list["c2h",2,3,2] = 0.0 ; s_list["c2h",2,3,3] =-1.0 ;
  #                                                                  
  s_list["c2h",3,1,1] =-1.0 ; s_list["c2h",3,1,2] = 0.0 ; s_list["c2h",3,1,3] = 0.0 ;
  s_list["c2h",3,2,1] = 0.0 ; s_list["c2h",3,2,2] =-1.0 ; s_list["c2h",3,2,3] = 0.0 ;
  s_list["c2h",3,3,1] = 0.0 ; s_list["c2h",3,3,2] = 0.0 ; s_list["c2h",3,3,3] =-1.0 ;
  # D2h
  symop_count["d2h"] = 7 ;
  s_list["d2h",1,1,1] =-1.0 ; s_list["d2h",1,1,2] = 0.0 ; s_list["d2h",1,1,3] = 0.0 ;
  s_list["d2h",1,2,1] = 0.0 ; s_list["d2h",1,2,2] =-1.0 ; s_list["d2h",1,2,3] = 0.0 ;
  s_list["d2h",1,3,1] = 0.0 ; s_list["d2h",1,3,2] = 0.0 ; s_list["d2h",1,3,3] = 1.0 ;
  #                                                                     
  s_list["d2h",2,1,1] = 1.0 ; s_list["d2h",2,1,2] = 0.0 ; s_list["d2h",2,1,3] = 0.0 ;
  s_list["d2h",2,2,1] = 0.0 ; s_list["d2h",2,2,2] = 1.0 ; s_list["d2h",2,2,3] = 0.0 ;
  s_list["d2h",2,3,1] = 0.0 ; s_list["d2h",2,3,2] = 0.0 ; s_list["d2h",2,3,3] =-1.0 ;
  #                                                                     
  s_list["d2h",3,1,1] =-1.0 ; s_list["d2h",3,1,2] = 0.0 ; s_list["d2h",3,1,3] = 0.0 ;
  s_list["d2h",3,2,1] = 0.0 ; s_list["d2h",3,2,2] =-1.0 ; s_list["d2h",3,2,3] = 0.0 ;
  s_list["d2h",3,3,1] = 0.0 ; s_list["d2h",3,3,2] = 0.0 ; s_list["d2h",3,3,3] =-1.0 ;
  #                                                                     
  s_list["d2h",4,1,1] = 1.0 ; s_list["d2h",4,1,2] = 0.0 ; s_list["d2h",4,1,3] = 0.0 ;
  s_list["d2h",4,2,1] = 0.0 ; s_list["d2h",4,2,2] =-1.0 ; s_list["d2h",4,2,3] = 0.0 ;
  s_list["d2h",4,3,1] = 0.0 ; s_list["d2h",4,3,2] = 0.0 ; s_list["d2h",4,3,3] =-1.0 ;
  #                                                                     
  s_list["d2h",5,1,1] =-1.0 ; s_list["d2h",5,1,2] = 0.0 ; s_list["d2h",5,1,3] = 0.0 ;
  s_list["d2h",5,2,1] = 0.0 ; s_list["d2h",5,2,2] = 1.0 ; s_list["d2h",5,2,3] = 0.0 ;
  s_list["d2h",5,3,1] = 0.0 ; s_list["d2h",5,3,2] = 0.0 ; s_list["d2h",5,3,3] =-1.0 ;
  #                                                                     
  s_list["d2h",6,1,1] = 1.0 ; s_list["d2h",6,1,2] = 0.0 ; s_list["d2h",6,1,3] = 0.0 ;
  s_list["d2h",6,2,1] = 0.0 ; s_list["d2h",6,2,2] =-1.0 ; s_list["d2h",6,2,3] = 0.0 ;
  s_list["d2h",6,3,1] = 0.0 ; s_list["d2h",6,3,2] = 0.0 ; s_list["d2h",6,3,3] = 1.0 ;
  #                                                                     
  s_list["d2h",7,1,1] =-1.0 ; s_list["d2h",7,1,2] = 0.0 ; s_list["d2h",7,1,3] = 0.0 ;
  s_list["d2h",7,2,1] = 0.0 ; s_list["d2h",7,2,2] = 1.0 ; s_list["d2h",7,2,3] = 0.0 ;
  s_list["d2h",7,3,1] = 0.0 ; s_list["d2h",7,3,2] = 0.0 ; s_list["d2h",7,3,3] = 1.0 ;
  #
  }
#
#  Find the start of the $DATA section
#
/^ \$DATA *$/ { 
  datline = NR ; 
  }
(NR<datline) { 
  next ; 
  }
#
# Everything below here appears after the $DATA line
#
(NR>=datline) && (NR<=datline+1) { 
  print ;  # First two lines are simply copied to the output
  next ; 
  }
/^ \$END *$/ {  
  print ;  # Copy the last line, and stop processing
  datline = infty ; 
  exit ;
  }
# Symmetry line: C1 symmetry
(NR==datline+2) && (/^C1 *0 *$/||/^C1 *$/) {
  symcode = "c1" ;
  print ; # Copy the line
  geoline = NR + 1 ; # First atom is on the next line
  next ;
  }
# Symmetry line: something other than C1
(NR==datline+2) {
   symline = $0 ;
   }
(NR==datline+2) && /^DN  *2 *$/ { symcode = "d2" ; }
(NR==datline+2) && /^CS  *0 *$/ { symcode = "cs" ; }
(NR==datline+2) && /^CN  *2 *$/ { symcode = "c2" ; }
(NR==datline+2) && /^CNV  *2 *$/ { symcode = "c2v" ; }
(NR==datline+2) && /^CNH  *2 *$/ { symcode = "c2h" ; }
(NR==datline+2) && /^DNH  *2 *$/ { symcode = "d2h" ; }
# General-case handling for non-C1 symmetry
(NR==datline+2) {
  # Did we recognize the symmetry code?
  if (symcode=="") {
    printf "Can't handle symmetry '%s'\n", symline > "/dev/stderr" ;
    exit (1) ;
    }
  print "C1" ; # Transform it to C1
  geoline = NR + 2 ; # First atom starts two lines down
  next ;
  }
# Skip until the first geometry line
(NR<geoline) {
  next ;
  }
(geostate==geo_atom){
  at_name = $1 ;
  at_znuc = $2 ;
  at_x    = $3 ;
  at_y    = $4 ;
  at_z    = $5 ;
  geostate=geo_basis ;
  basis   = "" ;
  next ;
  }
(geostate==geo_basis)&&/^ *$/ {
  report_atom() ;
  geostate=geo_atom ;
  next ;
  }
(geostate==geo_basis) {
  basis = basis $0 "\n" ;
  }
#
#  Replicate atom as necessary by symmetry
#
function report_atom (nrep,irep) {
  nrep  = symop_count[symcode] ;
  nuniq = 0 ;
  for (irep=1;irep<=nrep;irep++) {
    replicate_atom(irep,nuniq) ;
    }
  report_single_atom(at_x,at_y,at_z) ;
  }
function replicate_atom (irep, rep_x,rep_y,rep_z,iuniq) {
  # Apply symmetry operation to the unique atom
  rep_x = s_list[symcode,irep,1,1]*at_x + s_list[symcode,irep,1,2]*at_y + s_list[symcode,irep,1,3]*at_z ;
  rep_y = s_list[symcode,irep,2,1]*at_x + s_list[symcode,irep,2,2]*at_y + s_list[symcode,irep,2,3]*at_z ;
  rep_z = s_list[symcode,irep,3,1]*at_x + s_list[symcode,irep,3,2]*at_y + s_list[symcode,irep,3,3]*at_z ;
  # The replica must not be the same as the original atom, or any of the prior replicas
  if (same_atom(at_x,at_y,at_z,rep_x,rep_y,rep_z)) return ;
  for (iuniq=1;iuniq<=nuniq;iuniq++) {
    if (same_atom(uniq_x[iuniq],uniq_y[iuniq],uniq_z[iuniq],rep_x,rep_y,rep_z)) return ;
    }
  ++nuniq ;
  uniq_x[nuniq] = rep_x ;
  uniq_y[nuniq] = rep_y ;
  uniq_z[nuniq] = rep_z ;
  report_single_atom(rep_x,rep_y,rep_z) ;
  }
function same_atom(x1,y1,z1,x2,y2,z2, r2) {
  r2  = (x2-x1)**2 ;
  r2 += (y2-y1)**2 ;
  r2 += (z2-z1)**2 ;
  return r2<=eps ;
  }
#
#  Print a single atom
#
function report_single_atom (x,y,z) {
  printf "%-8s %9.6f %16.12f %16.12f %16.12f\n", at_name, at_znuc, x, y, z ;
  printf "%s\n", basis ;
  }
