#!/bin/bash
out="$1"
ref="$2"
#
if [ ! -r "$out" ] ; then
  echo "$out is not readable"
  exit 1
fi
if [ ! -r "$ref" ] ; then
  echo "$ref is not readable"
  exit 1
fi
#
ref_over="$(awk '/^Total overlap is /{print $4}' "$ref")"
out_over="$(awk '/^Total overlap is /{print $4}' "$out")"
#
ref_dip="$(awk '/^Transition dipole is /{print substr($0,22)}' "$ref")"
out_dip="$(awk '/^Transition dipole is /{print substr($0,22)}' "$out")"
#
ref_oper="$(awk '/^Matrix element of .* operator is:$/{ do { print getline ; } while (NF>=1&&NF<=3) }' "$ref")"
out_oper="$(awk '/^Matrix element of .* operator is:$/{ do { print getline ; } while (NF>=1&&NF<=3) }' "$out")"
#
ref_rdm="$(awk '/^ \$RDMPCE/,/^ \$END/' "$ref")"
out_rdm="$(awk '/^ \$RDMPCE/,/^ \$END/' "$out")"
#
ref_dysnrm="$(awk '/^<psid|psid>/{print $3;exit;}' "$ref")"
out_dysnrm="$(awk '/^<psid|psid>/{print $3;exit;}' "$out")"
#
ref_dyson="$(awk '/^ \$VECDYS/,/^ \$END/' "$ref")"
out_dyson="$(awk '/^ \$VECDYS/,/^ \$END/' "$out")"
#
function chop_dyson () {
  echo "$1" | awk \
'/^[ 0-9.+-EDed]*$/{ 
  printf "%d %d",substr($0,1,2),substr($0,3,3); 
  for (i=0;i<5;i++) printf " %.7f", substr($0,6+15*i,15) ;
  printf "\n" ;
  }' - | awk '{gsub(" [+-]0?.0000000 "," 0.0000000 ");
               gsub(" [+-]0?.0000000$"," 0.0000000");print}'
  }

if [ -n "$ref_dyson" ] ; then
  ref_dyson_chop="$(chop_dyson "$ref_dyson")"
  out_dyson_chop="$(chop_dyson "$out_dyson")"
  if [ "$ref_dyson_chop" = "$out_dyson_chop" ] ; then
    cmp_dyson="ok"
  else
    cmp_dyson="notok"
  fi
else
  cmp_dyson="ok"
fi
#
if [ \( \
       -n "$ref_over"   -o \
       -n "$ref_dip"    -o \
       -n "$ref_oper"   -o \
       -n "$ref_rdm"    -o \
       -n "$ref_dysnrm" -o \
       -n "$ref_dyson"     \
     \) -a \
     "$ref_over"   = "$out_over"   -a \
     "$ref_dip"    = "$out_dip"    -a \
     "$ref_oper"   = "$out_oper"   -a \
     "$ref_rdm"    = "$out_rdm"    -a \
     "$ref_dysnrm" = "$out_dysnrm" -a \
     "$cmp_dyson"  = "ok"             \
   ] ; then
  echo "OK: $out matches $ref"
  exit 0
else
  if [ -z "$ref_over" -a -z "$ref_dip" -a -z "$ref_oper" -a -z "$ref_rdm" -a -z "$ref_dysnrm" ] ; then
    echo "NOTOK: Can't check output - $ref does not contain any useful data"
  else
    echo "NOTOK: $out differs from $ref"
    if [ "$ref_over" != "$out_over" ] ; then
      echo "overlap ref: $ref_over"
      echo "overlap out: $out_over"
    fi
    if [ "$ref_dip" != "$out_dip" ] ; then
      echo "dipole ref: $ref_dip"
      echo "dipole out: $out_dip"
    fi
    if [ "$ref_oper" != "$out_oper" ] ; then
      echo "operator ref: $ref_oper" 
      echo "operator out: $out_oper" 
    fi
    if [ "$ref_rdm" != "$out_rdm"  ] ; then
      echo "1-RDM singulal values ref: $ref_rdm" 
      echo "1-RDM singular values out: $out_rdm" 
    fi
    if [ "$ref_dysnrm" != "$out_dysnrm"  ] ; then
      echo "Dyson norm ref: $ref_dysnrm" 
      echo "Dyson norm out: $out_dysnrm" 
    fi
    if [ "$cmp_dyson" != "ok" ] ; then
      echo "Dyson+property orbital ref:"
      echo "$ref_dyson" 
      echo "Dyson+property orbital out:"
      echo "$out_dyson" 
    fi
  fi
  exit 1
fi
