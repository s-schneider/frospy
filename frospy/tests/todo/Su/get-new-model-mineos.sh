#!/bin/bash

if [ $# -ne 3 ]
then
        printf "USAGE: ./get-new-model-mineos.sh location variable scaling_factor\n"
        printf "location: IC, OC, LM, UM, LLT, ULT\n"
        printf "variable: vs, vp, rho, Qm, Qk, R (only IC)\n"
        printf "scaling_factor: in percentage i.e. 0.02 instead of 2%\n"
        printf "changes to R are calculared w/ node levels i.e 1,2... \n"
        exit
fi

loc=$1; var=$2;
#if [[ $var == R ]]; then # changes calculated in R levels
#  fac=$3
#else
#  fac=$(printf '%f\n' "$( printf '%f + %f\n' "1" "$3" | bc )")
#fi
fac=$(printf '%f\n' "$( printf '%f + %f\n' "1" "$3" | bc )")
in=prem_ocean.bak
out=prem_ocean.txt
rm $out

# $1=r, $2=rho, $3=vpv, $4=vsv, $5=qkappa, $6=qshear, $7=vph, $8=vsh, $9=eta
if [[ $loc == IC ]]; then
  echo "IC"
  if [[ $var == vs ]]; then # iso = 1
    awk -v f="$fac" '{if (NR>=4&&NR<=36) $4=$4*f; if (NR>=4&&NR<=36) $8=$8*f; print}' $in > $out
  elif [[ $var == vp ]]; then # iso = 1
    awk -v f="$fac" '{if (NR>=4&&NR<=36) $3=$3*f; if (NR>=4&&NR<=36) $7=$7*f; print}' $in > $out
  elif [[ $var == rho ]]; then
    awk -v f="$fac" '{if (NR>=4&&NR<=36) $2=$2*f; print}' $in > $out
  elif [[ $var == Qm ]]; then
    awk -v f="$fac" '{if (NR>=4&&NR<=36) $6=$6*f; print}' $in > $out
  elif [[ $var == Qk ]]; then
    awk -v f="$fac" '{if (NR>=4&&NR<=36) $5=$5*f; print}' $in > $out
  elif [[ $var == R ]]; then
    cp $in $out
    R=1221500.
    #This only works decreasing the R up to 3% and increasing R up to 5%
    awk -v f="$fac" '{if (NR>=36&&NR<=37) $1=$1*f; print}' $in > $out
#    if (($fac > 0)); then # incresing R
#      count=i++
#      start=37
#      end="i<="$((36+$fac))
#      n=1 # counter to grab correct radius
#      lref=36 # line with variables to be copies to other lines
#      Rnew=$(awk -v f="$fac" -v s="$start" 'NR==s+f {print $1}' $in)
#    else # decreasing R
#      count=i--
#      start=36
#      end="i>"$((36+$fac))
#      n=-1 # counter to grab correct radius
#      lref=37 # line with variables to be copies to other lines
#      Rnew=$(awk -v f="$fac" -v s="$start" 'NR==s+f {print $1}' $in)
#    fi
#
#    fnew=$(printf '%f\n' "$( printf '%f / %f\n' "$Rnew" "$R"| bc -l)")
#    nic=$((33+$fac))
#    echo "Levels="$fac "factor="$fnew"%" "New R="$Rnew "Original R="$R "nic="$nic
#    echo $count $start $end $n $nic
#    for (( i=$start; $end; $count )); do
#      r=$(awk -v i="$i" -v n="$n" 'NR==i+n {print $1}' $in)
#      l=$(awk -v r="$r" -v lref="$lref" 'NR==lref {$1=r; print}' $in)
#      echo $i $l
#      sed -i.bak "${i}s/.*/$l/g" $out && rm $out.bak # mac os
#      #sed -i "${n}s/.*/$l/g" $out # linux
#    done
#    sed -i.bak "3s/33/$nic/g" $out && rm $out.bak # mac os
#    #sed -i "${n}s/.*/$l/g" $out # linux
  fi
elif [[ $loc == OC ]]; then
  echo "OC"
  if [[ $var == vs ]]; then # iso = 1
    awk -v f="$fac" '{if (NR>=37&&NR<=69) $4=$4*f; if (NR>=37&&NR<=69) $8=$8*f; print}' $in > $out
  elif [[ $var == vp ]]; then # iso = 1
    awk -v f="$fac" '{if (NR>=37&&NR<=69) $3=$3*f; if (NR>=37&&NR<=69) $7=$7*f; print}' $in > $out
  elif [[ $var == rho ]]; then
    awk -v f="$fac" '{if (NR>=37&&NR<=69) $2=$2*f; print}' $in > $out
  elif [[ $var == Qm ]]; then
    awk -v f="$fac" '{if (NR>=37&&NR<=69) $6=$6*f; print}' $in > $out
  elif [[ $var == Qk ]]; then
    awk -v f="$fac" '{if (NR>=37&&NR<=69) $5=$5*f; print}' $in > $out
  elif [[ $var == R ]]; then
    echo "R changes need to be implemented"
  fi
elif [[ $loc == LM ]]; then
  echo "LM"
  if [[ $var == vs ]]; then # iso = 1
    awk -v f="$fac" '{if (NR>=70&&NR<=133) $4=$4*f; if (NR>=70&&NR<=133) $8=$8*f; print}' $in > $out
  elif [[ $var == vp ]]; then # iso = 1
    awk -v f="$fac" '{if (NR>=70&&NR<=133) $3=$3*f; if (NR>=70&&NR<=133) $7=$7*f; print}' $in > $out
  elif [[ $var == rho ]]; then
    awk -v f="$fac" '{if (NR>=70&&NR<=133) $2=$2*f; print}' $in > $out
  elif [[ $var == Qm ]]; then
    awk -v f="$fac" '{if (NR>=70&&NR<=133) $6=$6*f; print}' $in > $out
  elif [[ $var == Qk ]]; then
    awk -v f="$fac" '{if (NR>=70&&NR<=133) $5=$5*f; print}' $in > $out
  elif [[ $var == R ]]; then
    echo "R changes need to be implemented"
  fi
elif [[ $loc == UM ]]; then
  echo "UM"
  if [[ $var == vs ]]; then # iso = 1
    awk -v f="$fac" '{if (NR>=134&&NR<=151) $4=$4*f; if (NR>=134&&NR<=151) $8=$8*f; print}' $in > $out
  elif [[ $var == vp ]]; then # iso = 1
    awk -v f="$fac" '{if (NR>=134&&NR<=151) $3=$3*f; if (NR>=134&&NR<=151) $7=$7*f; print}' $in > $out
  elif [[ $var == rho ]]; then
    awk -v f="$fac" '{if (NR>=134&&NR<=151) $2=$2*f; print}' $in > $out
  elif [[ $var == Qm ]]; then
    awk -v f="$fac" '{if (NR>=134&&NR<=151) $6=$6*f; print}' $in > $out
  elif [[ $var == Qk ]]; then
    awk -v f="$fac" '{if (NR>=134&&NR<=151) $5=$5*f; print}' $in > $out
  elif [[ $var == R ]]; then
    echo "R changes need to be implemented"
  fi
elif [[ $loc == LLT ]]; then
  echo "LLT"
  if [[ $var == vs ]]; then # iso = 1
    awk -v f="$fac" '{if (NR>=152&&NR<=157) $4=$4*f; if (NR>=152&&NR<=157) $8=$8*f; print}' $in > $out
  elif [[ $var == vp ]]; then # iso = 1
    awk -v f="$fac" '{if (NR>=152&&NR<=157) $3=$3*f; if (NR>=152&&NR<=157) $7=$7*f; print}' $in > $out
  elif [[ $var == rho ]]; then
    awk -v f="$fac" '{if (NR>=152&&NR<=157) $2=$2*f; print}' $in > $out
  elif [[ $var == Qm ]]; then
    awk -v f="$fac" '{if (NR>=152&&NR<=157) $6=$6*f; print}' $in > $out
  elif [[ $var == Qk ]]; then
    awk -v f="$fac" '{if (NR>=152&&NR<=157) $5=$5*f; print}' $in > $out
  elif [[ $var == R ]]; then
    echo "R changes need to be implemented"
  fi
elif [[ $loc == ULT ]]; then
  echo "ULT"
  if [[ $var == vs ]]; then # iso = 1
    awk -v f="$fac" '{if (NR>=158&&NR<=184) $4=$4*f; if (NR>=158&&NR<=184) $8=$8*f; print}' $in > $out
  elif [[ $var == vp ]]; then # iso = 1
    awk -v f="$fac" '{if (NR>=158&&NR<=184) $3=$3*f; if (NR>=158&&NR<=184) $7=$7*f; print}' $in > $out
  elif [[ $var == rho ]]; then
    awk -v f="$fac" '{if (NR>=158&&NR<=184) $2=$2*f; print}' $in > $out
  elif [[ $var == Qm ]]; then
    awk -v f="$fac" '{if (NR>=158&&NR<=184) $6=$6*f; print}' $in > $out
  elif [[ $var == Qk ]]; then
    awk -v f="$fac" '{if (NR>=158&&NR<=184) $5=$5*f; print}' $in > $out
  elif [[ $var == R ]]; then
    echo "R changes need to be implemented"
  fi
fi
