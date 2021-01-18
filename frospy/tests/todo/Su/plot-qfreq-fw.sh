#! /bin/bash
# Plotting Q vs freq
# Copy color.dat to plot
rm *.ps

echo "Type modename folder, followed by [ENTER]:"
read modename

#modename=0-3mHz # 01s05-02s04 01s06-02s05 02s10-04s05 03s07-05s05 02s12 01s10 05s03
#PWD1=$modename-self/dat-matrix
PWD2=$modename-cross/dat-matrix
#PWD3=$modename-msrd/dat-matrix

declare -a PWD=($PWD2) # $PWD2) #$PWD3)
counter=${#PWD[@]}
omega=omega

echo "Type fmin of the freqency window, followed by [ENTER]:"
read fmin

echo "Type fmax of the frequency window, followed by [ENTER]:"
read fmax

## Getting the minimum and maximum for all files
#f0=0; f1=10
Q0=0; Q1=10000
for (( i=0; i<${counter}; i++ )); do
 for file in ${PWD[$i]}/*$omega*; do
  Q2=$(tr -s ' ' < $file | awk -v fmin="$fmin" -v fmax="$fmax" '(NR>1) && ($1 > fmin) && ($1 < fmax)' | cut -d' ' -f3 | sort -g | tail -1)
  Qmax=$(python -c "print max($Q0,$Q2)"); Q0=$Q2
  Q2=$(tr -s ' ' < $file | awk -v fmin="$fmin" -v fmax="$fmax" '(NR>1) && ($1 > fmin) && ($1 < fmax)' | cut -d' ' -f3 | sort -g | head -1)
  Qmin=$(python -c "print min($Q1,$Q2)"); Q1=$Q2
 done
done

# Setting limits
fmin=$(python -c "print $fmin*0.999")
fmax=$(python -c "print $fmax*1.001")
stepf=$(python -c "print round(($fmax-$fmin)/3,3)")

#Qmin=$(python -c "print $Qmin*0.99") #Qmin=0
#Qmax=$(python -c "print $Qmax*1.01") #Qmax=600
Qmin=234.11
Qmax=543.11
stepQ=$(python -c "print round(($Qmax-$Qmin)/4,0)")
x0=$(python -c "print $fmax*1.001")
x1=$(python -c "print $x0*0.9999")
y0=$Qmax

echo 'fmin='$fmin ' fmax='$fmax ' fstep='$stepf
echo 'Qmin='$Qmin ' Qmax='$Qmax ' Qstep='$stepQ

## Plot
gmtset PAPER_MEDIA=a4
gmtset X_ORIGIN=1i
gmtset Y_ORIGIN=1i
gmtset HEADER_FONT_SIZE=15p
gmtset HEADER_OFFSET=8p
gmtset ANOT_FONT_SIZE=15p
gmtset LABEL_FONT_SIZE=15p
gmtset PAGE_ORIENTATION  landscape

echo 0 0  | psxy -R1/2/1/2 -JX1/2 -Sp -K > out.ps
c=1
for (( i=0; i<${counter}; i++ )); do
 for file in ${PWD[$i]}/*$omega*; do
  color=$(ls -l | sed -ne "$((c))p" color.dat )
  if (($i == 0))
  then
   psxy $file -R$fmin/$fmax/$Qmin/$Qmax -JX6.5/6.5 -Y0 -Ba$stepf:"f[mHz]":/a$stepQ:"Q"::."Normal mode frequencies":WSne -Sd0.07 -G$color -K -O >> out.ps
  else
   psxy $file -R -J -Sd0.07 -G$color -K -O >> out.ps
  fi

  ## Legend
  y0=$(python -c "print $y0*0.975")
  legend1=$(echo $file | cut -d'R' -f2- | rev | cut -d'.' -f2- | rev)
  legend2=$(echo $file | cut -d'/' -f1 | rev| cut -d'-' -f1 | rev)
  echo -e "$x0 $y0 12 0 0 LM R=$legend1, $legend2" | pstext -N -R -J -K -O >> out.ps
  echo -e "$x1 $y0" | psxy -N -R -J -Sd0.07 -G$color -K -O >> out.ps
  ((c++))
 done
done

# Mode name
echo -e "$x0 $Qmax 13 0 0 LM $modename" | pstext -N -R -J -K -O >> out.ps

echo 0 0  | psxy -R1/2/1/2 -JX1/2 -Sp -O >> out.ps
gmtset PAGE_ORIENTATION  portrait
#ps2pdf out.ps Q-vs-freq.pdf
gv out.ps
