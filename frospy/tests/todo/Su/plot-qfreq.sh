#! /bin/bash
# Plotting Q vs freq
# Copy color.dat to plot
rm *.ps

#PWD1=02s00/060994A/02s00-07s02
PWD2=02s00/060994A/07s02
PWD3=02s00/060994A/02s00

#PWD1=03s00/060994A/03s00-06s05
#PWD2=03s00/060994A/03s00-00s24
#PWD3=03s00/060994A/06s05
#PWD4=03s00/060994A/00s24
#PWD5=03s00/060994A/03s00

declare -a PWD=($PWD1 $PWD2 $PWD3 $PWD4 $PWD5)
counter=${#PWD[@]}

## Getting the minimum and maximum for all files
f0=0; f1=10; Q0=0; Q1=10000
for (( i=0; i<${counter}; i++ )); do
 for file in ${PWD[$i]}/*omega*; do
  f2=$(tr -s ' ' < $file | cut -d' ' -f2 | sort -g |tail -1)
  fmax=$(python -c "print max($f0,$f2)"); f0=$f2
  f2=$(tr -s ' ' < $file | cut -d' ' -f2 | sort -g |head -1)
  fmin=$(python -c "print min($f1,$f2)"); f1=$f2

  Q2=$(tr -s ' ' < $file | cut -d' ' -f3 | sort -g |tail -1)
  Qmax=$(python -c "print max($Q0,$Q2)"); Q0=$Q2
  Q2=$(tr -s ' ' < $file | cut -d' ' -f3 | sort -g |head -1)
  Qmin=$(python -c "print min($Q1,$Q2)"); Q1=$Q2
 done
done

# Setting limits
fmin=$(python -c "print $fmin*0.999")
fmax=$(python -c "print $fmax*1.001")
stepf=$(python -c "print round(($fmax-$fmin)/3,3)")

Qmin=$(python -c "print $Qmin*0.99")
Qmax=$(python -c "print $Qmax*1.01")
stepQ=$(python -c "print round(($Qmax-$Qmin)/4,0)")

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
 for file in ${PWD[$i]}/*omega*; do
  color=$(ls -l | sed -ne "$((c))p" color.dat )
  if (($i == 0))
  then
   psxy $file -R$fmin/$fmax/$Qmin/$Qmax -JX6.5/6.5 -Y0 -Ba$stepf:"f[mHz]":/a$stepQ:"Q"::."Normal mode frequencies":WSne -Sd0.07 -G$color -K -O >> out.ps
  else
   psxy $file -R -J -Sd0.07 -G$color -K -O >> out.ps
  fi
  ((c++))
 done
done

c=1
xstep=6.7
ystep=5.3
for (( i=0; i<${counter}; i++ )); do
 #legend
 color=$(ls -l | sed -ne "$((c))p" color.dat )
 echo -e "0.1 1.0 12 0 0 LM ${PWD[$i]}" | pstext -N -R0/1/0/1 -JX1/1 -Y$ystep -X$xstep -K -O >> out.ps
 echo -e "0.0 1.0 " | psxy -N -R -J -Sd0.07 -G$color -K -O >> out.ps
 ystep=-0.2; xstep=0
 ((c++))
done

echo 0 0  | psxy -R1/2/1/2 -JX1/2 -Sp -O >> out.ps
gmtset PAGE_ORIENTATION  portrait
#ps2pdf out.ps Q-vs-freq.pdf
gv out.ps
