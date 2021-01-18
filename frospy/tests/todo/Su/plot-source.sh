#! /bin/bash
# Plotting Q vs freq
# Copy color.dat to plot
rm *.ps

#event=060994A # Bolivia Mw 8.3
#event=100494B # Kuril Mw 8.3
#event=050306F
#event=061796A # Flores MW 7.8
event=031111B # Tohoku Mw 9.1
#event=070508A
#event=081507F
#event=081676B
#event=062301E # Peru Mw 8.4
#event=032598B
#event=073095A
#event=092503C # Hokkaido Mw 8.3
#event=032805D # Sumatra Mw 8.6
#event=081977B
#event=022710A # Chile Mw 8.8
#event=121279A
#event=021796B # Iran Mw 8.2
#event=100995C
#event=040107E
#event=111506F # Kuril Mw 8.3
#event=011307A
#event=030385A
#event=062277A
#event=120678A
#event=030994E

PWD01=01s00/01s00/$event
PWD02=02s00/02s00/$event
PWD03=03s00/03s00/$event
PWD04=04s00/04s00/$event
PWD05=05s00/05s00/$event
PWD06=06s00/06s00/$event
PWD07=fundamental/00s03 # should be 00s02
PWD08=fundamental/00s03
PWD09=fundamental/00s04
PWD10=fundamental/00s05
PWD11=fundamental/00s06
PWD12=fundamental/00s07
PWD13=fundamental/00s09
PWD14=fundamental/00s13
PWD15=fundamental/00s15
PWD16=fundamental/00s16
PWD17=fundamental/00s19
PWD18=fundamental/01s10
PWD19=fundamental/02s12
PWD20=fundamental/05s07
PWD21=fundamental/05s08
PWD22=fundamental/05s12
PWD23=fundamental/05s17
PWD24=fundamental/06s10
PWD25=fundamental/06s18
PWD26=fundamental/07s05
PWD27=fundamental/07s07
model=s20rts #prem

declare -a modes=( 1s0 2s0 3s0 4s0 5s0 6s0 \
0s3   0s3   0s4  0s5  0s6   0s7  0s9  0s13  0s15 0s16  0s19 \
1s10  2s12  5s7  5s8  5s12  5s17 6s10 6s18  7s5  7s7)

declare -a PWD=($PWD01 $PWD02 $PWD03 $PWD04 $PWD05 $PWD06 \
$PWD07 $PWD08 $PWD09 $PWD10 $PWD11 $PWD12 $PWD13 $PWD14 $PWD15 $PWD16 $PWD17 \
$PWD18 $PWD19 $PWD20 $PWD21 $PWD22 $PWD23 $PWD24 $PWD25 $PWD26 $PWD27)

counter=${#PWD[@]}

## Plot
gmtset PAPER_MEDIA=a4
gmtset X_ORIGIN=1i
gmtset Y_ORIGIN=1i
gmtset HEADER_FONT_SIZE=15p
gmtset HEADER_OFFSET=3p
gmtset ANOT_FONT_SIZE=15p
gmtset LABEL_FONT_SIZE=15p
gmtset PAGE_ORIENTATION  landscape
R=R0/7/-0.5/100

echo 0 0  | psxy -R1/2/1/2 -JX1/2 -Sp -K > out.ps
c=1
for (( i=0; i<${counter}; i++ )); do
 for file in ${PWD[$i]}/$event.$model.peaks; do
  if (($i == 0))
  then
   #ignore comments
   #awk 'NF && $1!~/^#/ {print($2,$8)}' $file | psxy -R$R -JX25/15 -Y0 -Ba0.5:"f[mHz]":/a5g10-0u:"Source size rel. CMT"::.$event:WSne -Sd0.07 -Gblack -K -O >> out.ps
   awk 'NF && $1!~/^#/ {print($2,$8)}' $file | psxy -R$R -JX25/15 -Y0 -Ba0.5:"f[mHz]":/a0.5g3-0.5u:"Source size rel. CMT"::.$event:WSne -Sd0.07 -Gblack -K -O >> out.ps
   mean=$(awk 'NF && $1!~/^#/  {print $8}' $file | python -c "import sys; tot = [float(n) for n in sys.stdin]; print(sum(tot)/len(tot))" )
   #mean=$(awk 'NF && $1!~/^#/  {print $8}' $file | python -c "import sys; tot = [float(n) for n in sys.stdin]; tot=[e for e in tot if e < 5.5]; print(sum(tot)/len(tot))")
   #max=$(awk 'NF && $1!~/^#/  {print $8}' $file | python -c "import sys; tot = [float(n) for n in sys.stdin]; print(max(tot))" )
   min=$(awk 'NF && $1!~/^#/  {print $8}' $file | python -c "import sys; tot = [float(n) for n in sys.stdin]; print(min(tot)*0.9)" )
   f=$(sed -n 5p $file | cut -f2)
   echo $f $mean | psxy  -R -J -Sd0.25 -Gred -K -O >> out.ps

  else
   awk 'NF && $1!~/^#/ {print($2,$8)}' $file | psxy  -R -J -Sd0.07 -Gblack -K -O >> out.ps
   mean=$(awk 'NF && $1!~/^#/  {print $8}' $file | python -c "import sys; tot = [float(n) for n in sys.stdin]; print(sum(tot)/len(tot))")
   #mean=$(awk 'NF && $1!~/^#/  {print $8}' $file | python -c "import sys; tot = [float(n) for n in sys.stdin]; tot=[e for e in tot if e < 5.5]; print(sum(tot)/len(tot))")
   max=$(awk 'NF && $1!~/^#/  {print $8}' $file | python -c "import sys; tot = [float(n) for n in sys.stdin]; print(max(tot))" )
   min=$(awk 'NF && $1!~/^#/  {print $8}' $file | python -c "import sys; tot = [float(n) for n in sys.stdin]; print(min(tot)*0.9)" )
   f=$(sed -n 5p $file | cut -f2)

   # Fundamental and other modes
   echo $f $mean | psxy  -R -J -Sd0.25 -Gblue -K -O >> out.ps
   # Radial modes in blue
   if(($i <= 6)); then echo $f $mean | psxy  -R -J -Sd0.25 -Gred -K -O >> out.ps; fi
  fi

  #legend
  echo -e "$f $min 8 90 0 RB ${modes[$i]}" | pstext -N -R -J -K -O >> out.ps
  #echo -e "$f $max 8 90 0 LM ${modes[$i]}" | pstext -N -R -J -K -O >> out.ps
  echo $file $f ${modes[$i]} 'mean' $mean
  ((c++))
 done
done


echo 0 0  | psxy -R1/2/1/2 -JX1/2 -Sp -O >> out.ps
gmtset PAGE_ORIENTATION  portrait
#ps2pdf out.ps Q-vs-freq.pdf
gv out.ps
