#! /bin/bash
# Plotting Q vs freq
# Copy color.dat to plot
rm out.ps

#PWD07=fundamental/00s02
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
model=s20rts #s20rts #prem

declare -a modes=( 1s0 2s0 3s0 4s0 5s0 6s0 \
0s3   0s4  0s5  0s6   0s7  0s9  0s13  0s15 0s16  0s19 \
1s10  2s12  5s7  5s8  5s12  5s17 6s10 6s18  7s5  7s7) # 0s2


## Plot
gmtset PAPER_MEDIA=a4
gmtset X_ORIGIN=1i
gmtset Y_ORIGIN=1i
gmtset HEADER_FONT_SIZE=8p
gmtset HEADER_OFFSET=2p
gmtset ANOT_FONT_SIZE=5p
gmtset LABEL_OFFSET=1p
gmtset LABEL_FONT_SIZE=5p
gmtset TICK_LENGTH=0.1c
gmtset ANNOT_OFFSET_PRIMARY=0.25p
gmtset PAGE_ORIENTATION  landscape
#R=R0/7/0/21
R=R0/7/0/10
x=-1.5;y=14.5
j=1
echo 0 0  | psxy -R1/2/1/2 -JX1/2 -Sp -K > out.ps

#for event in 060994A 050306F 061796A 031111B 070508A 081507F 032598B 073095A 092503C 032805D  022710A 021796B 100995C 040107E 111506F 011307A 030994E 100494B 062301E; do
#for event in 060994A 030994E 061796A 070508A 100494B 050306F 031111B 062301E 073095A 032805D 022710A 111506F 011307A; do
for event in 060994A 031111B; do
PWD01=01s00/01s00/$event
PWD02=02s00/02s00/$event
PWD03=03s00/03s00/$event
PWD04=04s00/04s00/$event
PWD05=05s00/05s00/$event
PWD06=06s00/06s00/$event

declare -a PWD=($PWD01 $PWD02 $PWD03 $PWD04 $PWD05 $PWD06 \
$PWD08 $PWD09 $PWD10 $PWD11 $PWD12 $PWD13 $PWD14 $PWD15 $PWD16 $PWD17 \
$PWD18 $PWD19 $PWD20 $PWD21 $PWD22 $PWD23 $PWD24 $PWD25 $PWD26 $PWD27) # $PWD07

counter=${#PWD[@]}
c=1
for (( i=0; i<${counter}; i++ )); do
 for file in ${PWD[$i]}/$event.$model.peaks; do
  if (($i == 0))
  then
   #ignore comments
   awk 'NF && $1!~/^#/ {print($2,$8)}' $file | psxy -R$R -JX5/3 -X$x -Y$y -Ba1:"f[mHz]":/a2.5g2.5:"Source size rel. CMT"::.$event:WSne -Sd0.02 -Gblack -K -O >> out.ps
   mean=$(awk 'NF && $1!~/^#/  {print $8}' $file | python -c "import sys; tot = [float(n) for n in sys.stdin]; print(sum(tot)/len(tot))" )
   min=$(awk 'NF && $1!~/^#/  {print $8}' $file | python -c "import sys; tot = [float(n) for n in sys.stdin]; print(min(tot)*0.9)" )
   f=$(sed -n 5p $file | cut -f2)
   echo $f $mean | psxy  -R -J -Sd0.15 -Gred -K -O >> out.ps

  else
   awk 'NF && $1!~/^#/ {print($2,$8)}' $file | psxy  -R -J -Sd0.02 -Gblack -K -O >> out.ps
   mean=$(awk 'NF && $1!~/^#/  {print $8}' $file | python -c "import sys; tot = [float(n) for n in sys.stdin]; print(sum(tot)/len(tot))")
   min=$(awk 'NF && $1!~/^#/  {print $8}' $file | python -c "import sys; tot = [float(n) for n in sys.stdin]; print(min(tot)*0.9)" )
   f=$(sed -n 5p $file | cut -f2)

   # Fundamental and other modes
   echo $f $mean | psxy  -R -J -Sd0.1 -Gblue -K -O >> out.ps
   # Radial modes in red
   if(($i < 6)); then echo $f $mean | psxy  -R -J -Sd0.15 -Gred -K -O >> out.ps; fi
  fi

  #legend
  #echo -e "$f $min 8 90 0 RB ${modes[$i]}" | pstext -N -R -J -K -O >> out.ps
  #echo -e "$f $max 8 90 0 LM ${modes[$i]}" | pstext -N -R -J -K -O >> out.ps
  #echo $file $f ${modes[$i]} 'mean' $mean
  ((c++))
 done
done

# Skipping line
if (($j % 5 != 0))
then
 x=5.8
 y=0
else
 x=-23.2
 y=-4.1
fi
((j++))

#x=5.8;y=0
done

echo 0 0  | psxy -R1/2/1/2 -JX1/2 -Sp -O >> out.ps
gmtset PAGE_ORIENTATION  portrait
#ps2pdf out.ps Q-vs-freq.pdf
gv out.ps
cp out.ps cmt.$model.ps
