#/bin/bash
# Plotting the spectrum of one station channels
rm *.grd
rm *.ps

echo "Type mode folder name, followed by [ENTER]:"
read modename

#modename=0-3mHz # 01s05-02s04 01s06-02s05 02s10-04s05 03s07-05s05 02s12 01s10 05s03
#PWD1=$modename-self/060994A/R-0.25
#PWD2=$modename-self/060994A/R0
#PWD3=$modename-self/060994A/R0.25
PWD4=$modename-cross/060994A/R-0.25
PWD5=$modename-cross/060994A/R0.25
PWD6=$modename-cross/060994A/R0
#PWD7=$modename-msrd/060994A/R-0.25
#PWD8=$modename-msrd/060994A/R0
#PWD9=$modename-msrd/060994A/R0.25
#modename=$(echo $PWD1 | cut -d'-' -f-2)

declare -a PWD=($PWD4 $PWD5 $PWD6) #$PWD1 $PWD2 $PWD3 $PWD4 $PWD5 $PWD6 $PWD7 $PWD8 $PWD9)
counter=${#PWD[@]}
ls ${PWD[0]}/*.pf | rev |cut -d'.' -f2 | rev | cut -d'_' -f1 > stations.dat

echo "Type fmin of the freqency window, followed by [ENTER]:"
read f0

echo "Type fmax of the frequency window, followed by [ENTER]:"
read f1

while read staname; do

#echo "Type the station name, followed by [ENTER]:"
#read staname

#echo "Type the station channel (N/E/Z/R/T), followed by [ENTER]:"
#read channel
channel=Z

gmtset PAPER_MEDIA=a4
gmtset X_ORIGIN=1i
gmtset Y_ORIGIN=5i
gmtset HEADER_FONT_SIZE=20p
gmtset HEADER_OFFSET=0p
gmtset ANOT_FONT_SIZE=15p
gmtset LABEL_FONT_SIZE=15p

gmtset PAGE_ORIENTATION  landscape
echo 0 0  | psxy -R1/2/1/2 -JX1/2 -Sp -K > out.ps

mfre=$(sed -n '2p' ${PWD[0]}/*${staname}*${channel}.window | cut -f1); mfre=$(python -c "print round($mfre,3)")
mfim=$(sed -n '2p' ${PWD[0]}/*${staname}*${channel}.window | cut -f2); mfim=$(python -c "print round($mfim,3)")

#f0=$(sed -n '3p' ${PWD[0]}/*${staname}*${channel}.window | cut -f1)
#f1=$(sed -n '3p' ${PWD[0]}/*${staname}*${channel}.window | cut -f2)

h0=$(sed -n '4p' ${PWD[0]}/*${staname}*${channel}.window | cut -f1)
h1=$(sed -n '4p' ${PWD[0]}/*${staname}*${channel}.window | cut -f2)

#### Plotting phase
a0=$(cut -f2 ${PWD[0]}/*${staname}*${channel}.pf | sort -g | tail -1)
a1=$(cut -f2 ${PWD[0]}/*${staname}*${channel}.pfsyn | sort -g | tail -1)
amax=$(python -c "print max($a0,$a1)*1.2")

if (( $(echo "$amax < 0" | bc) == '1' )); then # assuring the -R is correct in every case
 amax=$(python -c "print $amax*-1.0")
fi

awk '{print($1,$2)}' ${PWD[0]}/*${staname}*${channel}.pf > out.grd
psxy out.grd -X0 -Y0 -R$f0/$f1/-$amax/$amax -JX10i/0.75i -B/g$amax:"Phase":Wsne:."Station ${staname}-${channel}, Hann taper ${h0%??}-${h1%??} hrs": -Wthick,black -K -O >> out.ps

for (( i=0; i<${counter}; i++ )); do
 color=$(ls -l | sed -ne "$((i+1))p" color.dat )
 awk '{print($1,$2)}' ${PWD[$i]}/*${staname}*${channel}.pfsyn | psxy -R -J -Wthick,$color,-- -K -O >> out.ps
done

#### Plotting spectrum
a0=$(awk -v f0="$f0" -v f1="$f1" '(NR>1) && ($1 > f0) && ($1 < f1)' ${PWD[0]}/*${staname}*${channel}.af | cut -f2 | sort -g | tail -1)
a1=$(awk -v f0="$f0" -v f1="$f1" '(NR>1) && ($1 > f0) && ($1 < f1)' ${PWD[0]}/*${staname}*${channel}.afsyn | cut -f2 | sort -g | tail -1)
amax=$(python -c "print max($a0,$a1)*1.2")
x0=$f0; y0=$(python -c "print -$a0*0.35")
step=$(python -c "print round(($f1-$f0)/4,2)")

awk '{print($1,$2)}' ${PWD[0]}/*${staname}*${channel}.af | psxy -X0 -Y-2.5 -R$f0/$f1/0/$amax -JX10i/2.5i -Ba$step:"f[mHz]":/g$amax:"Amplitude":WSne -Wthick,black -K -O >> out.ps

for (( i=0; i<${counter}; i++ )); do
 # Synthetic spectrum
 #a1=$(awk -v f0="$f0" -v f1="$f1" '(NR>1) && ($1 > f0) && ($1 < f1)' ${PWD[i]}/*${staname}*${channel}.afsyn | cut -f2 | sort -g | tail -1)
 #amax=$(python -c "print max($amax,$a1)*1.2")
 color=$(ls -l | sed -ne "$((i+1))p" color.dat )
 awk '{print($1,$2)}' ${PWD[$i]}/*${staname}*${channel}.afsyn | psxy -R$f0/$f1/0/$amax -J -Wthick,$color,-- -K -O >> out.ps

 # Legend
 legend1=$(echo ${PWD[$i]} | cut -d'/' -f1 | rev | cut -d'-' -f1 | rev)
 legend2=$(echo ${PWD[$i]} | cut -d'R' -f2-)
 echo -e "$x0 $y0 13 0 0 LM R=$legend2, $legend1" | pstext -R$f0/$f1/0/$a0 -N -J -K -O >> out.ps
 echo -e "$x0 $y0" | psxy -R -N -J -Ss0.07 -G$color -K -O >> out.ps
 x0=$(python -c "print $x0+($f1-$f0)/$counter")
 #echo ${PWD[i]}/*${staname}*${channel}.afsyn $legend2 $color

done

# Mode name
x0=$(python -c "print $f0*1.005")
y0=$(python -c "print $a0*0.95")
echo -e "$x0 $y0 15 0 0 LM $modename" | pstext -N -R -J -K -O >> out.ps

# plotting modes in frequency band
rm modeplot.dat
a2=$(python -c "print $amax*0.9")
awk -v f0="$f0" -v f1="$f1" '(NR>1) && ($1 > f0) && ($1 < f1)' modesin/modes.dat | awk -v a2="$a2" -v LM="LM" '{print($1" "a2" "8" "0" "0" "LM" @-"$2"@-"$3"@-"$4"@-")}' > modeplot.dat
pstext modeplot.dat -R$f0/$f1/0/$amax -J -K -O >> out.ps

echo 0 0  | psxy -R1/2/1/2 -JX1/2 -Sp -O >> out.ps
gmtset PAGE_ORIENTATION  landscape
#ps2pdf out.ps all-stations-spectrum.pdf
gv out.ps
done < stations.dat
