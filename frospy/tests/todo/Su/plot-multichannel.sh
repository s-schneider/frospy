#/bin/bash
# Plotting the spectrum of one station channels
rm *.grd
rm *.ps

modename=05s03 # 01s05-02s04 01s06-02s05 02s10-04s05
#PWD1=$modename-self/060994A/R-0.25
#PWD2=$modename-self/060994A/R0
#PWD3=$modename-self/060994A/R0.25
PWD4=$modename-cross/060994A/R-0.25
PWD5=$modename-cross/060994A/R0
PWD6=$modename-cross/060994A/R0.25
#PWD7=$modename-msrd/060994A/R-0.25
#PWD8=$modename-msrd/060994A/R0
#PWD9=$modename-msrd/060994A/R0.25

#modename=$(echo $PWD1 | cut -d'-' -f-2)

declare -a PWD=($PWD4 $PWD5 $PWD6) #$PWD1 $PWD2 $PWD3 $PWD4 $PWD5 $PWD6 $PWD7 $PWD8 $PWD9)
counter=${#PWD[@]}

ls ${PWD[0]}/*.pf | rev |cut -d'.' -f2 | rev | cut -d'_' -f1 > stations.dat

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

f0=$(sed -n '3p' ${PWD[0]}/*${staname}*${channel}.window | cut -f1)
f1=$(sed -n '3p' ${PWD[0]}/*${staname}*${channel}.window | cut -f2)

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
psxy out.grd -X0 -Y0 -R$f0/$f1/-$amax/$amax -JX10i/0.75i -B/g$amax:"Phase":Wsne:."Station ${staname}-${channel}, Hann taper ${h0%??}-${h1%??} hrs": -Wthicker,black -K -O >> out.ps

for (( i=0; i<${counter}; i++ )); do
 color=$(ls -l | sed -ne "$((i+1))p" color.dat )
 awk '{print($1,$2)}' ${PWD[$i]}/*${staname}*${channel}.pfsyn | psxy -R -J -Wthicker,$color,-- -K -O >> out.ps
done

#### Plotting spectrum
a0=$(cut -f2 ${PWD[0]}/*${staname}*${channel}.af | sort -g | tail -1)
a1=$(cut -f2 ${PWD[0]}/*${staname}*${channel}.afsyn | sort -g | tail -1)
amax=$(python -c "print max($a0,$a1)*1.1")
x0=$f0
y0=$(python -c "print -$a0*0.7")

awk '{print($1,$2)}' ${PWD[0]}/*${staname}*${channel}.af | psxy -X0 -Y-2.5 -R$f0/$f1/0/$amax -JX10i/2.5i -B/g$amax:"Amplitude":Wsne -Wthicker,black -K -O >> out.ps

for (( i=0; i<${counter}; i++ )); do
 # Synthetic spectrum
 #a1=$(cut -f2 ${PWD[i]}/*${staname}*${channel}.afsyn | sort -g | tail -1)
 #amax=$(python -c "print max($amax,$a1)*1.1")
 color=$(ls -l | sed -ne "$((i+1))p" color.dat )
 awk '{print($1,$2)}' ${PWD[$i]}/*${staname}*${channel}.afsyn | psxy -R$f0/$f1/0/$amax -J -Wthicker,$color,-- -K -O >> out.ps

 # Legend
 legend1=$(echo ${PWD[$i]} | cut -d'-' -f3 | cut -d'/' -f1)
 legend2=$(echo ${PWD[$i]} | cut -d'R' -f2-)
 echo -e "$x0 $y0 10 0 0 LM R=$legend2, $legend1" | pstext -R$f0/$f1/0/$a0 -N -J -K -O >> out.ps
 echo -e "$x0 $y0" | psxy -R -N -J -Ss0.07 -G$color -K -O >> out.ps
 x0=$(python -c "print $x0+($f1-$f0)/$counter")
done

# Mode name
x0=$(python -c "print $f0*1.0005")
y0=$(python -c "print $a0*0.85")
echo -e "$x0 $y0 15 0 0 LM $modename  " | pstext -N -R -J -K -O >> out.ps

#### Plotting Singlets
amax=$(tr -s ' ' < ${PWD[0]}/*${staname}*.singlets | cut -d' ' -f3 | sort -g |tail -1)
amax=$(python -c "print $amax*1.1")
step=$(python -c "print round(($f1-$f0)/4,2)")

awk '{print($1,$2)}' ${PWD[0]}/*${staname}*.singlets > out.grd
awk '{print($1,0)}' ${PWD[0]}/*${staname}*.singlets >> out.grd # Adding zero amplitude data to plot columms
sort out.grd -o out.grd
psxy out.grd -X0 -Y-0.75 -R$f0/$f1/0/$amax -JX10i/0.75i -Ba$step:"f[mHz]":/g$amax:."":WSne -Sc0.001 -K -O >> out.ps

singlets=$(wc -l < out.grd)
for (( i=2; i<=$singlets; i=i+2 )); do
 ls -l | sed -ne "$((i-1)), ${i}p" out.grd | psxy -R -J -Wthicker,black -K -O >> out.ps
done

echo 0 0  | psxy -R1/2/1/2 -JX1/2 -Sp -O >> out.ps
gmtset PAGE_ORIENTATION  landscape
#ps2pdf out.ps all-stations-spectrum.pdf
gv out.ps
done < stations.dat
