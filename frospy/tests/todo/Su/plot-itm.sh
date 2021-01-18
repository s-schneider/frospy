#! /bin/bash
rm *.ps
rm temp.dat

RUNNAME=$(basename `pwd`)
RUNDIR=`pwd` #/net/home/talavera/radial-inversion/01s00/01s00/$RUNNAME

if [ $# -ne 3 ]
then
        printf "USAGE: ./plot-cst.sh start_iter_number end_iter_number plot_title \n"
        exit
fi

its=$1; ite=$2; title=$3

declare -a grays=(gray30 darkgray dimgray black)
declare -a color=(blue green purple brown hotpink darkgreen cyan indianred orange red salmon firebrick)
declare -a damping=(0.00001 0.0001 0.001 0.01 0.1 1 10 100 1000 10000)
xmax=$(printf "%e\n" "${damping[@]}" | sort -g | tail -1); xmax=$(python -c "print $xmax*10.0")
xmin=$(printf "%e\n" "${damping[@]}" | sort -g | head -1); xmin=$(python -c "print $xmin*0.10")

gmtset PAPER_MEDIA=a4
gmtset X_ORIGIN=0.7i
gmtset Y_ORIGIN=4.5i
gmtset HEADER_FONT_SIZE=15p
gmtset HEADER_OFFSET=0p
gmtset LABEL_OFFSET=0p
gmtset ANOT_FONT_SIZE=7
gmtset LABEL_FONT_SIZE=10p
gmtset PAGE_ORIENTATION landscape

echo 0 0  | psxy -R1/2/1/2 -JX1/2 -Sp -K > out.ps

############################
## Var. Reduction vs damp ##
############################

# min & max
ymin=2;ymax=-2
for (( iterno=$its; iterno<=$ite; iterno++ )); do
 for damp in ${damping[@]}; do
  file=$RUNDIR/d$damp/Iter_$iterno/Z/synseis/Total_VR.out 
  temp=$(awk -F":" '{print $2}' $file | head -2 | tail -1)
  temp=$(python -c "print (1.-$temp/100.)")
  damp=$(python -c "print('{0:f}').format($damp)")
  echo $damp $iterno $temp >> temp.dat
  ymin=$(python -c "print min($temp, $ymin)"); ymax=$(python -c "print max($temp,$ymax)")
 done
done

ymin=$(python -c "print round($ymin*0.99,4)"); ymax=$(python -c "print round($ymax*1.01,4)")
ystep=$(python -c "print round(($ymax-$ymin)/4,4)")
R=R$xmin/$xmax/$ymin/$ymax
i=0

# plot
echo 0 0 | psxy -$R -JX3.3il/3i -Ba1:"Damping":/a$ystep:"Misfit":WS -Sp -K -O >> out.ps
for (( iterno=$its; iterno<=$ite; iterno++ )); do
 c=(${color[$i]})
 awk "/ ${iterno} / { print \$1, \$3 }" temp.dat | psxy -R -J -Sd0.2 -G$c -K -O >> out.ps 
 awk "/ ${iterno} / { print \$1, \$3 }" temp.dat | psxy -R -J -Wthick,$c -K -O >> out.ps
 ((i++))
done

echo 'Var. Reduction vs damp: Done!'

#################################
## Var. Reduction vs iteration ##
#################################

# plot
R=R$(($its-1))/$(($ite+1))/$ymin/$ymax; i=0
echo 0 0 | psxy -X0 -Y-10 -$R -JX3.3i/3i -Ba1:"Iteration":/a$ystep:"Misfit":WS -Sp -K -O >> out.ps

for damp in ${damping[@]}; do
 c=(${color[$i]})
 damp=$(python -c "print('{0:f}').format($damp)")
 awk "/${damp}/ { print \$2, \$3 }" temp.dat | psxy -R -J -Sd0.2 -G$c -K -O >> out.ps # syntax to recognize pattern
 awk "/${damp}/ { print \$2, \$3 }" temp.dat | psxy -R -J -Wthick,$c -K -O >> out.ps
 ((i++))
done

rm temp.dat
echo 'Var. Reduction vs iteration: Done!'

##########################
### Model size vs damp  ##
##########################

## min & max
ymax=0
for (( iterno=$its; iterno<=$ite; iterno++ )); do
 for damp in ${damping[@]}; do
  #file=$RUNDIR/d$damp/Iter_$iterno/cst_solver/new_model.damp_$(python -c "print('{0:.2e}').format($damp)").cst
  file=$RUNDIR/d$damp/Iter_$iterno/cst_solver/mcst.dat
  msize=$(python -c "import numpy as np; cst=np.loadtxt('$file', skiprows=1);print(np.sum(cst**2))")
  damp=$(python -c "print('{0:f}').format($damp)")
  #awk -F":" '{printf "%f %d %f \n",  "'"$damp"'", "'"$iterno"'", "'"$msize"'"}' $file | head -1 >> temp.dat
  echo $damp $iterno $msize >> temp.dat
  ymax=$(python -c "print max($msize, $ymax)")
 done
done

ymax=$(python -c "print round($ymax*1.1,2)")
ystep=$(python -c "print round(($ymax+0.1)/4,2)")
R=R$xmin/$xmax/-0.10/$ymax
i=0

# plot
echo 0 0 | psxy -X9.7 -Y10 -$R -JX3.3il/3i -Ba1:"Damping":/a$ystep:"Model size"::."$title":WS -Sp -K -O >> out.ps
for (( iterno=$its; iterno<=$ite; iterno++ )); do
 c=(${color[$i]})
 awk "/ ${iterno} / { print \$1, \$3 }" temp.dat | psxy -R -J -Sd0.2 -G$c -K -O >> out.ps 
 awk "/ ${iterno} / { print \$1, \$3 }" temp.dat | psxy -R -J -Wthick,$c -K -O >> out.ps
 ((i++))
done
echo 'Model Size vs damp: Done!'

##############################
## Model size vs iteration  ##
##############################

R=R$(($its-1))/$(($ite+1))/-0.1/$ymax; i=0

# plot
echo 0 0 | psxy -X0 -Y-10 -$R -JX3.3i/3i -Ba1:"Iteration":/a$ystep:"Model size":WS -Sp -K -O >> out.ps
for damp in ${damping[@]}; do
 c=(${color[$i]})
 damp=$(python -c "print('{0:f}').format($damp)")
 awk "/${damp}/ { print \$2, \$3 }" temp.dat | psxy -R -J -Sd0.2 -G$c -K -O >> out.ps
 awk "/${damp}/ { print \$2, \$3 }" temp.dat | psxy -R -J -Wthick,$c -K -O >> out.ps
 ((i++))
done

rm temp.dat
echo 'Model Size vs iteration: Done!'

##############################
## Eff. eigenvalue vs damp  ##
##############################

#file=$RUNDIR/d${damping[0]}/Iter_1/cst_solver/new_model.damp_$(python -c "print('{0:.2e}').format(${damping[0]})").cst
file=$RUNDIR/d${damping[0]}/Iter_1/cst_solver/mcst.dat
cst=$(head -1 $file); cst=$(python -c "print round($cst*1.1,2)")
R=R$xmin/$xmax/-0.1/$cst; i=0

# plot
echo 0 0 | psxy -X9.5 -Y10 -$R -JX3.3il/3i -Ba1:"Damping":/a1:"Effective Eigenvalue":WS -Sp -K -O >> out.ps
for (( iterno=$its; iterno<=$ite; iterno++ )); do
 for damp in ${damping[@]}; do
  file=$RUNDIR/d$damp/Iter_$iterno/cst_solver/matrix_output/R_new_model.damp_$(python -c "print('{0:.2e}').format(1./$damp)").dat
  effeig=$(python -c "import numpy as np; R=np.fromfile('$file', dtype=np.dtype('f8')); print(np.sum(R))")
  damp=$(python -c "print('{0:f}').format($damp)")
  echo $damp $effeig >> temp.dat
 done
 c=(${color[$i]})
 cat temp.dat | psxy -R -J -Sd0.2 -G$c -K -O >> out.ps
 cat temp.dat | psxy -R -J -Wthick,$c -K -O >> out.ps
 rm temp.dat
 ((i++))
done

echo 'Effective Eigenvalue vs damp: Done!'

###################################
## Eff. eigenvalue vs iteration  ##
###################################

R=R$(($its-1))/$(($ite+1))/-0.1/$cst; i=0

# plot
echo 0 0 | psxy -X0 -Y-10 -$R -JX3.3i/3i -Ba1:"Iteration":/a1:"Effective Eigenvalue":WS -Sp -K -O >> out.ps
for damp in ${damping[@]}; do
 for (( iterno=$its; iterno<=$ite; iterno++ )); do
  file=$RUNDIR/d$damp/Iter_$iterno/cst_solver/matrix_output/R_new_model.damp_$(python -c "print('{0:.2e}').format(1./$damp)").dat
  effeig=$(python -c "import numpy as np; R=np.fromfile('$file', dtype=np.dtype('f8')); print(np.sum(R))")
  echo $iterno $effeig >> temp.dat
 done
 c=(${color[$i]})
 cat temp.dat | psxy -R -J -Sd0.2 -G$c -K -O >> out.ps
 cat temp.dat | psxy -R -J -Wthick,$c -K -O >> out.ps
 rm temp.dat
 ((i++))
done

echo 'Effective Eigenvalue vs iteration: Done!'

#############
## Legend  ##
#############

# iteration
i=0; xshift=5; yshift=17; totshift=0
#i=0; xshift=5; yshift=7; totshift=0
for (( iterno=$its; iterno<=$ite; iterno++ )); do
 c=(${color[$i]})
 echo -e "0.2 1.0 9 0 0 LM it=$iterno" | pstext -N -R0/1/0/1 -JX1/1 -Y$yshift -X$xshift -K -O >> out.ps
 echo -e "0.0 1.0 " | psxy -N -R -J -Sd0.2 -G$c -K -O >> out.ps
 yshift=-0.5; xshift=0
 ((i++))
done

# damping
xshift=1; yshift=$(python -c "print abs($yshift)*($i-1)"); i=0
for damp in ${damping[@]}; do
 c=(${color[$i]}); damp=$(python -c "print 1.0/${damping[i]}")
 echo -e "0.2 1.0 9 0 0 LM damp=$damp" | pstext -N -R0/1/0/1 -JX1/1 -Y$yshift -X$xshift -K -O >> out.ps
 echo -e "0.0 1.0 " | psxy -N -R -J -Sd0.2 -G$c -K -O >> out.ps
 yshift=-0.5; xshift=0
 ((i++))
done


echo 0 0  | psxy -R1/2/1/2 -JX1/2 -Sp -O >> out.ps
#ps2pdf out.ps Q-vs-freq.pdf
rm temp.dat
cp out.ps damping_$title.ps
gv out.ps

