#! /bin/bash
rm *.ps

RUNNAME=$(basename `pwd`)
RUNDIR=`pwd` #/net/home/talavera/radial-inversion/01s00/01s00/$RUNNAME
#RUNNAME=$(basename "$RUNDIR")

if [ $# -ne 3 ]
then
        printf "USAGE: ./plot-cst.sh start_iter_number end_iter_number mode \n"
        exit
fi

its=$1; ite=$2; mode=$3

declare -a grays=(gray30 darkgray dimgray black)
declare -a color=(blue green red  orange cyan brown hotpink darkgreen purple indianred salmon firebrick)
#declare -a damping=(1.00e+04 1.00e+03 1.00e+02 1.00e+01 1.00e+00 1.00e-01 1.00e-02 1.00e-03 1.00e-04)
declare -a damping=(10000 1000 100 10 1 0.1 0.01 0.001 0.0001 0.00001)
#declare -a damping=(1.00e+03) # selected damping 

gmtset PAPER_MEDIA=a3
gmtset X_ORIGIN=1i
gmtset Y_ORIGIN=5i
gmtset HEADER_FONT_SIZE=15p
gmtset HEADER_OFFSET=0p
gmtset ANOT_FONT_SIZE=10p
gmtset LABEL_FONT_SIZE=10p

gmtset PAGE_ORIENTATION landscape
echo 0 0  | psxy -R1/2/1/2 -JX1/2 -Sp -K > out.ps

########################
### cst vs iteration ###
########################

## min & max; inversion
min=1000;max=-1000
for damp in ${damping[@]}; do # inversion
 for (( iterno=$its; iterno<=$ite; iterno++ )); do
  file=$RUNDIR/d$damp/Iter_$iterno/cst_solver/new_model.damp_$(python -c "print('{0:.2e}').format($damp)").cst
  temp_min=$(awk 'NR == 0; NR > 1 {print $1 | "sort -g"}' $file | head -1); min=$(python -c "print min($temp_min,$min)")
  temp_max=$(awk 'NR == 0; NR > 1 {print $1 | "sort -g"}' $file | tail -1); max=$(python -c "print max($temp_max,$max)")
 done
done

## min & max; previous models
cst=$(head -1 $file)
for model in HT QM1; do
  file=$RUNDIR/$mode-$model.cst
  temp_min=$(awk '{print $1 | "sort -g"}' $file | head -1); min=$(python -c "print min($temp_min,$min)")
  temp_max=$(awk '{print $1 | "sort -g"}' $file | tail -1); max=$(python -c "print max($temp_max,$max)")
 ((i++))
done

## plotting iteration=its...ite-1
min=$(python -c "print round($min*1.3,2)"); max=$(python -c "print round($max*1.25,2)")
step=$(python -c "print round(($max-$min)/4,2)"); R=R0/$(($cst+1))/$min/$max; i=0

echo 0 0 | psxy -X0 -Y0 -$R -JX15i/5i -Ba1/a$step:"cst[@%12%\155@%%Hz]"::."$mode with $RUNNAME; it=$its...$ite":WSne -Sp -K -O >> out.ps
for damp in ${damping[@]}; do
 c=(${color[$i]})
 for (( iterno=$its; iterno<=($ite-1); iterno++ )); do
  file=$RUNDIR/d$damp/Iter_$iterno/cst_solver/new_model.damp_$(python -c "print('{0:.2e}').format($damp)").cst
  awk 'NR == 0; NR > 1 {printf("%1d %s\n", NR-1, $1)}' $file | psxy -R -J -W0.1p,$c,- -K -O >> out.ps
  awk 'NR == 0; NR > 1 {printf("%1d %s\n", NR-1, $1)}' $file | psxy -R -J -Sd0.25 -G${c} -K -O >> out.ps #-Gp1/1:B${c}F-
 done
 ((i++))
done

## plotting iteration=ite as thicker line
i=0
for damp in ${damping[@]}; do
 c=(${color[$i]})
 file=$RUNDIR/d$damp/Iter_$iterno/cst_solver/new_model.damp_$(python -c "print('{0:.2e}').format($damp)").cst
 awk 'NR == 0; NR > 1 {printf("%1d %s\n", NR-1, $1)}' $file | psxy -R -J -Wthicker,$c -K -O >> out.ps
 awk 'NR == 0; NR > 1 {printf("%1d %s\n", NR-1, $1)}' $file | psxy -R -J -Sd0.5 -G$c -K -O >> out.ps
 ((i++))
done

## plotting previous models
i=0
for model in HT QM1; do
 c=(${grays[$i]})
 file=$RUNDIR/$mode-$model.cst
 awk '{printf("%1d %s\n", NR, $1)}' $file | psxy -R -J -Wthicker,$c -K -O >> out.ps
 awk '{printf("%1d %s\n", NR, $1)}' $file | psxy -R -J -Sd0.5 -G$c -K -O >> out.ps
 ((i++))
done

#############
## Legend  ##
#############

## iterations
i=0; xstep=35; ystep=6
for damp in ${damping[@]}; do
 damp=$(python -c "print 1.0/$damp")
 c=(${color[$i]})
 echo -e "0.2 1.0 12 0 0 LM $damp" | pstext -N -R0/1/0/1 -JX1/1 -Y$ystep -X$xstep -K -O >> out.ps
 echo -e "0.0 1.0 " | psxy -N -R -J -Sd0.3 -G$c -K -O >> out.ps
 ystep=-0.5; xstep=0
 ((i++))
done

## previous models
i=0
for model in HT QM1; do
 c=(${grays[$i]})
 echo -e "0.2 1.0 12 0 0 LM $model" | pstext -N -R0/1/0/1 -JX1/1 -Y$ystep -X$xstep -K -O >> out.ps
 echo -e "0.0 1.0 " | psxy -N -R -J -Sd0.3 -G$c -K -O >> out.ps
 ystep=-0.5; xstep=0
 ((i++))
done


echo 0 0  | psxy -R1/2/1/2 -JX1/2 -Sp -O >> out.ps
gmtset PAGE_ORIENTATION landscape
#ps2pdf out.ps cst.pdf
gv out.ps

