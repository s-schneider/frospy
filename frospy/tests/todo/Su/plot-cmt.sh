#! /bin/bash
# Plotting Q vs freq
# Copy color.dat to plot
rm *.ps

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
pscoast -B60g60 -R0/360/-90/90 -JR180/25 -Dc -A10000 -W2 -Sskyblue -K > out.ps
psmeca cmt.dat -R -J -Sm0.2i -Gred -O -K >> out.ps

echo 0 0  | psxy -R1/2/1/2 -JX1/2 -Sp -O >> out.ps
gmtset PAGE_ORIENTATION landscape
#ps2pdf out.ps Q-vs-freq.pdf
gv out.ps
