reset
set terminal pngcairo size 800,800
set output '3d-polar.png'

set lmargin at screen 0.05
set rmargin at screen 0.85
set bmargin at screen 0.1
set tmargin at screen 0.9

set pm3d map
unset key

set multiplot

# plot the heatmap
set parametric
set isosamples 500

unset border
unset xtics
unset ytics

set angles degree
R1 = 0.5
R2 = 2.9
set urange[R1:R2] # radius
set vrange[0:360] # angle
set xrange[R1:R2]
set yrange[R1:R2]
set colorbox user origin 0.9,0.1 size 0.03,0.8
splot ("./out.txt")

# now plot the polar grid only
