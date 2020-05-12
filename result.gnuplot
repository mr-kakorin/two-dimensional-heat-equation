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
r = 0.05 # CHANGED THIS
set urange[0.5:2.9] # radius
set vrange[0:360] # angle
set xrange[0.5:2.9]
set yrange[0:6.28]
set colorbox user origin 0.9,0.1 size 0.03,0.8
set dgrid3d         # ADDED THIS
splot "./cmake-build-debug/out.txt" using 1:2:3

# now plot the polar grid only
set style line 11 lc rgb 'white' lw 2
set grid polar ls 11
set polar
set rrange[0:r]
unset raxis
set rtics format '' scale 0
unset parametric
set for [i=0:330:30] label at first (r+0.35)*cos(i), first (r+0.35)*sin(i)\
center sprintf('%d', i)
plot NaN w l
unset multiplot
