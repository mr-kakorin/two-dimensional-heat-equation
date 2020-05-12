reset
set term png
set output "./out.png"
plot ("./cmake-build-debug/out.txt") with points palette
