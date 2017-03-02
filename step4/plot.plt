set term png
set out "trajectory.png"

set style data line
set xla "time"
set yla "x"
p [][0:10] "test.dat" u 1:2 lw 3 t "Atom-1","test.dat" u 1:3 lw 3 t "Atom-2"

set out "energy.png"
unset key
set yla "Total Energy"
p [][0:] "test.dat" u 1:4 
