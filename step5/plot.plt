set term png
set out "energy.png"

set xla "time"
set yla "Total Energy"
p "test.dat" pt 6 ,"../step4/test.dat" u 1:4

set out "trajectory.png"

p "trac00.dat" t "Process 0", "trac01.dat" t "Process 1"\
,"../step4/test.dat" u 1:2 w l lt 1 lc rgb "black" t "Serial"\
,"../step4/test.dat" u 1:3 w l lt 1 lc rgb "black" t ""\


