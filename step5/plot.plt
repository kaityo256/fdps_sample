set term png
set out "energy.png"

set xla "time"
set yla "Total Energy"
p "serial.dat" t "Serial","test.dat" pt 6 t "2-MPI"

set out "trajectory.png"

p "trac00.dat" t "Process 0", "trac01.dat" t "Process 1"\
,"serial_trac.dat" u 1:2 w l lt 1 lc rgb "black" t "Serial"\
,"serial_trac.dat" u 1:3 w l lt 1 lc rgb "black" t ""\


