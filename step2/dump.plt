set term png
set out "dump.png"

unset key
set pointsize 2
p "dump0.dat"\
, "dump1.dat"\
, "dump2.dat"\
, "dump3.dat"
