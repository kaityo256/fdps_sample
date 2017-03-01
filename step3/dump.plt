unset key
set pointsize 2
set size square
do for [i=0:99] {
  set term png
  set out sprintf("plot%03d.png",i)
  plot sprintf("dump00_%03d.dat",i)\
  ,sprintf("dump01_%03d.dat",i)\
  ,sprintf("dump02_%03d.dat",i)\
  ,sprintf("dump03_%03d.dat",i)
}
