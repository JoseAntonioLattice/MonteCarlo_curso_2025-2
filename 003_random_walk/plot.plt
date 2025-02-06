reset

f(x) = 1 - b*sqrt(x)
#f(x) = a - b*x**c
#f(x) = a - b*log(c*x+d)
d = 2
fit [0:0.01]f(x) 'probabilities.dat' u (1/$1):2:3 via b

set title sprintf('{%d}d random walk',d) font ',20'
set xlabel '1/N' font ',16'
set ylabel 'probability of return' font ',16'
set key font ',10' at 0.1,0.7
#set label at 0.033,0.32 font ',14' sprintf('{/Symbol c}^2/d.o.f = %1.4f',FIT_STDFIT**2)
p 'probabilities.dat' u (1/$1):2:3 w e notitle lw 2 pt 6 lc 7, f(x) lw 2 lc 8 t sprintf("(%1.4f+/-%1.4f) - (%1.4f+/-%1.4f)N^{-%1.4f+/-%1.4f}\n{/Symbol c}^2/d.o.f = %1.4f",b,b_err,FIT_STDFIT**2)
