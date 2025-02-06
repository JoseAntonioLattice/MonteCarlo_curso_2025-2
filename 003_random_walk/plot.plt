reset

#f(x) = a+b*x  
#f(x) = 1 - b*sqrt(x)
f(x) = a - b*x**c
#f(x) = a - b*log(c*x+d)


d = 1

set terminal epslatex color colortext standalone
set output sprintf('%dd_random_walk.tex',d)


filename = sprintf('%dd_probabilities.dat',d)
fit [0:0.05]f(x) filename u (1/$1):2:3 via a,b,c

set title sprintf('%d$d$ random walk',d) font ',20'
set xlabel '$1/N$' font ',16'
set ylabel 'probability of return' font ',16'
set key font ',12' at graph 1,0.9

p filename u (1/$1):2:3 w e notitle lw 2 pt 6 lc 7, f(x) lw 2 lc 8 t sprintf("($ %1.4f \\pm %1.4f) + (-%1.4f \\pm %1.4f) N^{-%1.4f\\pm%1.4f}$\n\n$\\chi^2/d.o.f = %1.4f$",a,a_err,b,b_err,c,c_err,FIT_STDFIT**2)
unset output
