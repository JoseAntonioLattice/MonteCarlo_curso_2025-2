reset

f(x) = 1 - a*x**0.5

fit f(x) 'probabilities.dat' u (1/$1):2:3 via a

set title '1d random walk' font ',20'
set xlabel '1/N' font ',16'
set ylabel 'probability of return' font ',16'
set key font ',14'
p 'probabilities.dat' u (1/$1):2:3 w e notitle lw 2 pt 6 lc 7, f(x) lw 2 lc 8 t sprintf("1 - %1.4f x^{1/2}",a)
