set term pngcairo size 1800,600
set output "solutions_all.png"

# IMPOSTAZIONE COLORI ARCOBALENO
# 33,13,10 è la formula standard per l'arcobaleno (Blue-Green-Yellow-Red)
set palette rgbformulae 33,13,10

set pm3d map
set view map
set size ratio -1

# Ignora righe vuote
set datafile commentschar "#"

set multiplot layout 1,3 title "Soluzioni Poisson 2D (Jacobi, Gauss-Seidel, SOR)"

#######################
#      JACOBI
#######################
set title "Jacobi"
set xlabel "x"
set ylabel "y"
splot "solution.dat" using 1:2:3 with pm3d notitle

#######################
#      GAUSS-SEIDEL
#######################
set title "Gauss-Seidel"
set xlabel "x"
set ylabel "y"
splot "solution.dat" using 1:2:4 with pm3d notitle

#######################
#      SOR
#######################
set title "SOR (omega)"
set xlabel "x"
set ylabel "y"
splot "solution.dat" using 1:2:5 with pm3d notitle

unset multiplot