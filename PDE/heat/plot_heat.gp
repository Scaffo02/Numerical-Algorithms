# Gnuplot script per salvare le soluzioni con il massimo ingrandimento possibile
# RAPPORTO ASSI X:Y = 2:1 con palette arcobaleno

# Impostazioni generali
set title font ",14"
set xlabel "x"
set ylabel "y"
set view 0, 0
set style data pm3d
set style function pm3d
set pm3d map

# Margini ridotti per massimizzare l'area di plottaggio
set lmargin 5  
set rmargin 1  
set bmargin 4  
set tmargin 1  

# --- IMPOSTA LA PALETTE ARCOBALENO PURO ---
set palette defined ( 0 "blue", 1 "cyan", 2 "yellow", 3 "red" ) 

# Impostazioni range e etichette
set cbrange [0:2.0]
set cblabel "{\/Symbol y}" 

# L'area del plot occuperà quasi tutto il terminale, che è 600x300 (ratio 2:1)
set size 1.0, 1.0
set origin 0.0, 0.0


# --- 1. Grafico Jacobi ---
set terminal pngcairo size 600, 300 
set output 'soluzione_jacobi_RAINBOW.png'
set title "Metodo Jacobi (\psi)"
splot "solution.dat" using 1:2:3 notitle


# --- 2. Grafico Gauss-Seidel ---
set terminal pngcairo size 600, 300
set output 'soluzione_gauss_seidel_RAINBOW.png'
set title "Metodo Gauss-Seidel (\psi)"
splot "solution.dat" using 1:2:4 notitle


# --- 3. Grafico SOR ---
set terminal pngcairo size 600, 300
set output 'soluzione_sor_RAINBOW.png'
set title "Metodo SOR (\psi)"
splot "solution.dat" using 1:2:5 notitle