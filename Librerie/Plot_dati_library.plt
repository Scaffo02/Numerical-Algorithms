
#=============================================================================================
#
#                                          PLOT_DATI_LIBRARY.PLT
#
# Script per la generazione di grafici a partire da dati sperimentali
# contenuti in un file .dat .
#
# Per eseguire questo script, utilizzare il comando: 
#
# gnuplot Plot_dati_library.plt
#
#=============================================================================================




# Imposta terminale grafico
set terminal pngcairo size 800,600 enhanced font 'Helvetica,12'
set output 'grafico.png'

# Titolo e assi
set title "Grafico dei dati"
set xlabel "x"
set ylabel "y"

# Griglia
set grid

# Legenda in alto a destra
set key top right

# Plotta i dati dal file
plot 'trajectory.dat' using 1:2 with linespoints lw 2 lc rgb "blue" title "dati sperimentali"
