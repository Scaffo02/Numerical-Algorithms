# ----------------------------------------------------------------------
# Script Gnuplot per visualizzare la gittata in funzione dell'angolo
# ----------------------------------------------------------------------

# Impostazioni dell'output (scegli il formato che preferisci, ad esempio png)
set terminal pngcairo enhanced font "Arial,12"
set output 'gittata_vs_angolo.png'

# Titolo e etichette
set title "Gittata in funzione dell'Angolo di Tiro (v0=200 m/s, B=0.0005)"
set xlabel "Angolo di Tiro [Gradi]"
set ylabel "Gittata [metri]"

# Griglia
set grid

# Imposta l'asse X da 0 a 90 gradi
set xrange [0:90]

# Disegna i dati
plot 'gittata_vs_angolo.dat' using 1:2 with lines linewidth 2 title "Gittata Calcolata", \
     (1000.0) with lines dashtype 2 linewidth 1 title "Distanza Bersaglio (1000 m)"