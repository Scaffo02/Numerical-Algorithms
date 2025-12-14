# ----------------------------------------------------------------------
# Script Gnuplot per visualizzare le traiettorie di tiro basso e alto
# Questo script salva il risultato in un file PNG.
# ----------------------------------------------------------------------

# Impostazioni di output per salvare l'immagine in PNG
set terminal pngcairo enhanced font "Arial,12"
set output 'traiettorie_tiro.png'

# Rimuovi/commenta l'opzione per la visualizzazione a schermo:
# set terminal wxt size 800, 600

# Titolo e etichette
# Attenzione: V0 nel codice precedente era 200 m/s, lo script è stato aggiornato per coerenza.
set title "Confronto Traiettorie con Attrito dell'Aria (V0=200m/s, Target=1000m)"
set xlabel "Distanza Orizzontale X (m)"
set ylabel "Altezza Y (m)"

# Griglia
set grid

# Impostazioni del range degli assi
set xrange [0:1100]  # Da 0 metri a leggermente oltre il bersaglio
set yrange [0:*]     # Da 0 metri in su (y positiva)

# Plot dei dati
plot \
    'traj_bassa.dat' using 1:2 with lines linewidth 2 lc rgb "red" title "Tiro Basso", \
    'traj_alta.dat' using 1:2 with lines linewidth 2 lc rgb "blue" title "Tiro Alto", \
    'target.dat' using 1:2 with points pointtype 3 pointsize 2 lc rgb "black" title "Bersaglio (1000m)"

# Rimuovi/commenta la pausa, non necessaria quando si salva su file:
# pause -1 "Premi invio per uscire."