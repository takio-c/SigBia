#plot "data.txt" using 1 title "sensor" with line
plot "data.txt" using 2 title "signal" with line
replot "data.txt" using 3 title "bias" with line
replot "data.txt" using 4 title "status signal" with line
replot "data.txt" using 5 title "P for signal" with line
replot "data.txt" using 6 title "P for bias" with line

