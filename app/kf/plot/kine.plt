 head 2 1
 x_f 0 1 6 "Time"
 y_f 1 7 8 "North (m)"
 poly_un day mm 1 1000
 new 20 20 600 900
 view 0.15 0.85 0.1 0.9
 read
 
 pen 10
 poi 1
 draw
 fit 0 0 A
 pen 11
 fit 0
 pdr

 pen 17
 label %5 97 1 0 :F
 label %5 95 1 0 :P2
 pen 1
 xmn -1 1 "Time"
 ymn -1 1
 xmx -1 1 "Time"
 ymx -1 1 "dNorth (mm)"

 key 
