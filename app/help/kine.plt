 head 2 1
 x_f 0 1 6 "Time"
 y_f 1 7 8 "North (m)"
 poly_un day mm 1 1000
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
 axes

 key 

*** EAST
 new_win 0
 y_f 1 9 10 "East (m)"
 poly_un day mm 1 1000
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
 axes

 key

*** UP
 new_win 0
 y_f 1 11 12 "Up (m)"
 poly_un day mm 1 1000
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
 axes

 key



 
