  view 0.1  0.95   0.3  0.95

  font 1
  x_field 1 1 0
  y_field 1 6 0
  read
  line 1 
  point 0
  draw
  xmx -1 2 "Time (sec)"
  xmn -1 0
  ymx -1 0 
  ymn -1 1 "RESIDUAL (microsec)  Cubic fit"
**** the following lines will  restart scaling !!!! 
  reset_SC on
  read
   y_scale  0 .5
  view 0.1 0.95 0. 0.6
  label 0 0.02  0 0 \h10
  label 0 0.04  0 0 \h9
  label 0 0.06  0 0 \h8
  label 0 0.08  0 0 \h7
  label 0 0.10  0 0 \h6
  label 0 0.12  0 0 \h5
  label 0 0.14  0 0 \h4
  label 0 0.16  0 0 \h3
  label 0 0.18  0 0 \h2
  label 0 0.20  0 0 \h1
