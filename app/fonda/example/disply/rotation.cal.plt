* make rotation map in California network by kplot

 head 3 1 
 file /data10/tah/jplbrad/vg5_res.str 
.
 mproj me 35 240 0
 mset  co 32  238 38 245  
 line 21
 mdraw  ou us sa 3 la n gr  2
 char 2.0

*prepare to read strains and omega
 sfield 2 1 3 4 5 6 7 8 9 10
 x_field 1 1 0
 sread
 line 10
 sdraw 2 5

 char 4
 line 21
 label 239.5 37.8  1  0 "Residual rotation Field"
 
*put a scaled rotation fan
 char 2.0
 file /home/chandler/dong/fonda/velo_model/scale.rot
 head 4
 sfield 2 1 3 4 5 6 7 8 9 10
 x_field 1 1 0
 sread
 sdraw 2  5
 label  238.3 32.8 1 0 "1.e-7rad/yr"

* draw network
 head 4
 file /data10/tah/jplbrad/vg5site.map
 line 21
 v_field 1 2 3 4 5 6 7 0
 x_field 1 1 0
 vread
 line 1
 network

* Now put some Faults on the plot
 head 3 1
 file /home/chandler/tah/kf/plot/faults.line
 x_f 1 2 0
 y_f 1 1 0
 p_f   3 0
 read
 line 3  1  
 poi 0
 draw


* Label commands for site names
 char 2.0
 line 1
 font 1
 label  243.112  35.332 1 0 "MOJM"
 label  242.270  33.400 1 0 "NIGU"
*label  241.829  34.205 1 0 "JPL1"
 label  241.706  37.233 1 0 "OVRO"
*label  241.596  33.744 1 0 "PVER"
*label  241.595  33.3   1 0 "BRSH"
 label  241.381  32.800 1 0 "BLUF"
*label  241.399  34.330 1 0 "SAFE"
*label  241.331  34.496 1 0 "LOVE"
*label  241.214  34.086 1 0 "CATO"
*label  241.150  34.358 1 0 "HAPY"
*label  241.134  34.478 1 0 "HOPP"
*label  241.001  34.388 1 0 "SNPA"
*label  240.990  34.440 1 0 "SNP2"
*label  240.961  34.326 1 0 "SCLA"
*label  240.846  34.120 1 0 "COTR"
*label  240.699  34.636 1 0 "MUNS"
*label  240.657  34.298 1 0 "SOLI"
 label  240.606  35.420 1 0 "FIBR"
 label  240.400  33.132 1 0 "TWIN"
*label  240.286  34.494 1 0 "LACU"
*label  240.247  33.850 1 0 "CENT"
*label  239.933  35.076 1 0 "MADC"
*label  239.801  34.502 1 0 "GAVI"
*label  239.701  35.400 1 0 "POZO"
*label  239.586  34.731 1 0 "GRAS"
*label  239.394  34.894 1 0 "LOSP"
 label  239.050  34.556 1 0 "VNDN"
*label  239.168  35.259 1 0 "BLHL"
 label  238.400  35.565 1 0 "BLAN"
 label  237.900  36.570 1 0 "FTOR"


 end
