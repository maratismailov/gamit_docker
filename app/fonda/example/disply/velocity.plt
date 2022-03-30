*plot velocity field in Southern California
 mproj me 33.5 243.3 0
 mset co 32.4 241.6 34.6 245
* mset co 32.5 242 34.5 245
 mdraw ou us la n
* mreset
 
 head 4 
 file /home/chandler/dong/fonda/agnewd/testb.map
* view .1 .95 .1 .95
 x_f 1 1 0 "longitude"
* x_scale 242.0 245.0
 y_f 1 2 0 "latitude"
* y_scale 32.5 34.5
 font 3
 charsz 2.0 3.0
 v_f 1 2 3 4 5 6 7 8
 reset off
 vread
 vdraw 1.00 0.393
 line 1
 

*xmn -1 1 "longitude"
*ymn -1 1 "latitude"

 label 243.8 34.4 1 0 \h1
 label 243.8 34.3 1 0 \h2
 label 243.8 34.2 1 0 \h3
x label 1.8 3.1 1 0 \h4
 label 243.8 34.1 1 0 "  1 sigma ellipse"
* axes
 label 0 -3 1 0 \t

 head 2 
 file /home/chandler/dong/fonda/agnewd/scale.vel
 x_f 1 1 0
 v_f 1 2 3 4 0 0 0 8
 reset off
 line 1
 vread
 vdraw 1.00 0.393

 head 2 
 file /home/chandler/dong/fonda/agnewd/vlbi.vel
 x_f 1 1 0
 v_f 1 2 3 4 5 6 7 8
 reset off
 line 2
 vread
 vdraw 1.00 0.393

.

