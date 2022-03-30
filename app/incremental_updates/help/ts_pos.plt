*
* CPLOTX command file to plot .pos files.
*
* Run with:
* % cplotx ~/gg/help/ts_pos.plt '' <tssum.pos file> 37 1 
*
* VERSION: Plot of single series. Use ts_pos2.plt to plot overlay.
*
 new_wind 20 20 600 800
 ps_size 200 275
 head 37 1
 file #1
 x_f 3 3 5   "Date"
 y_f 1 16 19  "dNorth (m)"
 read 
*
* Plot detrended data (to fit 0 1 A to remove only mean)
*
 fit 1 1 A 
 view 0.125 0.90 0.65 0.95
 axes -x

 pen 10
 poi 3
 char 3
 draw
 fit 0
 pen 1
 pdraw
 label %5 100.5 1 0 :F
 label %5  96   1 0 :H4
 label %5   5   1 0 :P2
*
*----------------------------------------------
* Read east component
 
 y_f 1 17 20 "dEast (m)"
 file #1
 read 
* Plot detrended data (to fit 0 1 A to remove only mean)
 fit 1 1 A 
 view 0.125 0.90 0.35 0.65 
 pen 1
 axes -x 

 pen 10
 poi 3
 char 3
 draw
 stat
 fit 0
 pen 1
 pdraw
 label %5   5   1 0 :P2

*----------------------------------------------
* Read height component
 y_f 1 18 21 "dHeight (m)"
 file #1
 read 
* Plot detrended data (to fit 0 1 A to remove only mean)
 fit 1 1 A 
 view 0.125 0.90 0.05 0.35
 pen 1
 axes

 pen 10
 poi 3
 char 3
 draw
 fit 0
 pen 1
 pdraw
 label %5   5   1 0 :P2
*
* use END or RET commands to finish. Other
* cplotx commands can not be used on the last 
* frame plotted.
 key 
