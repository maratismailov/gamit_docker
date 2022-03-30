CTITLE G1000_BD

      block data g1000_bd

*     Block data to initialize the g1000 common block

      include 'g1000.h'

      data gvl, gvr, gvb, gvt / 0.05, 0.955, 0.05, 0.95 /
      data gwl, gwr, gwb, gwt / 0.00, 1.00, 0.00, 1.00 /

      data jcenter / 0 /
      data jor_deg / 0 /
      data jsze    / 1 /
      data jdash   / 0 /
      data jwidth  / 1 /
      data jfnt    / 2 /

      data jmode   / .false. /

*     Thats all
      end

