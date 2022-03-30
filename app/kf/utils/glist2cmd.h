*     Include file for the glist2cmd program

*  Common for edits         
      integer*4 maxedit,minobs,numu,numx
      parameter(maxedit=100)            
      character*8 ulist,xlist
      real*4 span,minlat,maxlat,minlon,maxlon      
      common /edits/minobs,span
     .            , ulist(maxedit),numu,xlist(maxedit),numx
     .            , minlat,maxlat,minlon,maxlon 


