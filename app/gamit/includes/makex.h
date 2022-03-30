c     makex.h
c       Declarations for modules that use dimensions for epochs
c       or observation types:  /lib, /clean, /makex, /clean, /ctox, /fica, /utils  
c       Last changed by M. Floyd 4 April 2019

c     file names
      character*120     fsited,fcoord,fclock,fsp3,frinex,fbatch,
     &                 fsceno,fsvclk,fxfile,finfor,fficaf,fnav,
     &                 fextr3,fanthi 

c     logical units for each file
      integer*4  usited,ucoord,uclock,usp3,urinex,ubatch,
     &           usceno,usvclk,uxfile,uinfor,uficaf,unav,
     &           uextr3,uscren,uanthi
                                
c     flags as to whether or not to use these files .true. if open
      logical    qsited,qcoord,qclock,qsp3,qrinex,
     &           qsceno,qsvclk,qxfile,qinfor,qficaf,qnav,
     &           qextr3,qanthi
                                        

C     max number of channels in receiver
      integer*4   maxchn                        
C MND TAH 122720: Increased from 72 to 80.
* MOD TAH 201215: Increased from 80 to 99                      
      PARAMETER   (MAXCHN=99)                
C     max number of epochs in X or C file
      integer*4   maxepc
      PARAMETER   (MAXEPC=2880)             
C     max elements in fica arrays
      integer*4   maxfic 
      PARAMETER   (MAXFIC=500)              
C     max RINEX observation types
      integer*4   maxobt 
C MOD TAH 120201: Increased to 27, max allowed by current code              
C MOD MAF 190404: Increased to 63, max allowed by 1024-character line length in gamit/lib/rrinex.f
      PARAMETER   (MAXOBT=63)                
C     max lines in X-file header
      integer*4   maxlin
      PARAMETER   (MAXLIN=40)               
c     max RINEX or FICA observation files
      integer*4   maxfil
      PARAMETER   (MAXFIL=100)  
c     max segments in an SV clock fit 
      integer*4   maxseg
      PARAMETER   (MAXSEG=10)

      common /fileinfo/fsited,fcoord,fclock,fsp3,frinex,fbatch,
     &                 fsceno,fsvclk,fxfile,finfor,fficaf,fnav,
     &                 fextr3,fanthi,
     &                 usited,ucoord,uclock,usp3,urinex,ubatch,
     &                 usceno,usvclk,uxfile,uinfor,uficaf,unav,
     &                 uextr3,uscren,uanthi,
     &                 qsited,qcoord,qclock,qsp3,qrinex,
     &                 qsceno,qsvclk,qxfile,qinfor,qficaf,qnav,
     &                 qextr3,qanthi 
