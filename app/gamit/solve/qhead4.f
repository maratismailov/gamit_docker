      Subroutine qhead4(jb,jp,usig,ssig,symbl,lwl,lw,mode)
c
c     Called by GETWL to print the resolution of wide-lane ambiguities
          
c       jb : bias number (index )
c       jp: parameter number (index)

c       usig : unscaled sigma
c       ssig : scaled sigma 

c       symbl : asterisk or blank indicating estimated or fixed parametes

c       lwl :  total number of WL biases
c       lw  :  WL bias number (index)

c       mode = 1: Unestimated biases, line by line (called before and after fixing) 
c       mode = 2: Estimated biases, line by line (called before and after fixing)
c       mode = 3: Bias list from AUTCLN N-file

      implicit none
    
c  Nearly global common blocks

      include '../includes/dimpar.h'
      include 'solve.h'
      include 'parameters.h'
                                       
c   Local 

      integer jb,jp,lwl,lw,mode 

      real*8 usig,ssig 

      character*1 symbl


c     Print the WL estimates, line-by-line

      if (mode.eq.1) then
c       unestimated biases (dependent or no data)
         IF (jb.GT.lwl) THEN
            if( logprt ) WRITE (6,10)
     .         jp,SYMBL,RLABEL(jp),ADJUST(jp),DDWL(LW),DDWV(LW)
            WRITE (10,10)
     .         jp,SYMBL,RLABEL(jp),ADJUST(jp),DDWL(LW),DDWV(LW)
         else
            if( logprt ) WRITE (6,10)
     .         jp,SYMBL,RLABEL(jp),ADJUST(jp)
            WRITE (10,10)
     .         jp,SYMBL,RLABEL(jp),ADJUST(jp)
         endif
 10   format (1x,i4,a1,a20,2x,f11.1,22x,2f10.3)

      elseif (mode.eq.2) then
c       estimated biases
         if (jb.gt.lwl) then
            if( logprt ) write (6,20)
     .        jp,symbl,rlabel(jp),adjust(jp),usig,ssig,ddwl(lw),ddwv(lw)
            write (10,20)
     .        jp,symbl,rlabel(jp),adjust(jp),usig,ssig,ddwl(lw),ddwv(lw)
         else
            if( logprt ) write (6,20)
     .         jp,symbl,rlabel(jp),adjust(jp),usig,ssig
            write (10,20)
     .         jp,symbl,rlabel(jp),adjust(jp),usig,ssig
         endif
 20   format (1x,i4,a1,a20,2x,f11.3,2f10.3,2x,2f10.3)

          
      elseif (mode.eq.3) then 
c       biases passed from AUTCLN   
        if( logprt ) write(6,30)  jp,rlabel(jp),wlval(jb),symbl
        write(10,30) jp,rlabel(jp),wlval(jb),symbl   
 30     format (1x,i4,1x,a20,6x,f11.3,1x,a1)

      else    
        call report_stat('FATAL','SOLVE','qhead4',' '
     .          , 'Invalid mode in bias printout',0) 

      endif 

      return
      end


