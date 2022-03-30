Copyright (c) Massachusetts Institute of Technology,1986. All rights reserved.
      Subroutine start_int ( s,l2 )
c
c           A.Rasinski, Dec.1964 (rev. 1024-68) Subroutine NINT in PEP
c           revised July 1965.  branch for FN added oct 67, nov 69.
c           revised Sept 1969   variable tabular interval (itype=0)
c           revised Dec  1979   step-size control   RDR/KCL
c           Rick Abbot Oct 1984  modification to do GPS satellites only 
c           renamed start_int Oct 96  R. King

      implicit none

      include '../includes/dimpar.h'   
      include '../includes/arc.h'
        
      integer*4 jj,idim,l2,i,j

      real*8 hr,tdif,sbfn,s,dir
            
      logical debug/.false./

      save dir

c        l1
c        -1  first step of continuing procedure  (rwk 1101014: this seems never to happen)
c         0  starting procedure
c         1  continuing procedure
c
c        k2      count of steps -- used to find end of
c                starting procedure at k2=19 in eval
c                and used to control accuracy tests
c        l4      if .lt. 0 step size not to be changed
c                used to prevent frequent step size changes
c        l5      used in tabular output mode
c                -1 do not change step size
c                set to -1 in tabular output mode except
c                for first step after each tabular point.
c                this prevents the need for a fractional
c                step to get to the next tabular point.
c        n1, n2  range of equation number for integration.
c                first (1 to meq), second(meq+1 to neq).
c        s       integration time variable, seconds from start
c                    
      if(debug) print *,'START INIT enter l1 l2',l1,l2
      if (l1)  20,10,21
c
c        starting procedure
   10 dir = tstop - tint0
      hc = dsign(hc,dir)
      l4 =-30
      k2 = 0
      n1 = 1
      n2 = neq
      l2 = 0
      s  = tint0
      if(debug) print *,
     .   'START_INT l1 l2 l4 tstop hc s ',l1,l2,l4,tstop,hc,s
      do i=1,neq
        anord(i) = 0.0d0
        bnord(i) = 0.0d0
        cnord(i) = 0.0d0
        dnord(i) = 0.0d0
        enord(i)=0.d0
        y(i,3) = satics(i)
        y(i,4) = satics(i) 
      enddo
      idim = 3
      do i=1,neq
        fnord(i,3) = sbfn(i,idim,s)
        if(debug) print *,'fnord3 i ',fnord(i,3) 
      enddo              

c        integration started with 20 steps, reversing sense after
c        every 5 steps
      do 18 j=1,4
       do 15 jj=1,5
        if(debug)  print *,'START_INT calling EVAL l2 j jj s ',l2,j,jj,s
        call eval (s,l2)
        if(debug) print *,'after EVAL s l2 ',s,l2 
        if (l2) 10,15,999
   15  continue  
       hc = - hc 
       if(debug) print *,'START_INT l2=0 s hc tint0 ',s,hc,tint0
       if ((s - tint0).eq.0.0d0) then
         do i=1,neq
           y(i,3) = satics(i)
           y(i,4) = satics(i)
         enddo
       endif
   18 continue
c
c        first step of continuing procedure
   20 tab=tint0+hmx     
      if(debug) print *, 'START_INT 1st step nsign dir hc hr hmx '
     .                  , nsign,dir,hc,hr,hmx  
      l1 =1
      l5 =0
      dir=nsign
      do  i=1,neq
        y(i,3)=satics(i)
        y(i,4)=satics(i)
      enddo
      idim = 3
      do i=1,neq
        fnord(i,3)=sbfn(i,idim,s)
        dydt(i,1)=fnord(i,3)
        if(debug) print *,'dydt1 i ',dydt(i,1) 
      enddo      
      
c        continuing procedure
  21  hr=0.0d0      
      if(debug) print *,'START_INT continue step nsign dir hc hr hmx '
     .                  , nsign,dir,hc,hr,hmx 
c
c        from here to call to eval is logic to keep
c        from passing a tab point
c
   24 n1 = 1
      n2 = meq
      l2 = 0
      tdif = tab - (s+hc)         
cd      print *,'At 24 s hc tab tdif ',s,hc,tab,tdif
c
c     ignore the following conditional if dir<=0 and tdif<0, or
c     dir>0 and tdif>0.  (note that the conditional is entered
c     whenever tdif=0.)
c
      if ((dir.gt.0.0d0 .and. tdif.le. 1.0d-16) .or.
     1    (dir.le.0.0d0 .and. tdif.ge.-1.0d-16)) then
c
c      a step of hc will land on or pass the tab point
         if (l4.ge.0) l4 = -1                     
         if(debug) print *,'step will land or pass '
         if (dabs(tdif).ge.1.0d-16) then    
c
c      a step of hc will pass the tab point
           hr = hc
           hc = tab - s
           if(debug) print *,'START_INT step will pass s hc tab '
     .           ,s,hc,tab
         endif
      endif
       
      if(debug) print *,'START_INT continuing EVAL s l2 s l2 hc hr '
     .                  ,s,l2,hc,hr
      call eval (s,l2)
                   
      if(debug) print *,'START_INT tab check l2 s tab ',l2,s,tab 
      if (l2.lt.0) goto 24
      if (l2.eq.0) then
         if (tab.ne.s) goto 24
c
c        have just hit a tab point.  undo the special
c        conditions
c
         tab = tab + hmx
         if (hr.ne.0.0d0) hc = hr
         if (l5.eq.0) l5 = -1
         if (l5.lt.0) l5 = +1 
         if(debug) print *,'hit a tab l5 s hmx hc hr ',l5,s,hmx,hc,hr 
      endif
c
  999 return
c
      end
