      i2 = 0
      call getcmd(5,'aprior',wcmd,lcmd,2)
      do 170 i0 = 1,1000
         call getcmd(5,'tight_apr_sv',wcmd,lcmd,3)
         if (lcmd.le.0) goto 180
         i2 = i2+1
c        decompose command line
         ic = count_arg(wcmd)
c        not enough arguments
         if (ic.le.1) goto 170
c        pointer to
         ib = lift_arg(wcmd,code,1)
         if (ib.le.0) goto 170
         type = 1
         if (upperc(code(1:4)).eq.'ALL_') type = 2
         if (type.eq.1) snam = upperc(code(1:4))
         read (wcmd(ib+1:lcmd),*,err=500,end=500) (temp(i),i=1,norb) 
c*       read (wcmd(ib+1:lcmd),*) (temp(i),i=1,3) 
         do i=1,norb
           if( temp(i).eq.0.d0 ) then
               call report_stat('FATAL','SOLVE','get_sat_apr',' ',
     .          'SV antenna tight constraints are zero in batch file',0)
               goto 160
           endif
         enddo
 160  do 165 i=1,nsat   
c           i1 is index of first SV ant offset in parameters array 
            i1=3*(i-1)+ (indorb+norb*nsat)
            if (type.eq.1) then
               j = mchkey(rlabel(i1),snam,20,4)
               if (j.le.0) goto 165
            endif
            do j=1,3                 
c             convert units to km to agree with partials from MODEL
              sat_apr(i,norb+j) = temp(j)*.001.d0
            enddo
 165  continue
 170  continue
 180  if (i2.gt.0) then
         write( 6,185)
         write(10,185)
 185     format(/,' Apriori SV antenna offset errors (m)',/,
     .          ' Sat#    dX       dY       dZ ',/)
         do i=1,nsat
           write( 6,'(i4,3x,3f8.3)') i,(sat_apr(i,norb+j),j=1,3)
           write(10,'(i4,3x,3f8.3)') i,(sat_apr(i,norb+j),j=1,3)
         enndo
      endif

