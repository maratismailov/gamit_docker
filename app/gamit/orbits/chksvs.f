Copyright (c) Massachusetts Institute of Technology,1987. All rights reserved.
      subroutine chksvs( nesat,iesat,nsat,ixsat,iexpnt )
C
C Written by R. King May 1989
C
C  Input:
C        nesat              number of satellites in the E-file
C        iesat              PRN numbers for E-file satellites
C        nsat               number of satellites in the X-file
C        ixsat(nsat)        PRNs numbers for X-file satellites
C
C  Output:
C        iexpnt(nsat         pointers from X-file satellites to positions in
C                             the E-file arrays:  ixsat(i)= iesat( iexpnt(i) )
C
C
      implicit none

      include '../includes/dimpar.h'

      integer*4 nesat,iesat,nsat,ixsat,iexpnt,i,j
      character*256 message

      dimension iesat(MAXSAT),ixsat(MAXSAT),iexpnt(MAXSAT)
C
C  Determine the X-file pointers to the E-file satellites

      do i=1,maxsat
         iexpnt(i) = 0
      enddo
      do i=1,nsat
         do j=1,nesat
           if( iesat(j).eq.ixsat(i) ) iexpnt(i) = j
         enddo
      enddo


C  Locate and issue warning for satellites missing on the E=file

      do i=1,nsat
         if( iexpnt(i).eq.0 ) then
           write(message,31) ixsat(i)
   31      format('Missing E-file elements for PRN '
     .           ,I3,' on X-file')
           call report_stat('WARNING','BCTOT','orbits/chksvs',' '
     .                     ,message,0)
         endif
       enddo


C  Locate and display  message for extra E-file satellites

      do 50 i=1,nesat
         do j=1,nsat
            if( iexpnt(j).eq.i ) goto 50
         enddo
      write(message,41) iesat(i)
   41 format('PRN ',I3,' on E-file but not on X-file')
      call report_stat('WARNING','BCTOT','orbits/chksvs',' ',message,0)     
   50 continue

      return
      end


