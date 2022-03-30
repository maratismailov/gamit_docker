      SUBROUTINE APRIOR(ISITE,NSAT,ISAT,NORB,
     $                  SITE,PREVAL,ITSAT,TOTSAT)
C   
C PURPOSE:
  
C     SELECT SATELLITES
C  
      implicit none
C
      include '../includes/dimpar.h'
C                  
      integer i,j,k,kk
      INTEGER ISITE
      INTEGER IELM
      INTEGER NORB 
      INTEGER ISAT(MAXSAT),TOTSAT(MAXSAT)
      INTEGER NSAT,ITSAT
      REAL*8  PREVAL(MAXPRM)
      REAL*8  APRCRD,APRCLK,APRORB,APRVAL,ADJUST,apreop,aprzen,aprgrd
     .     ,  aprant
      COMMON/MERCOM/APRVAL(MAXPRM),ADJUST(MAXPRM),APRCRD(maxsit,3),
     $              APRCLK(maxsit,3),APRORB(MAXORB,MAXSAT),
     $              aprant(3,maxsat),APREOP(6),APRZEN(maxsit),
     $              aprgrd(maxsit,2)
      CHARACTER*12 SITE

      LOGICAL oldslot 
c      logical piksat
c      DATA PIKSAT/.FALSE./

C.....APRIORI SITE COORDINATES (FIRST APPEARANCE ONLY)

         DO  I=1,3
            APRCRD(ISITE,I)=PREVAL(I)
         enddo
                     

C......A priori atmospheric zenith delays

         oldslot = .false.
         if( preval(7).ge.2.0d0 .and. preval(7).le.2.5d0 ) then
            oldslot = .true.
            aprzen(isite) = preval(7) 
          else
            aprzen(isite) = preval(4) 
          endif    

C.....A priori atmospheric gradients

      do i=1,2    
           aprgrd(isite,i) = 0.d0  
      enddo   
 

c.....A priori station clock epoch (rate and acceleration no longer estimatable)

      aprclk(isite,1) = preval(5)
                                  
C        
c..... A priori orbital elements and non-gravitational force parameters

         do i=1,itsat
            do j=1,nsat
              if ( isat(j).eq.totsat(i) ) then
                 k=0
                 do ielm=1,norb
                    k=k+1
                    aprorb(ielm,i) = preval(5+(j-1)*norb+k)
                 enddo
              endif
            enddo
          enddo


c..... A priori satellite antenna offsets
                       
          kk = 0      
          do i=1,itsat
             do j=1,nsat
                if ( isat(j).eq.totsat(i) ) then
                   do k=1,3
                     kk=kk+1
                     aprant(k,i) = preval(5+nsat*norb+kk)
                   enddo
                endif
             enddo
          enddo  

C.... A priori Earth orientation parameters
                          
          do  j=1,6
            apreop(j)= preval(5+(nsat*(norb+3))+j)
          enddo

      RETURN
      END
