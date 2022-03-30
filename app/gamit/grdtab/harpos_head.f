      Subroutine harpos_head( luh,message )      

c     Construct the header of a HARPOS file for OTL
c     R. King from Scherneck sample.   15 January 2008

      implicit none
                       
      integer*4 luh,i

      character*256 message                             
      character*62 harlines(11) 
 
      data harlines
     . /'H  m2         2.169437D+00   1.405189027044D-04   1.240D-23'
     . ,'H  s2         6.283185D+00   1.454441043329D-04  -7.350D-40'
     . ,'H  n2         6.097067D+00   1.378796996516D-04  -1.860D-23'
     . ,'H  k2         3.506941D+00   1.458423171028D-04   2.130D-24'
     . ,'H  k1         3.324267D+00   7.292115855138D-05   1.060D-24'
     . ,'H  o1         5.128356D+00   6.759774415297D-05   1.130D-23'
     . ,'H  p1         2.958919D+00   7.252294578148D-05  -1.060D-24'
     . ,'H  q1         2.772800D+00   6.495854110023D-05  -1.970D-23'
     . ,'H  mf         4.479096D+00   5.323414398410D-06  -1.030D-23'
     . ,'H  mm         5.497148D+00   2.639203052741D-06   3.100D-23'
     . ,'H  ssa        3.653480D-01   3.982127698995D-07   2.130D-24'/

                 
        write(luh,'(a)') 'HARPOS Format version of 2002.12.12'
cx        write(luh,'(2a)') '# Ocean tidal loading values for ',dfile  
        write(luh,'(2a)') '#',message  
        write(luh,'(a)') '#'
        write(luh,'(a)') '# Values included for each station:' 
        write(luh,'(a,/,a,/,a)') '#'
     .  ,'#  Harmonic   Phase          Frequency           Acceleration'
     .  ,'#'
      do i=1,11
        write(luh,'(a)') harlines(i)
      enddo
      write(luh,'(a)') '#'
      return
      end  

