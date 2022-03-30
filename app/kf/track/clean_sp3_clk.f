 
CTITLE CLEAN_SP3_CLK

      subroutine clean_sp3_clk
      
      implicit none

*     Routine to clean the 99999's from the SP3 clocks by linear
*     interpolation

      include 'track_com.h'
      
* LOCAL VARIABLES

* nb, nf  - Indices of god clock before and after bad clock
* i,j     - Loop over satelllites and epochs

      integer*4 nb, nf, i, j
      
* dt      - TIme difference from back point
* dt_sp3, dclk_sp3 - Interpolation values
      
      real*8    dt, dt_sp3, dclk_sp3
      
* bood, fgood - Indicates good clock found in back and forward
*               directions
      
      logical bgood, fgood
      

***** Loop over all the satellites accross all the sp3 epochs
      do i = 1, num_sat
          do j = 1, num_sp3
             if( sp3_clk(i,j).gt.99999.d-6 ) then
             
*                This clock epoch needs interpolation
                 nb = j
                 bgood = .false.
                 do while ( .not.bgood .and. nb.gt.1 )
                     nb = nb - 1
                     if( sp3_clk(i,nb).lt.99999.d-6 ) then
                         bgood = .true.
                     end if
                 end do
                 if( .not.bgood ) then
                     write(*,110) i,j
 110                 format('WARNING: Cant find first good clock',
     .                      ' for PRN ',i3.2,' at epoch ',i4)
                 end if
                 
****             now do foward search
                 nf = j
                 fgood = .false.
                 do while ( .not.fgood .and. nf.lt.num_sp3 )
                     nf = nf + 1
                     if( sp3_clk(i,nf).lt.99999.d-6 ) then
                         fgood = .true.
                     end if
                 end do
                 if( .not.fgood ) then
                     write(*,120) i,j
 120                 format('WARNING: Cant find last good clock',
     .                      ' for PRN ',i3.2,' at epoch ',i4)
                 end if
                                    
****             Now do interpolation
                 if( bgood .and. fgood ) then
                     dt_sp3 = sp3_time(nf) - sp3_time(nb)
                     dclk_sp3 = sp3_clk(i,nf) - sp3_clk(i,nb)
                     dt = sp3_time(j) - sp3_time(nb)
                     sp3_clk(i,j) =  sp3_clk(i,nb) +  
     .                              (dclk_sp3/dt_sp3)*dt
                     write(*,200) i, j,  sp3_clk(i,j)*1.d6
 200                 format(' Interpolating clock for PRN    ',i3.2,
     .                      ' at epoch ',i4,' Value ',f10.3,' usec')
                 else if ( bgood ) then
                     sp3_clk(i,j) =  sp3_clk(i,nb)
                     write(*,210) i, j,  sp3_clk(i,j)*1.d6
 210                 format(' Back estimate of clock for PRN ',i3.2,
     .                      ' at epoch ',i4,' Value ',f10.3,' usec')
                 else if ( fgood ) then
                     sp3_clk(i,j) =  sp3_clk(i,nf)
                     write(*,220) i, j,  sp3_clk(i,j)*1.d6
 220                 format(' Forw estimate of clock for PRN ',i3.2,
     .                      ' at epoch ',i4,' Value ',f10.3,' usec')
                 else
                     sp3_clk(i,j) = 0
                 end if
              end if
          end do
      end do

      
****  Thats all
      return
      end 
 
     
                      
