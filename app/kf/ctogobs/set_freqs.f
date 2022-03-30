CTITLE set_freqs

      subroutine set_freqs

      implicit none

*     This routine assigns the frequency quantities for each
*     used by autcln.  These are the same for all SVs except for GLONASS

* INCLUDES
                              
      include '../includes/kalman_param.h'
      include '../includes/cfile_def.h'      
      include 'ctogobs_com.h'  

* LOCAL 

      integer*4 i

*     Set Reference frequencies for GLONASS                   
      fL1_R0 = 1602.0d6  ! GLONASS L1 center Hz
      fL2_R0 = 1246.0d6  ! GLONASS L2 center Hz
      fL3_R0 = 1201.0d6  ! GLONASS L3 center Hz

*     For GLONASS: Set all frequenceies contrant valuw because we
*     re-map the data 

      do i=1,cf_nsat 

* MOD TAH 200511: This loop is now done while cf2 headers are read so
*       we can save values that are non-zero (issue with some GNSS)
C       fL1(i) = cf_fL1(i)
C       fL2(i) = cf_fL2(i)
* MOD TAH 180320: Save the original frequencies (in full GNSS processing
*       we map all constellations to the same frequencies.
* MOD TAH 180311: We remap GLONASS to constant frequencies so
*       set these values here.
* MOD TAH 200617: Saved the used frequencies (needed for Glonass processing)
        fL1u(i) = fL1(i)
        fL2u(i) = fL2(i)
        
        if( cf_gnss.eq.'R' ) then
           fL1(i) = fL1_R0
           fL2(i) = fL2_R0

* MOD TAH 180321: Make this optional with remap_glonass command.
*          Default is to do the mapping
* MOD TAH 200617: Only do this for Glonass.
           if( remap_glonass ) then 
              fL1u(i) = cf_fL1(i)
              fL2u(i) = cf_fL2(i)
           else  ! This option will keep ambiguities as integers rather
                 ! than scaled by original frequencies
              fL1u(i) = fL1(i)
              fL2u(i) = fL2(i)
           endif
        endif

        dfsf(i) = (fL1(i)-fL2(i))/(fL1(i)+fL2(i)) 
        sfdf(i) = (fL1(i)+fL2(i))/(fL1(i)-fL2(i)) 
        lcf1(i) = 1.d0/(1.d0 - (fL2(i)/fL1(i))**2) 
        lcf2(i) = -(fL2(i)/fL1(i))/(1.d0 - (fL2(i)/fL1(i))**2) 
        lgf1(i) = -fL2(i)/fL1(i)
        lgf2(i) = 1.d0 
        exf1(i) = 1.d0    
        exf2(i) = -fL1(i)/fL2(i) 
        pcf1(i) =  fL1(i)**2/(fL1(i)**2-fL2(i)**2) 
        pcf2(i) = -fL2(i)**2/(fL1(i)**2-fL2(i)**2) 
*         l1tecu - Change in range (m) per TECU at L1
*         l2tecu - Change in range (m) per TECU at L2
        l1tecu(i) = 40.3d0/fL1(i)**2*1.d16 
        l2tecu(i) = 40.3d0/fL2(i)**2*1.d16 

*****   Compute the multiplier factors for L1 ion delay in cycles
*       from L1 L2 phase values
        lif1(i) =   1.0d0/(1.0d0 - (fL1(i)/fL2(i))**2) 
        lif2(i) = -(fL1(i)/fL2(i))/(1.d0 - (fL1(i)/fL2(i))**2)
      
      enddo      

*     Need a scalar for clock routines
      if( cf_gnss.eq.'R' ) then        
c       don't have the Glonass center frequency or the dF on the c-file, so this is a kluge
        fClk = fL1_R0  
      else
        fClk = fL1(1)
      endif 

      return
      end   

                          



