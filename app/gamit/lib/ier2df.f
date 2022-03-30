      Subroutine ier2df( cf_ierfl, data_flag )

c     This routine will map the cfile error flag to the appriopriate
c     bits in data_flag (see Gobs_data.h) for definition.
c     Written by T. Herring for autcln; transferred to GAMIT /lib with
c     change of ier to I*4 and use of errflg.fti by R. King May 95.


c INCLUDES
c   erinclude '../includes/errflg.h'
      include '../includes/errflg.h'

c PASSED VARIABLES

c   cf_ierfl   - Original GAMIT C-file error flag  - numerical corrspondences in /lib/errflg.fti

c             a good observation
c              iggood
c             deleted observation
c              igchop
c             no observation
c              ignone
c             low elevations (potentially OK)
c              igloel
c             low amplitude (potentially OK)
c              iglamp
c             unweight
c              igunwt
c             outlier?
c              igoutl
c             not enough points for detection
c              ig2few
c             really OK
c              igisok
c             reweight
c              igrewt
c             add a bias here
c              igbias


c   data_flag   - AUTCLN Gobs_file data flag (bit mapped)


      integer*4 cf_ierfl,data_flag

c LOCAL VARIABLES
c     None

c     Test each of the possible cases
      data_flag = 0

      if( cf_ierfl.ne.0 ) then
c                                         ! Unweighted
          if( cf_ierfl.eq. igunwt ) then
             call sbit(data_flag, 14, 1)
c                                         ! Deleted data
          else if( cf_ierfl.eq.  igchop ) then
              call sbit(data_flag, 4,  1)
c                                         ! Low SNR
          else if( cf_ierfl.eq.  iglamp ) then
              call sbit(data_flag,  2, 1)
c                                         ! Low elev
          else if( cf_ierfl.eq.  igloel ) then
              call sbit(data_flag, 15, 1)
c                                         ! Bias flag
          else if( cf_ierfl.eq. igbias ) then
              call sbit(data_flag, 32, 1)
c                                         ! Unknown Cfile flag
          else
              call sbit(data_flag, 21, 1)
          end if
      end if

      return
      end

