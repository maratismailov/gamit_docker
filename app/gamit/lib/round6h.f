      subroutine round6h( day )

c     round a decimal day to the nearest 6hr interval
c     used for VMF and ATML grid files
                    
      real*4 day,day4

      day4 = nint(4*day)
      day = day4/4.
      return
      end

