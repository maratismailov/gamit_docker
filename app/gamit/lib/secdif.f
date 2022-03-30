      real*8 function secdif (iwkn2,sow2,iwkn1,sow1)
c
c     return difference in seconds between epoch 2 minus epoch 1
c     handles going over weeks and positive and negative diffs
c
C     GPS week numbers
      integer*4   iwkn1,iwkn2
C     GPS second of week
      real*8      sow1,sow2

C     seconds per week
      real*8      spw,
     .            hold1,hold2
      integer*4   ihold

      data        spw/604800.d0/


      ihold = min(iwkn1,iwkn2)

      hold2 = spw * dble(iwkn2 - ihold) + sow2
      hold1 = spw * dble(iwkn1 - ihold) + sow1

      secdif = hold2 - hold1

      return
      end
