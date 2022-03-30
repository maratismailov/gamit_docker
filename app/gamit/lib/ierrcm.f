      integer function ierrcm(in1,in2)
      integer ierr1,ierr2,in1,in2

      include '../includes/errflg.h'

c     Mod:  implement new error flag conventions -- MHM 880331
c
c*****************************************************************************
c Function IERRCM calculates the error flag of a combined observation based  *
c upon the error flags of the constituent observations according to the      *
c following scheme with:                                                     *

c           iggood    0      good observation                                *
C           ignone    1      no observation                                  *
c           igchop    2      bad, deleted observation                        *
c           igisok   98      decided good?                                   *
c           igunwt   -1      unweighted observation                          *
c           igrewt  -11      reweighted observation                          *
c           iglamp    3      low amplitude (potentially ok)                  *
c           igloel    4      low elevation, (potentially ok)                 *
c           ig2few    5      too few points to tell (potentially ok)         *
c           igbias           insert a new bias parameter here (ok)           *
c           igoutl           probably an outlier (prob. bad)                 *


c     The +200 codes are for internal use only.
c     They must be unmapped when writing the C-file.

      if (in1 .gt. 100) then
         ierr1 = in1 - 200
      else
         ierr1 = in1
      endif
      if (in2 .gt. 100) then
         ierr2 = in2 - 200
      else
         ierr2 = in2
      endif


c     default should be no data
      ierrcm = ignone + 200

      if(ierr1 .eq. iggood .or. ierr2 .eq. iggood) ierrcm = iggood + 200
      if(ierr1 .eq. igisok .or. ierr2 .eq. igisok) ierrcm = igisok + 200
      if(ierr1 .eq. igbias .or. ierr2 .eq. igbias) ierrcm = igbias + 200
      if(ierr1 .eq. igoutl .or. ierr2 .eq. igoutl) ierrcm = igoutl + 200
      if(ierr1 .eq. igrewt .or. ierr2 .eq. igrewt) ierrcm = igrewt + 200
      if(ierr1 .eq. igloel .or. ierr2 .eq. igloel) ierrcm = igloel + 200
      if(ierr1 .eq. igunwt .or. ierr2 .eq. igunwt) ierrcm = igunwt + 200
      if(ierr1 .eq. iglamp .or. ierr2 .eq. iglamp) ierrcm = iglamp + 200
cold  if(ierr1 .eq. igloel .or. ierr2 .eq. igloel) ierrcm = igloel + 200
      if(ierr1 .eq. ig2few .or. ierr2 .eq. ig2few) ierrcm = ig2few + 200
      if(ierr1 .eq. igchop .or. ierr2 .eq. igchop) ierrcm = igchop + 200
      if(ierr1 .eq. ignone .or. ierr2 .eq. ignone) ierrcm = ignone + 200
C
      return
      end
