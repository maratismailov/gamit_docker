      character*4 function errcod(ierror)

c     return a 4-letter code describing the error flag


      integer ierror

      include '../includes/errflg.h'

      errcod = '????'
      if (ierror .eq. iggood) errcod = 'good'
      if (ierror .eq. ignone) errcod = 'none'
      if (ierror .eq. igchop) errcod = 'chop'
      if (ierror .eq. igisok) errcod = 'isok'
      if (ierror .eq. igunwt) errcod = 'unwt'
      if (ierror .eq. igrewt) errcod = 'rewt'
      if (ierror .eq. iglamp) errcod = 'lamp'
      if (ierror .eq. igloel) errcod = 'loel'
      if (ierror .eq. ig2few) errcod = '2few'
      if (ierror .eq. igbias) errcod = 'bias'
      if (ierror .eq. igoutl) errcod = 'outl'
c      this not used (undefined in errflg.h)
c      if (ierror .eq. ighdwr) errcod = 'hdwr'


      return
      end

