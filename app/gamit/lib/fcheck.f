      logical function FCHECK (fname)
C
C     CHECK FOR REQUIRED INPUT FILES

      character*(*) fname
      logical okay
C

      inquire (file  = fname,
     .         exist = okay)

      fcheck = okay

      RETURN
      END
