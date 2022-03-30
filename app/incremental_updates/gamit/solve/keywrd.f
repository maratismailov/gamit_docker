      Subroutine KEYWRD( constraints, free_fix, phase_obs,numkey )

C     Sets keywords based on desired options, and allows for interactive
C     changes
C
C     MH Murray 881111
C
C     Add h-file mode      Dong 921229
C
C     KEYWORD (first column is "default" solution):
C        1  DEFLT                :  User code (defaults to DEFLT)
C        2  FULL  QUICK          :  Explicit or implicit bias solutions
C        3  DBLE  SINGL          :  Observable differences
C        4  LC    L1L2   L1      :  Phase observable
C        5  NOION ION    PSEUD   :  Ionospheric constraint
C        6  ATM   NOATM          :  Atmospheric constraint
C        7  FIXED FREE           :  Ambiguity resolution
C        8  STN   NOSTN          :  Station coordinate parameters
C        9  ORB   NOORB          :  Satellite orbit parameters
C        10 ZEN   NOZEN          :  Zenith delay parameters
C        11 NOCLK CLK            :  Station clock parameters
C        12 GCR GCX GLR GLX      :  H-file mode
C        13 EOP   NOEOP          :  Earth Orientation parameters
C        14 GRD   NOGRD          :  Atmospheric gradient parameters
C
C     Input:

C        islot1    :   Parameter slots (see CFMRG for documentation on correspondence)
C        free      :   Parameter estimation (1 = estimated, 0 = fixed)
c        constraints: 'tight' or 'loose'
c        free_fix  :  'free'  or 'fixd'
c        phase_obs :  'L1' 'L2' 'LC' 'L1L2' 'L12I'

C     Output:

C        nunmkey :   Number of keywords in array (set in data below)
C        keyword :   Keyword array (character*5 in common /keycom/
C
C------------------------------------------------------------------------------------

      implicit none

      include '../includes/dimpar.h'
      include 'solve.h'
      include 'parameters.h'

      integer*4 numkey,ipart,kl,isngle/0/,idble/1/

      character*2 alen
      character*4 free_fix, phase_obs
      character*5 constraints
      character*11 format

      keyword(1)='DEFLT'
C
C  Set up Variable format
C
      NUMKEY = 14
      WRITE(ALEN,'(I2)') NUMKEY
      FORMAT = '('//ALEN//'(1X,A5))'
C
C   LSQINT responses
C
      IF     (LQUICK.EQ.0) THEN
         KEYWORD(2) = 'FULL '
      ELSEIF (LQUICK.EQ.1) THEN
         KEYWORD(2) = 'QUICK'
      elseif (lquick.eq.2) then
         keyword(2) = 'IMPLI'
      ELSE
         KEYWORD(2) = 'UNK  '
      ENDIF
C
      IF     (ISNGLE.EQ.1 .OR.  IDBLE.EQ.1) THEN
         KEYWORD(3) = 'DBLE '
      ELSEIF (ISNGLE.EQ.0 .AND. IDBLE.EQ.0) THEN
         KEYWORD(3) = 'SINGL'
      ELSE
         KEYWORD(3) = 'UNK  '
      ENDIF
C
      IF ( phase_obs.eq.'L1  ' ) then
         KEYWORD(4) = 'L1   '
         KEYWORD(5) = 'NOION'
      ELSEIF (phase_obs.eq.'L2  ') then
         KEYWORD(4) = 'L2   '
         KEYWORD(5) = 'NOION'
      ELSEIF ( phase_obs.eq.'LC  ' ) then
         KEYWORD(4) = 'LC   '
         KEYWORD(5) = 'NOION'
      ELSEIF ( phase_obs.eq.'L1L2' .or. phase_obs.eq.'L12I') then
         KEYWORD(4) = 'L1L2 '
         KEYWORD(5) = 'ION  '
      ELSE
         KEYWORD(4) = 'UNK  '
         KEYWORD(5) = 'UNK  '
      ENDIF
C
C   LSQUAR responses
C
      IF     (IATCON.EQ.1) THEN
         KEYWORD(6) = 'ATM  '
      ELSEIF (IATCON.EQ.2) THEN
         KEYWORD(6) = 'NOATM'
      ELSE
         KEYWORD(6) = 'UNK  '
      ENDIF
C
C  Set Estimated parameter defaults
C
      KEYWORD(8)  = 'NOSTN'
      KEYWORD(9)  = 'NOORB'
      KEYWORD(10) = 'NOZEN'
      KEYWORD(11) = 'NOCLK'
      KEYWORD(13) = 'NOEOP'
      KEYWORD(14) = 'NOGRD'

C  Reset according to estimated parameters (see islot1 mapping in cfmrg/fills1.f)

      DO IPART=1,NTPART
         IF (FREE(IPART).EQ.1) THEN
            KL=ISLOT1(IPART)
            IF (KL.GT.0    .AND. KL.LE.300  ) KEYWORD(8) = 'STN  '
c old       IF (KL.GT.300  .AND. KL.LE.400  ) KEYWORD(10)= 'ZEN  '
            IF (KL.GT.21500.AND. KL.LE.24000) KEYWORD(10)= 'ZEN  '
            IF (KL.GT.400  .AND. KL.LE.500  ) KEYWORD(11)= 'CLK  '
            IF (KL.GT.500  .AND. KL.LE.2400 ) KEYWORD(9) = 'ORB  '
            IF (KL.GT.80000.AND. KL.LE.80006) KEYWORD(13)= 'EOP  '
            IF (KL.GT.24000.AND. KL.LE.29000) KEYWORD(14)= 'GRD  '
         ENDIF
      ENDDO

C
C     Check Biases based on modes instead of values
C
      KEYWORD(7) = 'FREE '
      IF (LQUICK.GE.1) KEYWORD(7) = 'NONE '
      IF ( free_fix.eq.'fixd' ) KEYWORD(7) = 'FIXED'
C
C     H-file mode
C
      KEYWORD(12) = ' GCR '
      if( constraints.eq.'tight' .and. free_fix.eq.'fixd' )
     .    KEYWORD(12) = ' GCX '
      if( constraints.eq.'loose' .and. free_fix.eq.'free' )
     .    KEYWORD(12) = ' GLR '
      if( constraints.eq.'loose' .and. free_fix.eq.'fixd' )
     .    KEYWORD(12) = ' GLX '

      RETURN
      END

