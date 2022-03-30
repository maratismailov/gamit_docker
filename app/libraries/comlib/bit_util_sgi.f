*	IRIX OS version.  Same as vlbi/HP1000/libhp1000/bit_util.f at CFA.
*	Routines have not been tested by me.  
*	Lada L. Dimitrova      June 30  1998
*
C@KBIT   .******** TRUE IF BIT ON ***********************************.
C
      LOGICAL FUNCTION KBIT(IARRAY,IBIT) 

      implicit none

C	
C     KBIT
C
C 1.  KBIT PROGRAM SPECIFICATION
C
C 1.1.   KBIT is true if the IBIT-th bit in IARRAY is set (1),
C        false if it is not (0).
C        KBIT is designed to complement SBIT which sets or resets
C        bits identified in the same way.
C        Bits are numbered starting with 1 as the lowest-order bit
C        in the first word of IARRAY, 16 is the sign bit in the first
C        word of IARRAY, 17 is the low-order bit in the second word
C        of IARRAY, etc.
C
C 1.2.   RESTRICTIONS - NONE
C
C 1.3.   REFERENCES - FILE "SOLV2 = SOLV2 PURPOSE AND OVERALL STRUCTURE
C
C 2.  KBIT INTERFACE
C
C 2.1.   CALLING SEQUENCE: CALL KBIT(IARRAY,IBIT)
C
C     INPUT VARIABLES:
C
C     DIMENSION IARRAY(1)
C     - may or may not be an array in calling program.
C     - Variable in which the flag bits are located.
C
C     IBIT = index of bit to test. Bits are numbered starting with
C            1 as the lowest order bit in IARRAY(1), 16 is the sign
C            bit in IARRAY(1), 17 is the lowest order bit in IARRAY(2),
C            etc.

      integer*4 iarray(*), ibit
C
C
C     OUTPUT VARIABLES:
C
C     KBIT = FUNCTION VALUE = TRUE if the indicated bit is 1
C                           = FALSE if the indicated bit is 0
C
C
C     CALLING SEQUENCE VARIABLE SPECIFICATION STATEMENTS:
C
C
C     IA - array element with desired bit
C     IB - Bit number with in that element referred to 0 as least
C          significant bit

      integer*4 IA, IB

C     BIT - Fortran (Sun) function to return state of bit
C     BTEST - Fortran (HP-UX) logical function to return the state of bit

C     Logical BIT
C
C     PROGRAM STRUCTURE
C
C     1. Decompose IBIT into an array index IA and a bit index IB.
C
      IA = (IBIT+31)/32
*     Remove extra 1 from bit count since function Bit starts with
*     bit 0
      IB = IBIT - (IA-1)*32 - 1
C
C
C     2. Test the appropriate bit
C
C     KBIT =   BIT(IB,IARRAY(IA)) 
      KBIT = btest(IARRAY(IA),IB) 
C
C
C     3. Return to calling program.
C
      END
 
 
C@SBIT   .****** SET N-TH BIT IN WORD OR ARRAY **********************.
C
      SUBROUTINE SBIT(IARRAY,IBIT,IVALUE)
C
C     SBIT
C
C 1.  SBIT PROGRAM SPECIFICATION
C
C 1.1.   SBIT ... sets a bit in IARRAY indexed by IBIT to IVALUE (0 or 1).
C        Logical function KBIT can then test the bit (true if bit is 1).
C        The indexing scheme is described below under INPUT VARIABLES.
C
C 1.2.   RESTRICTIONS - IF IVALUE IS NOT 0 IT IS ASSUMED TO BE 1.
C                     - IF IVALUE IS OMMITTED IN CALLING SEQUENCE, 1 IS ASSUMED.
C
C 1.3.   REFERENCES - FILE "SOLV2 = SOLV2 PURPOSE AND OVERALL STRUCTURE
C
C 2.  SBIT INTERFACE
C
C 2.1.   CALLING SEQUENCE: CALL SBIT(IARRAY,IBIT,IVALUE)
C        NOTE THAT IVALUE IS OPTIONAL ARG, IF OMITTED DEFAULTS TO 1.
C
C     INPUT VARIABLES:
C
C     -Array containing the bits to be set or reset.
C     - may or may not really be dimensioned in calling program.
C
C     IBIT = index of bit to be set or reset.
C            Bit 1 is the lowest order bit of IARRAY(1), bit 16 is the sign
C            bit of IARRAY(1), bit 17 is the lowest order bit of IARRAY(2),
C            etc. This is the same indexing scheme employed in func. KBIT.
C
C     IVALUE = 0 or 1 = value to which bit is to be set;
C            = 1 if IVALUE is ommitted from calling sequence.
C            IF NON-ZERO, 1 IS ASSUMED.

      integer*4 iarray(1), ibit, ivalue
C
C
C     OUTPUT VARIABLES:
C
C
C     CALLING SEQUENCE VARIABLE SPECIFICATION STATEMENTS:
C
C
C 3.  LOCAL VARIABLES
C
C     IA - Array element with bit to be set
C     IB - Bit number with in element.  Counts from zero

      integer*4 IA, IB
C
C END EQUIVALENCE
C
C     PROGRAM STRUCTURE
C
C
C     1. Decompose IBIT into an array index IA and a bit index IB.
C
      IA = (IBIT+31)/32
      IB = IBIT - (IA-1)*32 - 1
C
C
C     2. IF IVALUE is 1, force the appropriate bit to 1.
C        BIC - Sun fortran function to clear bit
C        BIS - Sun fortran function to set bit
C        IBSET - HP-UX fortran function to set bit
C        IBCLR - HP-UX fortran function to clear bit
 
      if( ivalue.eq.1) then
C         call bis(ib, iarray(ia))             
          iarray(ia) = ibset(iarray(ia),ib)
      else
C         call bic(ib, iarray(ia))
          iarray(ia) = ibclr(iarray(ia),ib)
      end if
C
C
C     4. Return to calling program.
C
      END
 
