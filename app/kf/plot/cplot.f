      Program plot

      implicit none 
c
c
c     Program to produce publication quality plots using the HP
c     GRAPHICS 1000/II plotting package.
c
c     The program can be executed in batch mode by passing all
c     of the file names needed through the runstring
c     :RU,plot,<input control file>,<plot LU>,<file to be plotted>,<headers>
c     where
c     <input control file> is the name of file from which
c        commands should be read.  May be a terminal lu number
c     <plot lu> is the name of metafile to be created by plot
c     <file to be plotted> is the name of the data file
c     <headers> is the number of header records in the file
*     <skip colmn 1> if value is non-zero then a character in column
*        one will be treated as commented.
c
c     The control file format is basically free format with the
c     following rules:
c      - a nonblank character in column one denotes a comment line
c      - each command should be padded so that it is at least
c        8 characters long.
c      - The available commands are given in &PLOBD
c
c     Th format of the input data file is determined by commands
c     in the control file.
 
c Include files
c -------------
*                         ! the PLOT parameter file
      include 'plot_param.h'
c
*                         ! the PLOT_COM common block
      include 'plot_com.h'
c
*                         ! the ema declaration
      include 'plot_ema.h' 

c      external plot_bd
 
c Variable declarations
c ---------------------
c
c User buffer for reading comands
c ------------------------------
c buffer -- user buffer
c ierr  -- an error flags
c i     -- loop counter
c
 
      integer*4 ierr, i 
 
c
c
c Functions used
c -------------- 
c rcpar   -- HP runstring utility 
c lgbuf   -- HP utilty to increase buffer size
c loglu   -- return user lu
c id      -- a dummy argument
c
 
      integer*4 loglu, id
 
c
c
c Scratch blank common area
c -------------------------
c scr_common -- an dummy array used for space
c
 
      integer*4 scr_common(scr_size)
 
c
 
      common scr_common
 
c
 
      real*4 rlist(2)
 
c
 
      integer*4 ilist(5)
 
c
c................................
c Start execution of program    .
c................................
c
c.... Get the runstring.  Use the FTN77 routine RCPAR which will
c     return the runstring parameters in character strings
c     Initialize the runstrings incase their was nothing passed
c
      do i = 1,5
         runstring(i) = '  '
      end do
c
*                             ! get the user terminal
      termlu = loglu(id)
      pushlu = 0              ! Set default
c
c.... Now get the runstring
      do i = 1,5
         call rcpar(i,runstring(i))
      end do
c
c     Decode the runstring
c     ====================
 
      call PLTS8
c
c     Set up the plot lu
c     ==================
c
*                     ! just initalize
      pel = -1
*                                             ! initialize with file read
      if( input_file(1:1).ne.' ' ) pel = -3
c
      call PLTS1
c
c.... See program to abort if error
      call report_error('SEGLD',ierr,'schedul','PLTS1',0,'plot')
c
c.... Enquire about the work station
      call jiws(1,254,0,2,ilist,rlist)
c
c
c     Now start reading commands from userlu
c     ======================================
c
c     Schedule the next segment
c
      pel = 0
      do while (pel.ne.1)
c
         call PLTS2
c
c....    See if we need PLTS1 segment
         if( pcontrol.eq.1 ) then
            call PLTS1
         end if
c
c....    See if we need PLTS3 segment
*                                         ! AXES and LABELS
         if( pcontrol.eq.3 ) then
            call PLTS3
         end if
c
*        See if pltsp needed (Polynomial fitting segment)
*                                         ! Polynomial fitting
         if( pcontrol.eq.4 ) then
            call PLTSP
         end if
 
c....    See if we need PLTS9 segment
*                                     ! HELP
         if( pcontrol.eq.9 ) then
            call PLTS9
         end if
 
c....    See if any labels have to output
*                                       ! Yes there are
         if( num_labels.gt.0 ) then
            call PLTSA
         end if
c
      END DO
c
c.... Thats all
      close(userlu)
c
c.... End the pltting package
      pel = -2
      call PLTS1

      end

*     Include the block data in the main program since some compilers/linkers (specifically OSX) will not correctly 
*     load the initialized block data routines from library archives.
      include 'plot_bd.f' 
