       function trapez_intg(id,restart,integrand,step)
c
c** Function to perform one step trapezoidal integration. When 'restart'
c   is true it starts from zero. If 'restart' is not true it will assume
c   that the integration starts at (integrand=) zero. This integrator is
c   designed to handle calls from 'ncalls' different subroutines, each
c   requesting a different integral. The function accumulate each integral
c   based on the 'id' of the caller. So it is crucial that every caller keep
c   its 'id' consistently.
c
c   Yoaz Bar-sever. April 1994.
c
c   yeb@cobra.jpl.nasa.gov              MS 238-600
c   voice: (818) 354-2665               Jet Propulsion Laboratory
c   FAX:   (818) 393-4965               4800 Oak Grove Dr.
c                                       Pasadena, California  91109
c
c**  Copyright (C) 1994, California Institute of Technology
c
       implicit none
c
c** Input variables
       integer id             ! the id number of the caller subroutine.
                              ! if you want to perform several integrations
                              ! at the same time, give each integration a
                              ! separate id number. No more than 'ncalls'
                              ! integrations are allowed.
       logical restart        ! if true assumes the beginning of the integral
       double precision integrand  ! the new value of the integrand
       double precision step       ! the step size
c
c** Output variables
       double precision trapez_intg ! the updates value of the integral
c
c** Local variables
       integer ncalls          ! maximum number of different integrations that
       parameter (ncalls=9)    ! this subroutine can accumulate concurrently
       double precision last_integral(ncalls) ! last value of the integral
       data last_integral/ncalls*0.0d0/
       double precision last_integrand(ncalls) ! last value of the integrand
       data last_integrand/ncalls*0.0d0/
c
       save last_integral,last_integrand
c
c** Check if id is within [1,ncalls].
       if (id .lt. 1 .or. id .gt. ncalls) then
          write(6,*) ' '
          write(6,*) 'Error in function trapez_intg!'
          write(6,*) 'The id number is ilegal.'
          write(6,*) 'It must be an integer between 1 and',ncalls
          write(6,*)'If you need to do more than',ncalls,' simultaneous'
          write(6,*) 'integrations - contact cognizant engineer.'
          write(6,*) 'Program stopped.'
          stop
       endif
c
c** Intialize if necessary.
       if (restart) then
          last_integral(id)=0.0d0
          last_integrand(id)=integrand
          trapez_intg=0.0d0
          return
       endif
c
c** Update the value of the integral.
       trapez_intg=last_integral(id)+
     &             0.5d0*(last_integrand(id)+integrand)*step
c
c** Save the new values
       last_integral(id)=trapez_intg
       last_integrand(id)=integrand
c
c** Bye.
       return
       end
