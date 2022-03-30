      subroutine headers(ih,isnx,in_h,idat,istn,use_sol)

c  Purpose: To read the hfile header and write the sinex header
c
c  by Paul Tregoning
c  14th July 1995
c
c  IN:
c        ih      - unit number of hfile                       I*4
c        isnx    - unit number of sinex file                  I*4
c        in_h    - hfile name                                 C*10
c        idat    - unit number of sinex.dat file              I*4
c        istn    - unit number of station.info                I*4
c        use_sol - soln to read from hfile                    C*3

      implicit none

      include '../includes/dimpar.h'

      integer ih,isnx,idat,istn,i,qnum_sites,qnum_sats,prns(maxsat)
     .        ,nzen(maxsit),num_par,num_orb,qnum_par
     .        ,qparn_sites(3,maxsit),qparn_svs(maxorb,maxsat)
     .        ,qparn_eop(2,3),qnum_eop
     .        ,est_type(maxsit*3+maxsat*maxorb+6)

      character*3 use_sol
      character*4 aphs_mod(maxsit),qsite_code(maxsit),sol_type
      character*10 in_h,cfiles(maxsit)
      character*12 times(6)
      character*16 rcvr(maxsit,3),antenna(maxsit,2)
      character*20 qsite_names(maxsit)
      character*48 solv_vers
      character*50 agency

      real*8 elev_cutoff(maxsit),ant_off(maxsit,3),phs_off(maxsit,6)
     .      ,apr(maxsit*3+maxsat*maxorb+6),adj(maxsit*3+maxsat*maxorb+6)
     .      ,cov_parm(maxsit*3+maxsat*maxorb+6,maxsit*3+maxsat*maxorb+6)
     .      ,cov_apr(maxsit*3+maxsat*maxorb+6,maxsit*3+maxsat*maxorb+6)
     .      ,gps_utc,fix_eop(6),constr(maxsit*3+maxsat*maxorb+6)

      do i=1,6
        fix_eop(i) = 0.d0
      enddo

c  read the hfile header and return the information in arrays - sinex needs it in
c  different orders to what the hfile currently has
      call rdhhed(ih,idat,times,qnum_sites,qnum_sats,qsite_code
     .            ,qsite_names,cfiles,rcvr,antenna,aphs_mod,elev_cutoff
     .            ,nzen,ant_off,phs_off,prns,num_par,solv_vers,gps_utc
     .            ,sol_type,use_sol)

c  now read in the site/sat apriori coordinates and adjustments
      call get_estimates(ih,qnum_sites,qnum_sats,apr,adj,num_orb
     .                   ,qnum_par,qparn_sites,qparn_svs,qparn_eop
     .                   ,qnum_eop,gps_utc,fix_eop)

c  set the flag for whether eop parameters have been estimated
      if(qnum_eop.gt.0)then
        sol_type(3:3) = 'E'
      else
        sol_type(3:3) = ' '
      endif



c PT950815: Kluge to temporarily exclude the sat params from the total number of parameters
      num_par = qnum_sites*3 + qnum_eop

c  begin writing the sinex file header
      call wrtsnxhd(isnx,in_h,times,num_par,idat,solv_vers,sol_type
     .              ,agency)

c  write the input file information to the sinex file
      call write_input(isnx,in_h,times,num_par,sol_type,agency)

c  now create the SITE/ID block
      call site_id(isnx,idat,qnum_sites,qsite_code,qsite_names,apr )

c now create the SITE/RECEIVER and SITE/ANTENNA blocks
      call wrtsit_info(isnx,istn,qnum_sites,qsite_code,rcvr,antenna
     .                 ,times(5) )

c now create the SITE/GPS_PHASE_CENTER and SITE/ECCENTRICITY blocks
      call wrtant_info(isnx,qnum_sites,antenna,ant_off,phs_off
     .                 ,aphs_mod,qsite_code,times(2),times(3) )

c now create the SOLUTION/EPOCHS block
      call wrtepo_info(isnx,qnum_sites,qsite_code,times)

c now create the SOLUTION/ESTIMATE and SOLUTION/MATRIX blocks
      call soln_info(ih,isnx,qnum_sites,qnum_sats,qsite_code
     .               ,apr,adj,est_type,times(6),num_orb,qnum_par
     .               ,qparn_sites,qparn_svs,qparn_eop,cov_parm
     .               ,qnum_eop,fix_eop,constr,cov_apr,sol_type)


c and finally close the sinex file
      write(isnx,'(a)')'%ENDSNX'
      close(ih)
      close(isnx)

      return
      end
