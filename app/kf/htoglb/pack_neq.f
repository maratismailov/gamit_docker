CTITLE PACK_NEQ

      subroutine pack_neq(np, cov_parm, sol_parm)
 
      implicit none

*     Routine to remove the rows and colums from the normal
*     equations and remap the stations for stations not used.
*     All needed to handle 3-day COD SINEX files.

      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'
 
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'

* PASSED
  
      integer*4 np   ! Number of parameters

* cov_parm(np,np+1) - Covarince matrix (add extra column 
*                   because sol_parm is contiguous with
*                   cov_parm
* sol_parm(np)    - Solution vector
 
      real*8 cov_parm(np,np+1), sol_parm(np)
 
* LOCAL
      integer*4 i,j,k   ! Loop counters
      integer*4 c,r     ! Column/row counter (pack down columns)
      integer*4 iel     ! Remapped row/col of packed cov_patm (due
                        ! change in dimension need to address as column)
      integer*4 new_np  ! New number of parameters
      integer*4 new_qnum_sites  ! New number of sites
      integer*4 oton(max_glb_parn)  ! Mapping of old parameter number
                 ! to new number (values always reduce i.e. 3->1, 4->2
     .,         ntoo(max_glb_parn)  ! Reverse
      integer*4 ston(max_glb_sites) ! Mapping of old site number
                 ! to new number.
     .,         ntos(max_glb_sites) ! Reverse
      integer*4 pw, sw  ! Parameter and site counters.

****  This code does not work yet.  Complicated by cov_parm/sol_parm
*     being contigous.  Ultimate issue seems to be getting things
*     fully lined up for output.

      RETURN

****  First get the mapping from the current parameter 
*     number to the new parameter.
      pw = 0
      do i = 1, np
         if( .not. rm_param(i) ) then  ! Keeping this one
            pw = pw + 1
            oton(pw) = i    ! Old number to new
            ntoo(i)  = pw 
         endif
      enddo
      new_np = pw   ! New number of parameters

***   Now do same for sites
      sw = 0
      do i = 1, qnum_sites 
         if( .not. rm_site(i) ) then  ! Keeping this one
            sw = sw + 1
            ston(sw) = i    ! Old site to new number
            ntos(i)  = sw
         endif
      enddo
      new_qnum_sites = sw   ! New number of parameters

      write(*,120) np, new_np, qnum_sites, new_qnum_sites
 120  format('Reducing parameters from ',i4,' to ',i4,/,
     .       'Reducing sites      from ',i4,' to ',i4)


***** Start with normal equations and solution vector.  Tricky
*     bit is solution vector is contiguous with cov_parm and is
*     np+1 column.  So comparct cov_parm first
      do c = 1, new_np  ! Move accross column
         do r = 1, new_np   ! Run down row in columns (avoids 
                            ! overwriting
*           Compute the colum position of new value for
*           new_np dimension.
            iel = (c-1)*new_np+r
            cov_parm(iel,1) = cov_parm(oton(r),oton(c))
         enddo
      enddo

*     Now move sol_parm down to new_np+1 column (note: spl_parm
*     is total value (aproiri added at end of inv_normeq)
      do i = 1, new_np
         sol_parm(i) = sol_parm(oton(i))
         cov_parm(i,new_np+1) = sol_parm(oton(i))
      end do

****  Now compact the other arrays that are aligned to parameters.
      do i = 1, new_np
         qscale(i)    = qscale(oton(i))
         qref_ep(i)   = qref_ep(oton(i))
         qapr(i)      = qapr(oton(i))
         qapr_sig(i)  = qapr_sig(oton(i)) 
      end do

****  Remap qparn_sites (Only do sites positions (only in CODE sinex)
      do i = 1, new_qnum_sites
         do j = 1, 3
            qparn_sites(j,i) = ntoo(qparn_sites(j,ston(i)))
         enddo
      enddo

****  Now map PMU (only one will be left; center day)
      do i = 1,2      !  Value and rate
         do j = 1,3   !  X,Y and UT1
            qparn_mul_pmu(i,j,1) = ntoo(qparn_mul_pmu(i,j,2))
            qnum_mul_pmu(i,j) = 1
         enddo
      end do

****  Now map Translation
      do i = 1,3
         qparn_tran(i,1) = ntoo(qparn_tran(i,1))
      end do

****  Now mapping satellite antenna offset.
      do i = 1, qnum_svs
         do j = 1, 3      ! XYZ PCO values
            qparn_svs(max_svs_elem-3+j,i) = 
     .        ntoo(qparn_svs(max_svs_elem-3+j,i))
         end do
      end do 

      np = new_np
      qnum_parn  = np
      qnum_apr_codes = np


****  Now remap site names
      do i = 1, new_qnum_sites

         site_pos(:,i)  = site_pos(:,ston(i))

         qsite_names(i) = qsite_names(ston(i))
         if( ichar(qfull_names(ston(i))(1:1)).eq.0 ) then
            qfull_names(i) = qfull_names(ston(i)-1)  
         else
            qfull_names(i) = qfull_names(ston(i)) 
         endif 

         code_all(i)    = code_all(ston(i))
         pt_all(i)      = pt_all(ston(i))
         occ_all(i)     = occ_all(ston(i))

         qdata_st(i)    = qdata_st(ston(i))
         qdata_en(i)    = qdata_en(ston(i))
         qrecv_st(i)    = qrecv_st(ston(i)) 
         qrecv_en(i)    = qrecv_en(ston(i))
         qante_st(i)    = qante_st(ston(i))
         qante_en(i)    = qante_en(ston(i))
         qarp_ecc(:,i)  = qarp_ecc(:,ston(i)) 
         qL1a_ecc(:,i)  = qL1a_ecc(:,ston(i))
         qL2a_ecc(:,i)  = qL2a_ecc(:,ston(i))
         qelev_cut(i)   = qelev_cut(ston(i))
         qantdaz(i)     = qantdaz(ston(i))

         qant_mod(i)    = qant_mod(ston(i))
         qante_sn(i)    = qante_sn(ston(i))
         qrecv_sn(i)    = qrecv_sn(ston(i))
         qrecv_ty(i)    = qrecv_ty(ston(i)) 
         qante_ty(i)    = qante_ty(ston(i))
         qrecv_fw(i)    = qrecv_fw(ston(i))
         qradome_ty(i)  = qradome_ty(i)
      end do

****  Mow finally reset number of parameters and sites
      qnum_sites = new_qnum_sites
      cnum_sites = qnum_sites

***** Thats all
      RETURN
      end






         


 

      
      
