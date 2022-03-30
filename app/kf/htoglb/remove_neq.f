CTITLE REMOVE_NEQ

      subroutine remove_neq(np, cov_parm, sol_parm)
 
      implicit none

*     Routine to implicitly solve for parameters to be
*     removed from the Normal equations and the make the
*     modiications needed the estimated parameters.

      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'
 
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'

* PASSED
  
      integer*4 np   ! Number of parameters

* cov_parm(np,np) - Covarince matrix
* sol_parm(np)    - Solution vector
 
      real*8 cov_parm(np,np), sol_parm(np)
 
* LOCAL
      integer*4 i,j,k   ! Loop counters
      integer*4 P   ! Current parameter being tested
      integer*4 p1, p2  ! Parameter numbers for removing PCO constraint

***   MATLAB CODE That is implemented..
*     %% More general code: Here work on upper triangle using lower triangle.
*     P = 3;
*     for k = 1:NP
*         if k ~= P 
*             BVec(k) = BVec(k) - BVec(P)*NEQ(k,P)/NEQ(P,P);
*             for j = k:NP
*                 if j ~= P
*                     NEQ(j,k) = NEQ(j,k) - NEQ(j,P)*NEQ(P,k)/NEQ(P,P);
*                 end
*             end
*         end
*     end
*     for j = 1:NP
*         NEQ(P,j) = 0.0 ;
*         for k = 1:NP
*             NEQ(j,k) = NEQ(k,j);
*         end
*     end
*     BVec(P) = 0;

*     Loop over parameters to be eliminated (update upper-triangle
*     using values from lower triangle and then copy to lower).

      do P = 1,np
         if( rm_param(P) ) then    ! Pre-eliminate this parameter
            do k = 1,np
               if( k.ne.P ) then
                   sol_parm(k) = sol_parm(k) - sol_parm(P)*
     .                           cov_parm(k,P)/cov_parm(P,P)
                   do j = k,np
                      if( j.ne.P ) then
                          cov_parm(j,k) = cov_parm(j,k) - 
     .                        cov_parm(j,P)*cov_parm(P,k)/cov_parm(P,P)
                      endif
                   enddo
               endif
             enddo
*            Now clear the row and column eliminated copy upper triangle
*            lower triangle
             do j = 1,np
                cov_parm(P,j) = 0.0
                do k = 1,np
                   cov_parm(j,k) =  cov_parm(k,j)
                end do
             end do
             sol_parm(P) = 0.0
             cov_parm(P,P) = 1.d-6
         end if
      end do

****  Now see if we are to remove the satellites antenna constaints
*     if applied
* MOD TAH 210503: Constrants now handled in normal equation ADD_NEQ_CON
*     routine.
C     if( index(decon_str,'A').gt. 0 ) then
C        do i = 1,qnum_svs
C        do k = 1,3
C           p1 = qparn_svs(max_svs_elem-3+k,i)
C           if( p1.gt.0 ) then
C              do j = 1, np
C                 if( atoo(j).eq.p1 ) p2 = j
C              end do
C           end if
C           if ( p1.gt.0 .and. p2.gt.0 ) then
C               write(*,740) i, k, p1,p2, cov_parm(p1,p1), 
C    .                       qapr_sig(p2)
C740            format('Removing Sat PC Pos ',i3,i4,' NP ',2i5,
C    .              ' NEQ and Var ',2e15.6)
C               if( cov_parm(p1,p1).gt.qapr_sig(p2) )
C    .          cov_parm(p1,p1)=cov_parm(p1,p1)-qapr_sig(p2)
C           end if
C        end do
C        end do
C     end if

      RETURN
      end


