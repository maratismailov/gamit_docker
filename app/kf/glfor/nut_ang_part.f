CTITLE NUT_ANG_PART
 
      subroutine nut_ang_part( source, deriv, pos, comp, epoch)

      implicit none  
 
*     Routine to compute nutation angle partial derivatives
 
      include '../includes/const_param.h'
 
*   comp    - Component of position (1=RA, 2=dec)
*   source  - Global source number
 
      integer*4 comp, source
 
*   cosep, sinep- Cosine and sin of mean obliquity of ellciptic
*   deriv(1)    - Derivative with respect to DEPS and DPSI
*   epoch       - Epoch of experiment (JD)
*   pos(1)      - position of source (mas)
 
      real*8 sinep, deriv(2), epoch, pos(2) 
c     real*8 cosep  - not used
 
 
      data  sinep / 0.39777716d0 /
c     data  cosep / 0.9174820603d0 /  - not used
 
****  Use GOTO for RA or dec
 
      goto ( 100, 200 ) comp
 
*     RA partial
  100 continue
          deriv(1) = sinep*sin(pos(1)/rad_to_mas)*
     .                     tan(pos(2)/rad_to_mas)
 
          deriv(2) =      -cos(pos(1)/rad_to_mas)*
     .                     tan(pos(2)/rad_to_mas)
 
          return
 
*     Dec partial
  200 continue
          deriv(1) = sinep*cos(pos(1)/rad_to_mas)
          deriv(2) =       sin(pos(1)/rad_to_mas)
          return
 
****  Thats all
      end
 
