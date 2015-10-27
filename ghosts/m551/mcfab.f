*deck %W%  %G%
      subroutine mcfab(fab,nco,mrs,hess,lhcc)
C
C***Begin prologue
C***Date written       871022   (yymmdd)
C***Revision date      yymmdd   (yymmdd)
C
C***Keywords
C***Author             Lengsfield, Byron (BRL)
C***Source             %W%   %G%
C
C***Purpose
C
C***Description
C
C***References
C
C***Routines called    (none)
C
C***End prologue
C
      implicit real*8(a-h,o-z)
cc
cmp   extended dummy fab,hess
cc
      common / number / zero,pt5,one,two,four,eight
      dimension fab(2),hess(2)
c
c---------------------------------------------------
c    add fock operator contributions to the diagonal
c    core-core blocks of the hessian
c---------------------------------------------------
      lh=lhcc
c
      if(nco.eq.0) go to 11
c
      do 10 i=1,nco
c
         call mcstv(hess(lh),fab,one,mrs)
         lh=lh+(i+1)*mrs
c
 10   continue
 11   continue
c
      return
      end
