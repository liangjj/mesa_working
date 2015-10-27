      subroutine mcfab(fab,nco,mrs,hess,lhcc)
c
c***begin prologue
c***date written       871022   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             %W%   %G%
c
c***purpose
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue
c
c     implicit real*8(a-h,o-p,r-z),             integer*2(q)
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
