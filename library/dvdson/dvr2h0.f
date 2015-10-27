*deck dvr2h0.f
c***begin prologue     dvr2h0
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           transformation
c***author             schneider, barry (nsf)
c***source             
c***                   
c***references         
c
c***routines called    
c***end prologue       dvr2h0
      subroutine dvr2h0(resid,vec,u01,u02,u03,t1,t2,n1,n2,n3,n,m,dim)
      implicit integer (a-z)
      real*8 resid, vec, u01, u02, u03, t1, t2
      character*80 title
      dimension u01(n1,n1), u02(n2,n2), u03(n3,n3)
      dimension vec(n,m), resid(n,m), t1(*), t2(*)
      common/io/inp, iout
      if(dim.eq.1) then
         call ebtc(vec,u01,resid,n1,n,m)
      elseif(dim.eq.2) then
         call uta2d(resid,vec,u01,u02,t1,n1,n2,m)
      elseif(dim.eq.3) then
         call uta3d(resid,vec,u01,u02,u03,t1,t2,n1,n2,n3,m)
      endif
      return
      end       



