*deck ch02dvr.f
c***begin prologue     ch02dvr
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           transformation
c***author             schneider, barry (nsf)
c***source             
c***                   
c***references         
c
c***routines called    
c***end prologue       ch02dvr
      subroutine ch02dvr(resid,vec,u01,u02,u03,t1,t2,n1,n2,n3,n,m,dim)
      implicit integer (a-z)
      complex*16 resid, vec, u01, u02, u03, t1, t2
      character*80 title
      dimension u01(n1,n1), u02(n2,n2), u03(n3,n3)
      dimension vec(n,m), resid(n,m), t1(*), t2(*)
      common/io/inp, iout
      if(dim.eq.1) then
         call cebc(resid,u01,vec,n,n1,m)
      elseif(dim.eq.2) then
c
c         transform residual from h0 to dvr representation
c 
         call cua2d(vec,resid,u01,u02,t1,n1,n2,m)
      elseif(dim.eq.3) then
         call cua3d(vec,resid,u01,u02,u03,t1,t2,n1,n2,n3,m)
      endif
      return
      end       



