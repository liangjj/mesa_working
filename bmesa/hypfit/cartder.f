*deck cartder.f
c***begin prologue     cartder
c***date written       001230   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            derivatives of hyperspherical functions
c***                   wrt cartesian coordinates.
c***description
c***                   exdr1, exdr2 = exact derivatives wrt cartesian
c***                                  variables of a computed hyperspherical
c***                                  function
c***                   exdang, exdr = exact derivatives wrt hyperspherical
c***                                  variables of a computed hyperspherical
c***                                  function
c***                   apdr1, apdr2 = approximate derivatives wrt cartesian
c***                                  variables of dvr fitted function
c***references         
c
c***routines called    
c***end prologue       cartder
      subroutine cartder(apdr1,apdr2,exdang,exdr,exdr1,exdr2,
     1                   rho,ang,c,r1,r2,ply1,ply2,dply1,dply2,
     2                   type,n1,n2)
      implicit integer (a-z)
      real*8 apdr1, apdr2, exdang, exdr, exdr1, exdr2, rho, angrho, ang
      real*8 c, r1, r2, ply1, ply2, dply1, dply2
      character*(*) type
      dimension apdr1(n2,n1), apdr2(n2,n1), exdang(n2,n1), exdr(n2,n1)
      dimension exdr1(n2,n1), exdr2(n2,n1)
      dimension rho(n2,n1), ang(n2,n1), c(n2,n1)
      dimension ply1(n1,n1), dply1(n1,n1), ply2(n2,n2), dply2(n2,n2)
      dimension r1(n1), r2(n2)
      common/io/inp, iout
c
c     get exact derivatives of function wrt (r1,r2)
c     using the chain rule and the exact angular and
c     radial derivatives.
c
      call exdr1r2(exdr1,exdr2,exdang,exdr,rho,ang,type,n1,n2)
c
c     calculate the approximate derivatives wrt(r1,r2) and
c     then compare.
c
      call apdr1r2(apdr1,apdr2,c,ply1,ply2,dply1,dply2,n1,n2)
      call cmpr1r2(exdr1,apdr1,r1,r2,type,'r1',n1,n2)
      call cmpr1r2(exdr2,apdr2,r1,r2,type,'r2',n1,n2)
      return
      end       

















