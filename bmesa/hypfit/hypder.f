*deck hypder.f
c***begin prologue     hypder
c***date written       001230   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           
c***author             schneider, barry (nsf)
c***source             
c***purpose            derivatives of hyperspherical functions
c***                   wrt hyperspherical coordinates.
c***                   
c***references         
c
c***routines called    
c***end prologue       hypder
      subroutine hypder(apdr1,apdr2,exdang,exdr,apdang,apdr,rho,ang,c,
     1                  ply1,ply2,dply1,dply2,type,n1,n2)
      implicit integer (a-z)
      real*8 apdr1, apdr2, exdang, exdr, apdanr, apdr, rho, ang
      real*8 c, ply1, ply2, dply1, dply2
      character*(*) type
      dimension apdr1(n2,n1), apdr2(n2,n1)
      dimension exdang(n2,n1), exdr(n2,n1), apdang(n2,n1), apdr(n2,n1)
      dimension rho(n2,n1), ang(n2,n1), c(n2,n1)
      dimension ply1(n1,n1), dply1(n1,n1), ply2(n2,n2), dply2(n2,n2)
      dimension r1(n1), r2(n2)
      common/io/inp, iout
c
c     calculate the approximate derivatives wrt(r1,r2) and
c     then compare.
c
      call apdr1r2(apdr1,apdr2,c,ply1,ply2,dply1,dply2,n1,n2)
      if(type.eq.'Angular') then
         call derang(apdang,apdr1,apdr2,rho,ang,n1,n2)
         call cmphy(exdang,apdang,ang,'Angular',n1,n2)
      elseif(type.eq.'Radial') then
         call derr(apdr,apdr1,apdr2,ang,n1,n2)
         call cmphy(exdr,apdr,rho,'Radial',n1,n2)
      elseif(type.eq.'General') then
         call derang(apdang,apdr1,apdr2,rho,ang,n1,n2)
         call cmphy(exdang,apdang,ang,'General',n1,n2)
         call derr(apdr,apdr1,apdr2,ang,n1,n2)
         call cmphy(exdr,apdr,rho,'General',n1,n2)
      else
         call lnkerr('error in hypder')
      endif
      return
      end       

















