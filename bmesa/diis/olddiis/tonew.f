*deck tonew.f
c***begin prologue     tonew
c***date written       961209   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           hartree, non-linear
c***author             schneider, barry (nsf)
c***source             diis
c***purpose            
c***                   
c***                   
c***references         
c
c***routines called    
c***end prologue       tonew
      subroutine tonew(ham,psi,tmat,scr,n)
      implicit integer (a-z)
      real*8 ham, psi, tmat, scr
      dimension ham(n,n), tmat(n,n), psi(n), scr(n,n)
      common/io/inp, iout
      call ebtc(scr,tmat,psi,n,n,1)
      call copy(scr,psi,n)
      call ebc(scr,tmat,ham,n,n,n)
      call ebtc(ham,tmat,scr,n,n,n)
      return
      end       
