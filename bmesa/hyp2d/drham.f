*deck drham.f 
c***begin prologue     drham
c***date written       960718   (yymmdd)
c***revision date               (yymmdd)
c***keywords           
c***                   hamiltonian, dvr, hyperspherical
c***author             schneider, b. i.(nsf)
c***source             
c***purpose            explicit construction of total hyperspherical
c***                   hamiltonian in dvr basis.
c***                   
c
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       drham
      subroutine drham(pham,pv,pq,dscale,n,ngot)
c
      implicit integer (a-z)
      integer*8 pham, pv, pq
      real*8 hamrho, hamphi, v, ham, rho, phi, dscale
      dimension phamil(3), pq(2), n(3) 
      common/io/inp, iout      
      pointer (phrho,hamrho(1))
      pointer (phphi,hamphi(1))
      pointer (prho,rho(1))
      pointer (pphi,phi(1))
      pointer (pv,v(1))
      pointer (phtot,ham(1))
      pointer (pind,index(1))
      prho=pq(1)
      pphi=pq(2)
      phrho=pham(1)
      phphi=pham(2)
      need=wptoin(n(3)*n(3))
      call memory(need,phtot,ngot,'ham',0) 
      need=n(1)*n(2)
      call memory(need,pind,nwind,'index',0)
      call ind2d(index,n(1),n(2))
      hrho=1
      hv=hrho+n(1)*n(1)
      he=hv+n(1)
      hs=he+he+n(1)*n(1)
      heg=hs+n(1)
c
c     location of surface amplitude of the two end dvr primitives
c     this is needed to change the bloch operator
c
      srf=heg+n(1) 
      hphi=1
      h=1
      call ham2d(ham(h),v,index,hamrho(hrho),hamphi(hphi),rho,dscale,n)
c
      return
      end

