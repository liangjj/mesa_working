*deck hypfn.f 
c***begin prologue     hypfn
c***date written       960718   (yymmdd)
c***revision date               (yymmdd)
c***keywords           
c***                   hamiltonian, dvr, hyperspherical
c***author             schneider, b. i.(nsf)
c***source             
c***purpose            construction of adiabatic hyperspherical
c***                   harmonics in dvr basis.
c***                   
c
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       hypfn
      subroutine hypfn(pham,pq,pv,typpot,a,b,nword,n,prn)
c
      implicit integer (a-z)
      integer*8 pham, pq, pv
      real*8 a, b, had, hamphi, rho, phi, v, x, y
      character*16 fptoc
      character*80 title
      character*(*) typpot
      logical prn
      dimension pq(2), n(2), a(2,2), b(2,2) 
      common/io/inp, iout      
      pointer (phphi,hamphi(1))
      pointer (prho,rho(1))
      pointer (pphi,phi(1))
      pointer (pv,v(1))
      pointer (phad,had(1))
      pointer (px,x(1))
      pointer (py,y(1))
      pphi=pq(1)
      prho=pq(2)
      phphi=pham
      call memory(wptoin(n(1)),px,ndum,'x',0)
      call memory(wptoin(n(1)),py,ndum,'y',0)
      call memory(wptoin(n(1)*n(2)),pv,nword,'vad',0)
      call rzero(v,n(1)*n(2))
      h=1
      eig=h+n(1)*n(1)*n(2)
      work=eig+n(1)*n(2)
      need=wpadti(work+5*n(1))
      call memory(need,phad,ngot,'had',0) 
      ch=h
      ce=eig
      cr=1
      do 10 i=1,n(2)
         call vadiab(rho(i),phi,x,y,v(cr),typpot,a,b,n(1),2) 
         call setham(had(ch),hamphi,v(cr),rho(i),n(1))
         call dsyev('v','l',n(1),had(ch),n(1),had(ce),had(work),
     1               5*n(1),info)
         call vscale(had(ce),had(ce),1.d0/(rho(i)*rho(i)),n(1))
         title='adiabatic eigenvalues at rho = '//fptoc(rho(i))
         call prntfm(title,had(ce),n(1),1,n(1),1,iout)
         ch=ch+n(1)*n(1)
         ce=ce+n(1)
         cr=cr+n(1)
 10   continue   
      return
      end

