*deck tdrver.f 
c***begin prologue     tdrver
c***date written       970810   (yymmdd)
c***revision date               (yymmdd)
c***keywords           time dependent
c***                   
c***author             schneider, b. i.(nsf)
c***source             tprop
c***purpose            propagation of three dimensional, time-dependent
c***                   schroedinger equation.
c
c***description        the one, two or three dimensional, time-dependent
c***                   schroedinger equation is solved using a 
c***                   discrete variable representation in space and various
c***                   time-propagation schemes. 
c***references       
c
c***routines called    iosys, util and mdutil
c***end prologue       tdrver
      subroutine tdrver(p1,p2,p3,h01,h02,h03,eig1,eig2,eig3,v,
     1                  psi0,t0,tf,delt,eps,tpert,hbar,dim,nmax,n)
c
      implicit integer (a-z)
      real*8 p1, p2, p3, h01, h02, h03, eig1, eig2, eig3
      real*8 v, psi0, t0, tf, delt, eps 
      real*8 pi, fpkey, omega, efield, width, shift
      real*8 hbar
      real*8 z
      character*(*) tpert
      character*1600 card
      character*80 cpass
      character*8 chrkey, ptype
      logical dollar
      dimension nmax(3)
      dimension p1(nmax(1),nmax(1)), p2(nmax(2),nmax(2))
      dimension p3(nmax(3),nmax(3))
      dimension h01(nmax(1),nmax(1)), h02(nmax(2),nmax(2))
      dimension h03(nmax(3),nmax(3)), eig1(nmax(1)), eig2(nmax(2))
      dimension eig3(nmax(3)), v(n), psi0(n,2)
      dimension eps(2)
      pointer(p,z(1)),(p,iz(1))
      common/io/inp, iout      
      data pi/3.14159265358979323844d0/
      if(tpert.eq.'t') then
         write(iout,1)
      elseif(tpert.eq.'cosine') then
         if(dollar('$v(t)',card,cpass,inp) ) then
            omega=fpkey(card,'electric-field-frequency',10.d0,' ')
            efield=fpkey(card,'electric-field-magnitude',1.d0,' ')
            omega=2.d0*pi*omega   
            write(iout,2) efield, omega
         endif
      elseif(tpert.eq.'gaussian-pulse') then
         if(dollar('$v(t)',card,cpass,inp) ) then
            omega=fpkey(card,'electric-field-frequency',10.d0,' ')
            omega=2.d0*pi*omega   
            efield=fpkey(card,'electric-field-magnitude',1.d0,' ')
            width=fpkey(card,'width-of-gaussian-pulse',0.d0,' ')
            shift=fpkey(card,'shift-of-gaussian-pulse',0.d0,' ')
            write(iout,3) efield, omega, width, shift     
         endif
      endif
      npts=(tf-t0)/delt + 1
      t=1
      psi=t+npts
      vt=psi+2*n*npts
      scr=vt+n
      iscr=wpadti(scr+130+21*2*n)
      need=iscr+51
      call memory(need,p,ngot,'tdrver',0)
      call grid(z(t),t0,delt,npts)
      nstp=npts-1
      call deabm(z(t),z(psi),psi0,t0,tf,delt,info,eps(1),eps(2),
     1           2*n,nstp,h01,h02,h03,eig1,eig2,eig3,v,z(vt),
     2           tpert,efield,omega,width,shift,hbar,z(scr),iz(iscr),
     3           dim,nmax,n)
      call prdiff(p1,p2,p3,eig1,eig2,eig3,z(t),z(psi),dim,
     1            nmax(1),nmax(2),nmax(3),n,npts)
      call memory(-ngot,p,idum,'tdrver',idum)
      return
 1    format(/,5x,'constructing factors for a time-dependent '
     1            'perturbation linear in time')
 2    format(/,5x,'constructing factors for a time-dependent '
     1            'perturbation E*r cosine(omega*t) where:',/,5x,
     2 '                                            E     = ',e15.8,
     3                                                     /,5x,
     4 '                                            omega = ',e15.8)
 3    format(/,5x,'constructing factors for time-dependent '
     1            'perturbation E * r * cos(omega*t) * ',/,5x,
     2            'exp(-width*(t-shift)*(t-shift)) where:',/,5x,
     3 '                                   E     = ',e15.8,/,5x,
     4 '                                   omega = ',e15.8,/,5x,
     5 '                                   width = ',e15.8,/,5x,
     6 '                                   shift = ',e15.8)
      end
