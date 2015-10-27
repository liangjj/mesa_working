*deck nospac.f
c***begin prologue     nospac
c***date written       960723   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           plot, time-dependent, wavefunction
c***author             schneider, barry (nsf)
c***source             
c***purpose            plot or massage the complex time-dependent 
c***                   for analysis wavefunction.
c***                   
c***references         
c
c***routines called    
c***end prologue       nospac
      subroutine nospac(psi,psi0,pt,t,t0,n,m,nt,nc,pnch)
      implicit integer (a-z)
      real*8 psi, psi0, pt, t, t0, fac, q1, q2, q3
      real*8 scr
      complex*16 exact, approx, eye, cfac
      logical pnch
      dimension psi0(m,nc,2), psi(n,nc,2), pt(nt,nt), t(nt), fac(2)
      common/io/inp, iout
      common/punch/pun(4)
      data eye/(0.d0,1.d0)/
      pointer (pscr,scr(1))
      need=wptoin(nt*4)
      call getmem(need,pscr,got,'scr',0)
c
c     we have solved our problem without the inclusion of an overall
c     constant that depends on the initial state of the system at t=t0.
c     we need to include that in our comparison.
c
      write(iout,1) t0
c
c     total wavefunction at t=t0
c
      cfac=psi0(1,1,1)+eye*psi0(1,1,2)
      fac(1)=psi0(1,1,1)
      fac(2)=psi0(1,1,2)
      off=nt
      do 10 i=1,nt
c
c        exact wavefunction for all t
c
c         exact = ( exp(-eye*(t(i)*t(i) - t0*t0)*.5d0) )*cfac
         exact = exp(-.5d0*eye*t(i)*t(i))
         psi(i,1,1) = psi(i,1,1)*pt(i,i) + fac(1)
         psi(i,1,2) = psi(i,1,2)*pt(i,i) + fac(2) 	 
         approx = psi(i,1,1) + eye*psi(i,1,2)
c
c        approximate total wavefunction for all t
c
         write(iout,2) t(i)
         write(iout,3) approx, exact
         scr(i) =             real(exact) 
         scr(i+off) =         real(approx) 
	 scr(i+off+off) =     imag(exact) 
         scr(i+off+off+off) = imag(approx) 
 10   continue
      if(pnch) then
c         call plotit(scr,q1,q2,q2,'plot-1',t,0,0,n,nt,nc,pun(1))
c         call plotit(scr(off+1),q1,q2,q2,'plot-2',t,0,0,n,nt,
c     1               nc,pun(2))
c         call plotit(scr(off+off+1),q1,q2,q2,'plot-3',t,0,0,n,nt,
c     1               nc,pun(3))
c         call plotit(scr(off+off+off+1),q1,q2,q2,'plot-4',t,0,0,n,nt,
c     1               nc,pun(4))
         call plot1(scr,scr(1+off),scr(1+off+off),
     1              scr(1+off+off+off),t,nt,pun(1))
      endif
      call getmem(-got,pscr,dum,'scr',dum)
      return
 1    format(/,1x,'comparison between approximate and exact solution ',
     1       /,1x,'for linear potential in time',/,1x,
     2            'initial time = ',e15.8)
 2    format(/,1x,'time        = ',e15.8)
 3    format(/,10x,'approximate = ',e15.8,1x,e15.8,/,10x, 
     1             'exact       = ',e15.8,1x,e15.8)                   
      end       





