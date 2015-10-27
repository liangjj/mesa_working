*deck triply
      subroutine triply(p,dp,ddp,x,wt,a,b,endpts,nfix,lftbc,rtbc,
     1                  scr,n,npt)
c***begin prologue     triply
c***date written       940504   (yymmdd)
c***revision date               (yymmdd)
c***keywords
c***author             schneider, barry (nsf)
c***source             %W% %G% 
c***purpose            lanczos polynomials.
c***description
c***            
c               
c               
c***references
c
c***routines called
c
c***end prologue       triply
c
      implicit integer (a-z)
      real*8 p, dp, ddp, x, wt, a, b, endpts, scr, dum
      real*8 temp 
      dimension p(npt,n), dp(npt,n), ddp(npt,n), x(npt), wt(npt)
      dimension a(*), b(*), endpts(2), scr(*), temp(2)
      common /io/ inp, iout
c
c     get the points on the interval ( -1 , +1 )
c 
      temp(1) = -1.d0
      temp(2) = 1.d0
      if(nfix.eq.1) then
         temp(1) = 1.d0
      endif
      pleft=0
      pright=0
      if(lftbc.eq.0) then
         pleft=1
      endif
      if(rtbc.eq.0) then
         pright=1
      endif
      call gaussq('legendre',npt,0.d0,0.d0,nfix,temp,scr,x,wt,
     1             dum,.false.)
      call chnvar(x,wt,-1.d0,1.d0,endpts(1),endpts(2),dum,npt)
      call cpoly(p,x,wt,a,b,endpts(1),endpts(2),scr,n,npt,pleft,pright)
      call gpoly(p,dp,ddp,x,a,b,endpts(1),endpts(2),pleft,pright,
     1           n,npt,.false.)
      return
      end















