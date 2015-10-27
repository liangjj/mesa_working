*deck @(#)beckef.f	5.2  2/5/95
      subroutine beckef(ngrid,mxgrd,dengrida,dengridb,ga,gb,tmpx,tmpt,
     $     derivs,b,fout,calc,shell)
c***begin prologue     beckef.f
c***date written       930510     (yymmdd)
c***revisiondate       2/5/95
c
c***keywords           xm602, link 602, DFT, Becke, gradient corrected
c
c***author             RUSSO, thomas v.    (lanl)
c***source             @(#)beckef.f	5.2   2/5/95
c***description        computes the value of the Becke functional
c                      given the value of the density on a grid of points.
c***description
c
c***references
c        Becke, A. D., J Chem. Phys. 96(3),2155
c        Johnson, B. G. et al., J. Chem. Phys 98(7),5612
c***routines called
c                      
c
c***end prologue       beckef.f
c
c    ----------
c death to the bloated lackey of imperialist FORTRASH, implicit typing,
c the bane of mine existence and my mortal enemy
c
      implicit none
c
c --- input variables ---
c
c derivs>0 do functonal and derivs
c derivs=0, do only functional.  
c derivs<0 do only derivs
      integer ngrid,derivs,mxgrd,shell
      real*8 b
      character*(*) calc
c
c --- input arrays (unmolested) ---
c
      real*8 dengrida(ngrid),dengridb(ngrid),ga(ngrid),gb(ngrid)
c
c --- input arrays (modified) ---
c
c --- input arrays (scratch) ---
c
      real*8 tmpx(ngrid),tmpt(ngrid,4)
c
c --- output arrays ---
c fout is 5 or 1 depending on value of derivs
c
      real*8 fout(mxgrd,*)
c local
      real*8 two
      parameter (two=2.0d0)
      integer inp,iout
      common/io/inp,iout
c
c zero out the place where the functional goes
c
      if (derivs .ge.0) call rzero(fout,mxgrd)
c 
c if doing derivatives, zero that
c
      if (derivs .ne. 0) then
         call rzero(fout(1,2),mxgrd*4)
      endif
      call becke1(ngrid,mxgrd,dengrida,ga,tmpx,tmpt,fout,2,b,derivs)
      if (calc .eq. 'closed' .or. 
     $     ((calc .eq. 'open'.or.calc.eq.'general') .and. shell.eq.1))
     $     then
         call smul(fout,fout,two,ngrid)
      else if (calc .eq. 'uhf') then
         call becke1(ngrid,mxgrd,dengridb,gb,tmpx,tmpt,fout,4,b,derivs)
      endif

      return
      end

      subroutine becke1(ngrid,mxgrd,den,g,x,t,fout,drstrt,b,derivs)
      implicit none
      integer ngrid,drstrt,derivs,mxgrd
      real*8 den(ngrid),g(ngrid),x(ngrid),fout(mxgrd,*),b,t(ngrid,4)
      real*8 pi,one,half
c
c functions called --- asinh is NOT a standard library routine, need
c to declare it
c
      real*8 asinh
      integer inp,iout
      common/io/inp,iout

      integer i
      real*8 foofac,fr3,un3

      one=1.0d0
      half=.5d0
      pi=4.0d0*atan(one)
      fr3=4.0d0/3.0d0
      un3=1.0d0/3.0d0
      foofac=-1.5d0*(1.0d0/(fr3*pi))**un3

      call rzero(x,ngrid)
      call rzero(t,ngrid)
      do 10 i=1,ngrid
         t(i,4)=sqrt(g(i))
         if (den(i).ne.0) x(i)=t(i,4)/(den(i)**(fr3))
 10   continue
      do 20 i=1,ngrid
         t(i,1)=asinh(x(i))
         t(i,2)=1+6*b*x(i)*t(i,1)
         if (derivs.ne.0) t(i,3)=t(i,2)*t(i,2)
 20   continue 

c
c t=-bx**2/(1+6bx asinh(x))
c
      call vinv(t(1,2),t(1,2),ngrid)
      call vmul(t(1,2),t(1,2),x,ngrid)
      call vmul(t(1,2),t(1,2),x,ngrid)
      call smul(t(1,2),t(1,2),-b,ngrid)
c
c form g(x)=-1.5*(3/(4pi))**(1/3)-bx**2/(1+6bx asinh(x)) on grid
c 
      call sadd(t(1,2),t(1,2),foofac,ngrid)
c
c form becke functional for this piece of the density, add into fout.
c
      if (derivs .ne.0) then
         call vmove(fout(1,drstrt+1),t(1,4),ngrid)
      endif
      do 25 i=1,ngrid
         t(i,4)=den(i)**(fr3)
 25   continue 
      call vmul(t(1,2),t(1,2),t(1,4),ngrid)
c
c only blat into fout if necessary
c
      if ( derivs .ge. 0) then
         call vadd(fout,fout,t(1,2),ngrid)
      endif
c
c calculate stuff for fock matrix now if asked to
c
      if (derivs.ne.0) then
c
c recall that t(1,2) contains rho^(4/3) g(x), and that the first term of
c deriv w.r.t. rho is 4/3 rho^(1/3) g(x).  cheapest way to get that term
c is just 4/3 t(1,2)/rho
c
         do 45 i=1,ngrid
            if (den(i).ne.0) fout(i,drstrt)=t(i,2)/den(i)
 45      continue 
         call smul(fout(1,drstrt),fout(1,drstrt),fr3,ngrid)
c
c now we need to form g'(x)
         call rzero(t(1,2),ngrid)
         do 50 i=1,ngrid
            if (t(i,3).ne.0)
     $       t(i,2)=b*x(i)*(6*b*x(i)*(x(i)/sqrt(x(i)*x(i)+1)-t(i,1))-2)
     $           /t(i,3)
 50      continue 
c 
c finish dB/drho and do dB/dgamma
c note rho^(1/3)*x=sqgma/rho
c
         do 60 i=1,ngrid
            if (den(i).ne.0) fout(i,drstrt)=fout(i,drstrt)-
     $           fr3*(t(i,4)/den(i))*x(i)*t(i,2)
            if (t(i,2).ne.0) fout(i,drstrt+1)=t(i,2)/fout(i,drstrt+1)
 60      continue 
         call smul(fout(1,drstrt+1),fout(1,drstrt+1),half,ngrid)
      endif
      return
      end
