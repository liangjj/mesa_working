*deck chnvar
      subroutine chnvar(x,wt,lin,rin,lout,rout,a,npt)
c***begin prologue     chnvar
c***date written       940504   (yymmdd)
c***revision date               (yymmdd)
c***keywords
c***author             schneider, barry (nsf)
c***source             %W% %G% 
c***purpose            lagrange polynomials.
c***description
c***            
c               
c               
c***references
c
c***routines called
c
c***end prologue       chnvar
c
      implicit integer (a-z)
      real*8 x, wt, lin, rin, lout, rout
      real*8 a, b
      dimension x(npt), wt(npt)
      common /io/ inp, iout
c
c     generate actual variables and derivatives with respect to those
c     variables
c
      a = ( lout - rout ) / ( lin - rin )
      b = ( rin*lout - lin*rout ) / ( rin - lin )
      do 10 i=1,npt
         x(i) = a*x(i) + b
         wt(i) = a*wt(i)
 10   continue
      return
      end















