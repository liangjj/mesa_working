*deck %W% %G%
      subroutine dolsq(ndat,npol,x,f,m1,m2,v1,v2,ipvt)
c
c***begin prologue     %M%
c***date written       940104  (yymmdd)
c***revision date      %G%
c
c***keywords           least squares
c***author             russo, thomas (lanl)
c***source             %W% %G%
c***purpose            solves the least squares problem
c***description
c
c
c***references
c
c***routines called
c***end prologue       %M%
c

      implicit none

      integer ndat,npol,ione,i,j,mdim
      real*8 x(ndat),f(ndat),m1(ndat,npol+1),m2(npol,npol+1),
     $     v1(npol+1),v2(npol+1),rcond,discrim,xmin,fprime
      real*8 emin,fi
      integer ipvt(npol+1)
      integer iout,inp
      common /io/ inp,iout
      ione=1
      mdim=npol+1
c
c form m1(ij)=x(i)**(j-1)
c

      do  1 i=1,ndat
         m1(i,1)=1.0d0
         do 2 j=2,mdim
            m1(i,j)=m1(i,j-1)*x(i)
 2       continue 
 1    continue 

c form m2=m1^T m1 and v1=m1^T f

      call ebtc(m2,m1,m1,mdim,ndat,mdim)
      call ebtc(v1,m1,f,mdim,ndat,ione)

c now solve m2 v2=v1 for v2

      call sgeco(m2,mdim,mdim,ipvt,rcond,v2)
      if ((1.0d0+rcond).eq.1.0) then
         call lnkerr('xm666: Ill conditioned')
      endif

      call sgesl(m2,mdim,mdim,ipvt,v1,0)

      write (iout,*)' The coefficients are:'
      do 3 i=0,npol
         write(iout,*)i,v1(i+1)
 3    continue 

      if (npol.eq.2) then
         write(iout,*)' quadratic has minimum at x=',-v1(2)/(2*v1(3))
      endif
c
c     --- compute the energy and second derivative at the minimum

      if (npol.eq.3) then
         discrim=4*v1(3)**2-12*v1(2)*v1(4)
         if (discrim.gt.0) then
            xmin=(-2*v1(3)+sqrt(discrim))/(6*v1(4))
            write(iout,*)'First minimum of cubic=',xmin
            emin=v1(1)
            do 4 i=1,npol
               fi=float(i)
               emin=emin+v1(i+1)*(xmin**fi)
   4        continue
            write(iout,*) 'energy at minimum=',emin
            fprime=2*v1(3)+6*v1(4)*xmin
            write(iout,*)' Second derivative at minimum 1=',fprime
            xmin=(-2*v1(3)-sqrt(discrim))/(6*v1(4))
            write(iout,*)'second minimum of cubic=',xmin
            fprime=2*v1(3)+6*v1(4)*xmin
            write(iout,*)' Second derivative at minimum 2=',fprime
         endif
      endif
      return
      end
