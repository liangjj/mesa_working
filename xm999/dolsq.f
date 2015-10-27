*deck %W% %G%
      subroutine dolsq(ndat,npol,x,y,f,m1,m2,v1,v2,ipvt,terms)
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

      integer ndat,npol,ione,i,j,mdim,terms,ioff,jx,jy
      real*8 x(ndat),y(ndat),f(ndat),m1(ndat,terms),m2(terms,terms),
     $     v1(terms),v2(terms),rcond,discrim,xmin,fprime,determ,temp,
     $     ymin
      real*8 a1,b1,c1,d1,e1,f1
      integer ipvt(terms)
      integer iout,inp
      common /io/ inp,iout
      ione=1
      mdim=terms
c
c form m1(ij) where the elements are the 2-d polynomials
c

      do  1 i=1,ndat
         m1(i,1)=1.0d0
c for each polynomial order
         ioff=2
         do 2 j=1,npol
            do 3 jx=j,0,-1
               jy=j-jx
               m1(i,ioff)=(x(i)**jx)*(y(i)**jy)
               ioff=ioff+1
 3          continue 
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

      write(iout,*)"Constant term:",v1(1)
      ioff=2
c     for each polynomial order
      do 20 j=1,npol
         do 30 jx=j,0,-1
            jy=j-jx
            write(iout,*)"x^",jx,"*y^",jy," term = ",v1(ioff)
            ioff=ioff+1
 30      continue 
 20   continue 
 
      do 50 i=1,ndat
         temp=v1(1)
         ioff=2
c     for each polynomial order
         do 25 j=1,npol
            do 35 jx=j,0,-1
               jy=j-jx
               temp=temp+v1(ioff)*x(i)**jx*y(i)**jy
               ioff=ioff+1
 35         continue 
 25      continue 
         write(iout,*)"point ",i,'function=',f(i),' fit=',temp
 50   continue 

      if (npol.eq.2) then
         a1=v1(1)
         b1=v1(2)
         c1=v1(3)
         d1=v1(5)
         e1=v1(4)
         f1=v1(6)
         determ=4*e1*f1-d1**2
         if (determ .ne. 0) then
            xmin=(-2*f1*b1+d1*c1)/determ
            ymin=(-2*e1*c1+b1*d1)/determ
            write(iout,*)' quadratic surface has minimum at x=',
     $           xmin,' and y=',
     $           ymin
            write(iout,*)'at which point the x deriv is',
     $           b1+d1*ymin+2*e1*xmin,' and the y deriv is ',
     $           c1+d1*xmin+2*f1*ymin
            write(iout,*)'and the function itself is',
     $           a1+b1*xmin+c1*ymin+d1*xmin*ymin+
     $           e1*xmin*xmin+f1*ymin*ymin
         endif
      endif

      return
      end
