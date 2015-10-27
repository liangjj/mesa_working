*deck polyab.f
c***begin prologue     polyab
c***date written       940120   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           polyab
c***author             schneider, barry (lanl)
c***source             
c***purpose            orthogonal functions on interval (a,b)
c***                   with unit weight
c***                   a set of orthogonal polynomials on (a,b)
c***                   and their first and second derivatives
c***                   are computed using the relationship between
c***                   these and the legendre polynomials.  these
c***                   polynomials are not regular at the origin
c***                   if a=0 
c***references         none
c
c***routines called
c***end prologue       polyab
      subroutine polyab (ply,dply,ddply,y,x,a,b,npt,lmax,prnt)
      implicit integer (a-z)
      real*8 ply, dply, ddply, y, x, a, b, norm
      real*8 a1, a2, fac
      logical prnt
      character*80 title
      dimension ply(npt,0:lmax+1), dply(npt,0:lmax), ddply(npt,0:lmax)
      dimension y(npt), x(npt)
      common /io/ inp, iout
      ldim=lmax+2
      call rzero(ply,ldim*npt)
      a1=2.d0/(b-a)
      a2=(a+b)/(b-a)
      fac=sqrt(a1)
      do 10 i=1,npt
         x(i)=a1*y(i)-a2
   10 continue               
      do 20 i=1,npt
         ply(i,0)=1.d0
   20 continue
      if (lmax.gt.0) then
          do 30 i=1,npt
             ply(i,1)=x(i)
   30     continue
      endif   
      if (lmax.gt.1) then
          do 40 l=2,lmax+1
             lm1=l-1
             lm2=l-2
             lplm1=l+l-1
             do 50 j=1,npt
                ply(j,l)=( lplm1*x(j)*ply(j,lm1)-lm1*ply(j,lm2) )/l
   50        continue
   40     continue
      endif
c     derivative
      do 60 j=1,npt
         dply(j,0)=0.d0
         dply(j,1)=1.d0
         ddply(j,0)=0.d0
         ddply(j,1)=0.d0
   60 continue     
      do 70 l=2,lmax
         do 80 j=1,npt
            dply(j,l) = ( l+l-1 )*( ply(j,l-1) + x(j)*dply(j,l-1) )
     1                   - ( l-1 )*dply(j,l-2)
            dply(j,l)=dply(j,l)/l
   80    continue
   70 continue
c     second derivative
      do 90 l=2,lmax
         do 100 j=1,npt
            ddply(j,l) = ( l+l-1)*( 2.d0*dply(j,l-1) +
     1                              x(j)*ddply(j,l-1) )
     2                            - ( l-1 )*ddply(j,l-2)                      
            ddply(j,l) = ddply(j,l)/l
  100    continue
   90 continue
      do 110 l=0,lmax
         norm=sqrt((l+l+1)/2.d0)
         do 120 i=1,npt
            ply(i,l)=ply(i,l)*fac*norm
            dply(i,l)=dply(i,l)*fac*a1*norm
            ddply(i,l)=ddply(i,l)*fac*a1*a1*norm                             
  120    continue
  110 continue
      if (prnt) then
          title='polynomials'
          call prntrm(title,ply,npt,lmax+2,npt,lmax+2,iout)
          title='derivative polynomials'
          call prntrm(title,dply,npt,lmax+1,npt,lmax+1,iout)
          title='second derivative polynomials'
          call prntrm(title,ddply,npt,lmax+1,npt,lmax+1,iout)
      endif
      return
c
      end


