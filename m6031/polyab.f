*deck polyab.f
c***begin prologue     polyab
c***date written       940120   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           polyab
c***author             schneider, barry (lanl)
c***source             
c***purpose            orthogonal functions on interval (a,b)
c***                   with unit weight
c***references         none
c
c***routines called
c***end prologue       polyab
      subroutine polyab (ply,y,x,a,b,dfct,ddfct,npt,lmax,maxfac)
      implicit integer (a-z)
      real*8 ply, y, x, a, b, dfct, ddfct, f1, norm
      real*8 a1, a2, fac
      dimension ply(npt,0:lmax), y(npt), x(npt)
      dimension  dfct(0:maxfac), ddfct(0:maxfac)
      common /io/ inp, iout
      ldim=lmax+1
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
          lind=m+2
          n1=2
          n2=3
          n3=1
          do 40 i=2,lmax
             ii=i-1
             jj=i-2
             do 50 j=1,npt
                f1=n2*x(j)*ply(j,ii)-n3*ply(j,jj)
                f1=f1/n1
                ply(j,i)=f1
   50        continue
             n1=n1+1
             n2=n2+2
             n3=n3+1
   40     continue
      endif
c     normalize
      do 60 l=0,lmax
         norm=fac*sqrt((l+l+1)/2.d0)
         do 70 j=1,npt
            ply(j,l)=norm*ply(j,l)
   70    continue
   60 continue
      return
c
      end


