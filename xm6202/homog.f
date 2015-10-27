*deck homog
c***begin prologue     homog
c***date written       930524   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           yukawa, link m6203
c***author             schneider, barry (nsf)
c***source             m6203
c***purpose            calculate radial solutions to laplace equation
c***references         any e&m book
c***routines called
c***end prologue        homog
      subroutine homog(r,j,y,nr,lmax)
      implicit integer (a-z)
      real*8 r, j, y, tmp, tmp1
      dimension r(nr), j(nr,0:lmax), y(nr,0:lmax)
      common/io/ inp,iout
      do 10 pt=1,nr
         tmp=r(pt)
         tmp1=1.d0
         do 20 l=0,lmax
            j(pt,l)=tmp
            y(pt,l)=-tmp1
            tmp=tmp*r(pt)
            tmp1=tmp1/r(pt)
   20    continue     
   10 continue
      return
      end
