*deck mkthet.f
c***begin prologue     mkthet
c***date written       921223   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           mkthet, link 6201
c***author             schneider, barry (lanl)
c***source             m6050
c***purpose            form cos(theta)
c***references         none
c
c***routines called
c***end prologue       mkthet
      subroutine mkthet (xyz,costhe,npt)
      implicit integer (a-z)
      real*8 xyz, costhe
      dimension xyz(4,npt), costhe(npt)
         do 10 i=1,npt
            costhe(i) = xyz(3,i)/sqrt(xyz(1,i)*xyz(1,i) + 
     1                                xyz(2,i)*xyz(2,i) +
     2                                xyz(3,i)*xyz(3,i))
   10 continue
      return
c
      end
